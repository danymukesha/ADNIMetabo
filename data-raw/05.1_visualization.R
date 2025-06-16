library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(RColorBrewer)
library(viridis)

# Based on the previous code, in 05_normalization, variables should be present:
# - valid_data: the filtered data before normalization
# - normalized_data: data after normalization
# - correction_factors: plate-specific correction factors

# ---- 1. Visualization of Correction Factors across Plates ----

plot_correction_factors <- function(correction_factors) {
    if(!all(c("Plate", "Analyte", "CorrectionFactor") %in% colnames(correction_factors))) {
        stop("Correction factors data must contain Plate, Analyte, and CorrectionFactor columns")
    }

    ggplot(correction_factors, aes(x = Plate, y = Analyte, fill = CorrectionFactor)) +
        geom_tile() +
        scale_fill_viridis_c(option = "plasma") +
        labs(title = "Correction Factors by Plate and Analyte",
             subtitle = "Values > 1 indicate plate measurements higher than target",
             x = "Plate",
             y = "Analyte",
             fill = "Correction\nFactor") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ---- 2. Before vs After Normalization Comparison ----

plot_before_after <- function(normalized_data, analytes_to_plot = NULL) {
    required_cols <- c("Analyte", "Concentration", "NormalizedConcentration")
    if(!all(required_cols %in% colnames(normalized_data))) {
        stop("Data must contain Analyte, Concentration, and NormalizedConcentration columns")
    }

    if(is.null(analytes_to_plot)) {
        all_analytes <- unique(normalized_data$Analyte)
        analytes_to_plot <- sample(all_analytes, min(4, length(all_analytes)))
    }

    plot_data <- normalized_data %>%
        filter(Analyte %in% analytes_to_plot)

    ggplot(plot_data, aes(x = Concentration, y = NormalizedConcentration, color = Plate)) +
        geom_point(alpha = 0.7) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
        facet_wrap(~Analyte, scales = "free") +
        labs(title = "Before vs After Normalization",
             subtitle = "Points on the dashed line would be unchanged by normalization",
             x = "Original Concentration",
             y = "Normalized Concentration") +
        theme_bw() +
        theme(legend.position = "bottom")
}

# ---- 3. Distribution of Metabolite Concentrations ----

plot_distribution <- function(valid_data, normalized_data) {
    before_data <- valid_data %>%
        filter(`Sample Type` == "Sample") %>%
        select(RID, VISCODE, Plate, Analyte, Concentration) %>%
        mutate(Type = "Before Normalization")

    after_data <- normalized_data %>%
        filter(`Sample Type` == "Sample") %>%
        select(RID, VISCODE, Plate, Analyte, NormalizedConcentration) %>%
        rename(Concentration = NormalizedConcentration) %>%
        mutate(Type = "After Normalization")

    combined_data <- bind_rows(before_data, after_data)

    top_analytes <- combined_data %>%
        filter(Type == "Before Normalization") %>%
        group_by(Analyte) %>%
        summarize(MedianConc = median(Concentration, na.rm = TRUE)) %>%
        arrange(desc(MedianConc)) %>%
        head(8) %>%
        pull(Analyte)

    combined_data %>%
        filter(Analyte %in% top_analytes) %>%
        ggplot(aes(x = Analyte, y = Concentration, fill = Type)) +
        geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
        scale_fill_brewer(palette = "Set1") +
        labs(title = "Distribution of Metabolite Concentrations",
             subtitle = "Before and after normalization",
             x = "Analyte",
             y = "Concentration") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom")
}

# ---- 4. CV Analysis for QC Samples ----

plot_cv_analysis <- function(valid_data, normalized_data) {
    cv_before <- valid_data %>%
        filter(`Sample Type` == "QC") %>%
        group_by(Analyte) %>%
        summarize(CV = sd(Concentration, na.rm = TRUE) / mean(Concentration, na.rm = TRUE) * 100,
                  Type = "Before Normalization")

    cv_after <- normalized_data %>%
        filter(`Sample Type` == "QC") %>%
        group_by(Analyte) %>%
        summarize(CV = sd(NormalizedConcentration, na.rm = TRUE) /
                      mean(NormalizedConcentration, na.rm = TRUE) * 100,
                  Type = "After Normalization")

    cv_data <- bind_rows(cv_before, cv_after)

    ggplot(cv_data, aes(x = reorder(Analyte, CV), y = CV, fill = Type)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_hline(yintercept = 20, linetype = "dashed", color = "red") +
        scale_fill_brewer(palette = "Set1") +
        labs(title = "Coefficient of Variation (CV) of QC Samples",
             subtitle = "Dashed line represents 20% CV threshold",
             x = "Analyte",
             y = "CV (%)") +
        coord_flip() +
        theme_bw()
}

# ---- 5. Plate Effect Visualization ----

plot_plate_effect <- function(valid_data, normalized_data) {
    top_variable <- valid_data %>%
        filter(`Sample Type` == "Sample") %>%
        group_by(Analyte, Plate) %>%
        summarize(MedianConc = median(Concentration, na.rm = TRUE), .groups = "drop") %>%
        group_by(Analyte) %>%
        summarize(Variability = sd(MedianConc, na.rm = TRUE) / mean(MedianConc, na.rm = TRUE),
                  .groups = "drop") %>%
        arrange(desc(Variability)) %>%
        head(4) %>%
        pull(Analyte)

    before_plate <- valid_data %>%
        filter(`Sample Type` == "Sample", Analyte %in% top_variable) %>%
        mutate(DataType = "Before Normalization")

    after_plate <- normalized_data %>%
        filter(`Sample Type` == "Sample", Analyte %in% top_variable) %>%
        select(-Concentration) %>%  # Remove the existing Concentration column if it is not needed
        rename(Concentration = NormalizedConcentration) %>%
        mutate(DataType = "After Normalization")

    plate_data <- bind_rows(before_plate, after_plate)

    ggplot(plate_data, aes(x = Plate, y = Concentration, fill = Plate)) +
        geom_boxplot(alpha = 0.7) +
        facet_grid(DataType ~ Analyte, scales = "free_y") +
        scale_fill_viridis_d() +
        labs(title = "Plate Effect on Metabolite Concentrations",
             subtitle = "Before and after normalization for top variable analytes",
             x = "Plate",
             y = "Concentration") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none")
}

# ---- 6. Combined Dashboard ----

create_visualization_dashboard <- function(valid_data, normalized_data, correction_factors) {

    p1 <- plot_correction_factors(correction_factors)
    p2 <- plot_before_after(normalized_data)
    p3 <- plot_distribution(valid_data, normalized_data)
    p4 <- plot_cv_analysis(valid_data, normalized_data)
    p5 <- plot_plate_effect(valid_data, normalized_data)

    ((p1 + p2) / (p3 + p4) / p5) +
        plot_layout(heights = c(1, 1, 1.2)) +
        plot_annotation(
            title = "Metabolomics Data Normalization Analysis",
            subtitle = "Visualizing the effects of plate correction normalization",
            caption = "Generated from the provided metabolomics normalization script",
            theme = theme(plot.title = element_text(size = 16, face = "bold"),
                          plot.subtitle = element_text(size = 12),
                          plot.caption = element_text(hjust = 1))
        )
}

# ---- Function to run selected visualizations ----

run_visualization <- function(valid_data, normalized_data, correction_factors, selected_plots = NULL) {
    if (is.null(selected_plots)) {
        selected_plots <- c("correction_factors", "before_after", "distribution",
                            "cv_analysis", "plate_effect", "dashboard")
    }

    plots <- list()

    if ("correction_factors" %in% selected_plots) {
        plots$correction_factors <- plot_correction_factors(correction_factors)
        print(plots$correction_factors)
    }

    if ("before_after" %in% selected_plots) {
        plots$before_after <- plot_before_after(normalized_data)
        print(plots$before_after)
    }

    if ("distribution" %in% selected_plots) {
        plots$distribution <- plot_distribution(valid_data, normalized_data)
        print(plots$distribution)
    }

    if ("cv_analysis" %in% selected_plots) {
        plots$cv_analysis <- plot_cv_analysis(valid_data, normalized_data)
        print(plots$cv_analysis)
    }

    if ("plate_effect" %in% selected_plots) {
        plots$plate_effect <- plot_plate_effect(valid_data, normalized_data)
        print(plots$plate_effect)
    }

    if ("dashboard" %in% selected_plots) {
        plots$dashboard <- create_visualization_dashboard(valid_data, normalized_data, correction_factors)
        print(plots$dashboard)
    }

    # Return the plots invisibly
    invisible(plots)
}

# ---- To run the visualizations, add these lines to the script ----

# To run the visualizations (after your data processing is complete):
# run_visualization(valid_data, normalized_data, correction_factors)

# For specific plots only:
# run_visualization(valid_data, normalized_data, correction_factors,
#                  c("correction_factors", "before_after"))

# You can also save the plots:
# library(ggplot2)
# ggsave("correction_factors.png", plot_correction_factors(correction_factors),
#        width = 10, height = 8, dpi = 300)



# source code trash

cv_before <- transformed_data %>%
    dplyr::filter(`Sample Type` == "Sample") %>%
    dplyr::group_by(Analyte) %>%
    dplyr::summarize(
        CV = sd(Concentration, na.rm = TRUE) / mean(Concentration, na.rm = TRUE) * 100,
        .groups = "drop"
    )

cv_after <- normalized_data %>%
    dplyr::filter(`Sample Type` == "Sample") %>%
    dplyr::group_by(Analyte) %>%
    dplyr::summarize(
        CV = sd(NormalizedConcentration, na.rm = TRUE) /
            mean(NormalizedConcentration, na.rm = TRUE) * 100,
        .groups = "drop"
    )

cv_comparison <- cv_before %>%
    rename(CV_Before = CV) %>%
    left_join(cv_after %>% rename(CV_After = CV), by = "Analyte") %>%
    mutate(Improvement = CV_Before - CV_After)

cv_data <- cv_comparison %>%
    arrange(desc(Improvement)) %>%
    mutate(Analyte = factor(Analyte, levels = Analyte))

#' Plot CV improvement after normalization
ggplot(cv_data, aes(x = Analyte, y = Improvement)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "CV% Improvement After Normalization",
         subtitle = "Positive values indicate improved precision",
         x = "Analyte",
         y = "CV% Reduction")


before_norm <- normalized_data %>%
    filter(Analyte == Analyte, `Sample Type` == "Sample") %>%
    select(Plate, SampleID, Concentration) %>%
    mutate(Type = "Before Normalization")

after_norm <- normalized_data %>%
    filter(Analyte == Analyte, `Sample Type` == "Sample") %>%
    select(Plate, SampleID, NormalizedConcentration) %>%
    rename(Concentration = NormalizedConcentration) %>%
    mutate(Type = "After Normalization")

target <- target_values %>%
    filter(Analyte == Analyte) %>%
    pull(TargetValue)

# Combine data
plot_data <- bind_rows(before_norm, after_norm) %>%
    mutate(Type = factor(Type, levels = c("Before Normalization",
                                          "After Normalization")))

# Create plot
ggplot(plot_data, aes(x = Plate, y = Concentration, color = Type)) +
    geom_point() +
    geom_hline(yintercept = target, linetype = "dashed", color = "black") +
    facet_wrap(~Type, ncol = 1, scales = "free_y") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = paste("QC Sample Normalization for", Analyte),
         subtitle = "Dashed line represents target value",
         x = "Plate",
         y = "Concentration")

reformatted_data <- normalized_data %>%
    # Select the necessary columns for the pivot
    select(Plate, Sampl, SampleType, Analyte, NormalizedConcentration, any_of(id_cols)) %>%
    # Pivot wider to get metabolites back as columns
    pivot_wider(
        id_cols = c(Plate, SampleNumber, SampleType, any_of(id_cols)),
        names_from = Analyte,
        values_from = NormalizedConcentration
    ) %>%
    # Reconstruct the Plate Bar Code column
    mutate(`Plate Bar Code` = paste0(Plate, "-", SampleNumber)) %>%
    # Rename back to original column names
    rename(`Sample Type` = SampleType) %>%
    # Reorder columns to match original format
    select(all_of(orig_cols))

## toreview: reformate the normalized data into the original format `tbbl`.
reformatted_normalized_data <- normalized_data %>%
    # extract `SampleNumber` from SampleID (e.g., "QC_1" -> "1")
    dplyr::mutate(
        SampleNumber = stringr::str_extract(SampleID, "\\d+$"),
        `Plate Bar Code` = paste(Plate, SampleNumber, sep = "-")
    ) |>
    # Rename SampleType to match original column name
    dplyr::rename(`Sample Type` = SampleType) %>%
    # Select and order columns to match original tbbl structure
    dplyr::select(
        `Plate Bar Code`,
        `Sample Type`,
        Analyte,
        NormalizedConcentration
    ) %>%
    # Pivot back to wide format with analytes as columns
    tidyr::pivot_wider(
        names_from = Analyte,
        values_from = NormalizedConcentration
    ) %>%
    # Ensure the order matches the original tbbl (add other columns if needed)
    dplyr::select(`Plate Bar Code`, `Sample Type`, everything())

normalized_data |>
    dplyr::select(RID, VISCODE, Plate, `Sample Type`,
                  Analyte, NormalizedConcentration)


# View the reformatted data
head(reformatted_normalized_data)

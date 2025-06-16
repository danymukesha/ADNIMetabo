library(tidyverse)
library(ggpubr)
library(viridis)

# from 05_normalization.R
data_for_pca <- tbbl |>
    dplyr::select(-c(`Plate Bar Code`, `Sample Type`)) |>
    dplyr::filter(VISCODE != "-") |>
    dplyr::left_join(
        select(adnim, c("RID", "VISCODE", "DX_bl", "TAU")), by = c("RID", "VISCODE")
    ) |>
    dplyr::filter(VISCODE == "bl") |>
    relocate("DX_bl", .after = "VISCODE") |>
    dplyr::filter(DX_bl %in% c("AD", "EMCI", "LMCI", "CN"))

output <- runPCPR2(reformatted_data_all |>
                       dplyr::select(-c(RID, VISCODE)) |> as.matrix(),
                   reformatted_data_all |>
                       dplyr::select(c(RID, VISCODE)) |>
                       dplyr::left_join(
                           select(adnim, c("RID", "VISCODE", "DX_bl", "TAU",
                                           "Male", "AGE", "MMSE")),
                           by = c("RID", "VISCODE")),
                   pct.threshold = 0.8)
output$pR2

plot(output, col = "red",
     main = "Variability in metabolomics data explained by covariates")




# Create a clean theme function similar to the one in your example
clean_theme <- function() {
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title = element_text(face = "bold"),
        text = element_text(size = 12),
        plot.margin = margin(5, 10, 5, 5)
    )
}

# Select only the metabolites for PCA
# Excluding the ID columns (RID, VISCODE, DX_bl)
metabolite_data <- data_for_pca %>%
    select(-c(RID, VISCODE, DX_bl))

metabolite_data <- metabolite_data[, apply(metabolite_data, 2,
                                           function(col) sd(col) != 0)]

# Perform PCA
pca_result <- prcomp(metabolite_data, scale. = TRUE, center = TRUE)

# Extract scores
pca_scores <- as.data.frame(pca_result$x)

# Calculate standard deviations and proportion of variance
s_dev <- pca_result$sdev
var_explained <- s_dev^2 / sum(s_dev^2)

# Create data frame for plotting
plot_data <- data.frame(
    DX_bl = data_for_pca$DX_bl,
    score_PC1 = pca_scores[,1],
    score_PC2 = pca_scores[,2]
)

# Ensure DX_bl is a factor
plot_data$DX_bl <- as.factor(plot_data$DX_bl)

# Get colors for the diagnosis groups - using brighter colors similar to your example
# Define a color palette similar to the image
diagnosis_colors <- c(
    "CN" = "green",    # Bright red/pink
    "EMCI" = "purple",  # Teal
    "LMCI" = "blue",  # Bright yellow
    "AD" = "red"     # Purple
)

# If there are other diagnosis groups not covered, add more colors
all_groups <- levels(plot_data$DX_bl)
if(length(all_groups) > length(diagnosis_colors)) {
    additional_colors <- viridis(length(all_groups) - length(diagnosis_colors))
    names(additional_colors) <- all_groups[!(all_groups %in% names(diagnosis_colors))]
    diagnosis_colors <- c(diagnosis_colors, additional_colors)
}

# Create density plot - no title, clean look
plot_density <- ggplot(data = plot_data, aes(x = score_PC1, color = DX_bl, fill = DX_bl)) +
    geom_density(alpha = 0.25, size = 1.2) +
    theme_minimal() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none"
    ) +
    scale_color_manual(values = diagnosis_colors) +
    scale_fill_manual(values = diagnosis_colors)

# Create scatter plot - with grey reference lines and better formatting
plot_pca <- ggplot(data = plot_data, aes(x = score_PC1, y = score_PC2, color = DX_bl)) +
    geom_hline(yintercept = 0, color = "grey80", size = 0.5) +
    geom_vline(xintercept = 0, color = "grey80", size = 0.5) +
    geom_point(alpha = 0.75, size = 1.5) +
    labs(
        x = paste0("PC1 (", round(100 * var_explained[1], 1), "%)"),
        y = paste0("PC2 (", round(100 * var_explained[2], 1), "%)"),
        color = "Diagnosis"
    ) +
    scale_color_manual(values = diagnosis_colors) +
    theme_minimal() +
    theme(
        panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.line = element_line(color = "black", size = 0.5)
    )

# Combine the plots - density on top, PCA scatter below
combined_plots <- ggarrange(
    plot_density, plot_pca,
    ncol = 1, align = "v", heights = c(1, 5),
    common.legend = TRUE, legend = "right"
)

# Display the combined plots
print(combined_plots)

# Let's also create a loading plot to see which variables contribute most to the PCs
# Get top 20 loadings by absolute value for PC1 and PC2
loadings <- as.data.frame(pca_result$rotation)
loadings$metabolite <- rownames(loadings)

# For PC1
top_pc1_loadings <- loadings %>%
    select(metabolite, PC1) %>%
    arrange(desc(abs(PC1))) %>%
    head(20) %>%
    mutate(metabolite = factor(metabolite, levels = metabolite[order(PC1)]))

# For PC2
top_pc2_loadings <- loadings %>%
    select(metabolite, PC2) %>%
    arrange(desc(abs(PC2))) %>%
    head(20) %>%
    mutate(metabolite = factor(metabolite, levels = metabolite[order(PC2)]))

# Create loading plots
plot_loadings_pc1 <- ggplot(top_pc1_loadings, aes(x = metabolite, y = PC1)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = "Top 20 Metabolites Contributing to PC1", y = "Loading Value", x = "") +
    theme_classic() +
    clean_theme()

plot_loadings_pc2 <- ggplot(top_pc2_loadings, aes(x = metabolite, y = PC2)) +
    geom_bar(stat = "identity", fill = "coral") +
    coord_flip() +
    labs(title = "Top 20 Metabolites Contributing to PC2", y = "Loading Value", x = "") +
    theme_classic() +
    clean_theme()

# Combine the loading plots
loadings_plots <- ggarrange(
    plot_loadings_pc1, plot_loadings_pc2,
    ncol = 2, align = "h",
    common.legend = FALSE
)

# Display the loading plots
print(loadings_plots)

# Create a scree plot to show variance explained by each PC
var_explained_df <- data.frame(
    PC = 1:length(var_explained),
    VarExplained = var_explained,
    CumVarExplained = cumsum(var_explained)
)

scree_plot <- ggplot(var_explained_df[1:20,], aes(x = PC)) +
    geom_bar(aes(y = VarExplained), stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_line(aes(y = CumVarExplained), color = "red", size = 1) +
    geom_point(aes(y = CumVarExplained), color = "red", size = 3) +
    scale_y_continuous(
        name = "Proportion of Variance Explained",
        sec.axis = sec_axis(~., name = "Cumulative Proportion")
    ) +
    labs(title = "Scree Plot", x = "Principal Component") +
    theme_classic() +
    clean_theme()

# Display the scree plot
print(scree_plot)

# For more exploratory analysis, let's also create a biplot with ellipses
biplot_data <- plot_data %>%
    select(DX_bl, score_PC1, score_PC2)

# Calculate ellipse parameters for each group
ellipse_data <- biplot_data %>%
    group_by(DX_bl) %>%
    do({
        data.frame(
            ellipse::ellipse(
                cov(cbind(.$score_PC1, .$score_PC2)),
                centre = c(mean(.$score_PC1), mean(.$score_PC2)),
                level = 0.95
            ),
            group = .$DX_bl[1]
        )
    })

# Add ellipses to the PCA plot
enhanced_pca_plot <- ggplot() +
    geom_point(data = biplot_data, aes(x = score_PC1, y = score_PC2, color = DX_bl), alpha = 0.6, size = 1.5) +
    geom_path(data = ellipse_data, aes(x = x, y = y, color = group), size = 1.2) +
    geom_hline(yintercept = 0, color = "grey") +
    geom_vline(xintercept = 0, color = "grey") +
    labs(
        title = "PCA with 95% Confidence Ellipses",
        x = paste0("PC1 (", round(100 * var_explained[1], 1), "%)"),
        y = paste0("PC2 (", round(100 * var_explained[2], 1), "%)")
    ) +
    scale_color_manual(values = cols) +
    theme_classic() +
    clean_theme()

# Print the enhanced PCA plot
print(enhanced_pca_plot)

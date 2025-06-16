library(tidyverse)
library(ggpubr)
library(viridis)
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(ggsignif) # froo signif plottings

# from 05_normalization.R
adni_serum <- reformatted_data_all |>
    dplyr::left_join(dplyr::select(adnim, c(
        "RID", "VISCODE", "DX_bl", "APOE4", "AGE", "PTGENDER"
        ,"ABETA_bl", "TAU_bl", "MMSE_bl", "MOCA_bl"  # uncomment for descrip step down
    )), by = c("RID", "VISCODE")) |>
    dplyr::filter(VISCODE == "bl") |>
    relocate("DX_bl", .after = "VISCODE") |>
    relocate("AGE", .after = "VISCODE") |>
    relocate("PTGENDER", .after = "VISCODE") |>
    relocate("APOE4", .after = "VISCODE") |>
    dplyr::filter(DX_bl %in% c("AD", "EMCI", "LMCI", "CN")) |>
    dplyr::select(-c(VISCODE)) |>
    dplyr::rename("Gender" = "PTGENDER",
                  "Allgr" = "DX_bl",
                  "Age" = "AGE")

adni_serum_ADvsHC <- adni_serum |>
    dplyr::filter(Allgr %in% c("AD", "CN"))

write.csv(adni_serum_ADvsHC, file = "data-raw/adni_serum_ADvsHC.csv")


# descrip step
adni_serum_ADvsHC |>
    dplyr::select(RID, APOE4, Gender, Age, Allgr, TAU_bl, ABETA_bl) |>
    dplyr::mutate_at(vars(TAU_bl, ABETA_bl), as.numeric) |>
    dplyr::mutate_at(vars(APOE4), as.factor) |>
    DescrTab2::descr(group = "Allgr")

adni_serum_ADvsHC |>
    dplyr::mutate_at(vars(TAU_bl, ABETA_bl), as.numeric) |>
    dplyr::mutate_at(vars(APOE4), as.factor) |>
    ggplot2::ggplot(mapping = aes(x = APOE4, y = TAU_bl)) +
    geom_boxplot() +
    ggsignif::geom_signif(
        comparisons = list(c("0", "1"), c("0", "2"), c("1", "2")),
        map_signif_level = TRUE,
        textsize = 6,
        margin_top = 0.08,
        step_increase = 0.05,
        tip_length = 0.01,
    )

adni_serum_ADvsHC |>
    dplyr::mutate(TAU_bl = as.numeric(TAU_bl),
        ABETA_bl = as.numeric(ABETA_bl),
        APOE4 = as.factor(APOE4),
        MMSE_bl = as.numeric(MMSE_bl),
        MOCA_bl = as.numeric(MOCA_bl)) |>
    ggplot(aes(x = APOE4, y = MOCA_bl)) +
    geom_boxplot(fill = "coral3", color = "black", outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, color = "black", size = 1) +
    ggsignif::geom_signif(
        comparisons = list(c("0", "1"), c("0", "2"), c("1", "2")),
        map_signif_level = TRUE,
        textsize = 4.5,
        margin_top = 0.05,
        step_increase = 0.07,
        tip_length = 0.01) +
    labs(x = "APOE4 Genotype", y = "Baseline MOCA (MOCA_bl)",
        title = "Association between APOE4 Status and MOCA Levels") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey90"))


adni_serum_ADvsHC |>
    dplyr::select(MMSE_bl, MOCA_bl, Allgr) |>
    pivot_longer(cols = c(MMSE_bl, MOCA_bl), names_to = "Test",
                 values_to = "Score") |>
    filter(!is.na(Score)) |>
    ggplot(aes(x = Score, fill = Allgr)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~Test, scales = "free") +
    labs(title = "Cognitive test score distributions by diagnosis",
         x = "Score", y = "Density", fill = "Diagnosis") +
    theme_base() +
    theme(strip.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "top")

# quick check on apoe4 homozygous
adnim |> dplyr::filter(APOE4 == 2) |> dplyr::filter(DX_bl =="CN") |> dplyr::filter(RID == "520")

plot_data <- adnim %>%
    filter(APOE4 == 2, DX_bl == "CN") %>%
    mutate(
        DX = ifelse(DX == "", "Unknown", DX),
        VISCODE = factor(VISCODE, levels = paste0("m", seq(0, 96, by = 6)))
    ) |>
    dplyr::select(RID, PTID, VISCODE, DX_bl, DX, APOE4, AGE, PTGENDER
                  ,"ABETA_bl", "TAU_bl", "MMSE_bl", "MOCA_bl"
                  ) |> filter(if_all(VISCODE, ~ !is.na(.))) #|>
    #filter(if_all(everything(), ~ !is.na(.)))

baseline_labels <- plot_data %>%
    group_by(PTID) %>%
    slice_min(order_by = VISCODE, with_ties = FALSE) %>%
    mutate(
        biomarker_label = paste0(
            "Age: ", AGE, "\n"
            #,"ABETA: ", ABETA_bl, "\n", "TAU: ", TAU_bl, "\n"
            ,"MMSE: ", MMSE_bl, ", ", "MOCA: ", MOCA_bl
        )
    )

ggplot(plot_data, aes(x = VISCODE, y = PTID, group = PTID)) +
    geom_line(aes(color = DX), linewidth = 1.2) +
    geom_point(aes(shape = PTGENDER, color = DX), size = 3) +
    geom_text(
        data = baseline_labels,
        aes(label = biomarker_label),
        hjust = 0, vjust = -0.3,
        size = 3,
        color = "black"
    ) +
    #geom_text(
    #    data = plot_data %>%group_by(PTID) %>%
    #        slice_min(order_by = VISCODE, with_ties = FALSE),  # label only first visit
    #    aes(label = paste0("Age: ", AGE)),
    #    hjust = 0.2, vjust = -0.9, size = 3, color = "black"
    #) +
    scale_color_manual(values = c("CN" = "forestgreen", "MCI" = "orange", "Unknown" = "gray")) +
    scale_shape_manual(values = c("Male" = 16, "Female" = 17)) +
    labs(
        title = "Conversion from CN to MCI/AD/Unknown among APOE4 homozygous individuals",
        subtitle = "Baseline DX = CN | Tracking diagnostic change over time",
        x = "Months since baseline (VISCODE)",
        y = "Participant ID (PTID)",
        color = "Diagnosis",
        shape = "Gender",
        linetype = "Gender"
    ) +
    theme_classic() -> cn_apoe4_plot

ggsave("pictures/cn_apoe4_plot.tiff", cn_apoe4_plot,
       width = 15, height = 9, device = 'tiff', dpi = 1200)

long_time_data <- plot_data |>
    dplyr::mutate(TAU_bl = as.numeric(TAU_bl), ABETA_bl = as.numeric(ABETA_bl)) |>
    pivot_longer(cols = c(ABETA_bl, TAU_bl, MMSE_bl, MOCA_bl),
                 names_to = "Biomarker", values_to = "Value") |>
    mutate(Value = as.numeric(Value))

# Plot each biomarker over time per PTID
ggplot(long_time_data, aes(x = VISCODE, y = Value, group = PTID)) +
    geom_line(aes(color = DX), alpha = 0.7, linewidth = 1.1) +
    geom_point(aes(shape = PTGENDER, color = DX), size = 2.2) +
    facet_wrap(~ Biomarker, scales = "free_y") +
    scale_color_manual(values = c("CN" = "forestgreen",
                                  "MCI" = "orange",
                                  "Unknown" = "gray")) +
    scale_shape_manual(values = c("Male" = 16, "Female" = 17)) +
    labs(
        title = "Trajectory of Gold-standard Biomarkers Over Time",
        subtitle = "APOE4 homozygous individuals (Baseline CN)",
        x = "Months since baseline (VISCODE)",
        y = "Biomarker Value",
        color = "Diagnosis",
        shape = "Gender"
    ) +
    theme_bw()


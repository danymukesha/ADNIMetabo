# install.packages(c("tidyverse", "ggradar", "scales"))
library(tidyverse)
library(ggradar)
library(scales)

# ---- PREPARE YOUR TABLE ----
df <- summary_table_print %>%
    rownames_to_column("Model_ID") %>%
    separate(`95% Confidence Interval`,
        into = c("CI_low", "CI_high"),
        sep = ", ", remove = FALSE
    ) %>%
    mutate(
        CI_low = as.numeric(gsub("[()]", "", CI_low)),
        CI_high = as.numeric(gsub("[()]", "", CI_high))
    ) %>%
    mutate(AUC_color = AUC) %>%
    select(
        Model_ID, Accuracy, Sensitivity, Specificity, `F1 Score`,
        `Balanced Accuracy`, AUC, CI_low, CI_high
    )

# ---- NORMALIZE FOR RADAR ----
df_scaled <- df %>%
    mutate(across(
        c(
            Accuracy, Sensitivity, Specificity, `F1 Score`,
            `Balanced Accuracy`, AUC
        ),
        rescale
    )) # scales 0–1 for radar

# ggradar requires specific formatting
df_radar <- df_scaled %>%
    select(
        Model_ID, Accuracy, Sensitivity, Specificity,
        `F1 Score`, `Balanced Accuracy`, AUC
    )

# ---- RADAR PLOT ----
ggradar(df_radar,
    group.colours = scales::col_numeric(
        palette = c("navy", "cyan", "yellow", "orange", "red"),
        domain = df$AUC
    )(df$AUC),
    grid.min = 0, grid.mid = 0.5, grid.max = 1,
    axis.label.size = 4,
    group.line.width = 1.2,
    group.point.size = 2.5,
    background.circle.colour = "grey95",
    gridline.mid.colour = "grey80"
) +
    ggtitle("Model Performance Radar with AUC Heat Coloring\nAccuracy CI Shaded Separately") +
    theme(plot.title = element_text(size = 16, face = "bold"))



data0 |>
    left_join(ADNIMetabo::demographic_data |> unique() |>
        select(RID, PTMARRY), by = "RID") |>
    left_join(ADNIMetabo::baseline_clinical_data |> unique() |>
        select(RID, MMSE_bl, CDRSB_bl), by = "RID") |>
    dplyr::select(Gender, Age, Allgr, CDRSB_bl, PTMARRY) |>
    tibble() |>
    DescrTab2::descr(group = "Allgr")

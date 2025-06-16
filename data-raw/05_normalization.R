library(dplyr)
library(tidyr)
library(tidyverse)

# for test:
# Sample concentration data format (example)
# data <- data.frame(
#   Plate = c("P1", "P1", "P2", "P2", ...),
#   SampleID = c("QC2_1", "Sample1", "QC2_2", "Sample2", ...),
#   SampleType = c("QC", "Unknown", "QC", "Unknown", ...),
#   Analyte = c("MetA", "MetA", "MetA", "MetA", ...),
#   Concentration = c(1.2, 0.8, 1.5, 0.7, ...)
# )

# Target value table format (eg.for target normalization)
# target_values <- data.frame(
#   Analyte = c("MetA", "MetB", ...),
#   QC_Level = "QC2",
#   TargetValue = c(1.3, 2.1, ...)
# )

# FIA ====
## STEP 1: Preprocessing ====
# valid_data <- trans |>
#    dplyr::filter(Concentration > LOD_threshold, !is.na(Concentration))

tbbl <- ADNIMetabo::FIA_processed$qc_w_metabo_data
# todo: now FIA for test, then will create function
# to automatized the process with but FIA and UPLC

# check out the duplicated RID - VISCODE
tbbl |> dplyr::filter(`Sample Type` != "QC") |>
    dplyr::group_by(RID, VISCODE, `Plate Bar Code`) |>
    dplyr::filter(n() > 1 ) |>
    dplyr::summarize(n = n(),  .groups = "keep") -> duplicated_tbbl

duplicated_rows <- duplicated_tbbl$RID
tbbl <- tbbl |> dplyr::filter(!RID %in% c(duplicated_rows))
tbbl

orig_cols <- names(tbbl)

transformed_data <- tbbl |>
    tidyr::separate(`Plate Bar Code`,
                    into = c("Plate", "PlateID"),
                    sep = "-") |>
    dplyr::select(-PlateID) |>
    dplyr::mutate(SampleID = paste0(RID, "_", VISCODE)) |>
    tidyr::pivot_longer(
        cols = -c(RID, VISCODE, Plate, `Sample Type`, SampleID),
        names_to = "Analyte",
        values_to = "Concentration"
    )

id_cols <- names(transformed_data)

target_values <- tbbl |>
    dplyr::filter(`Sample Type` == "QC") |>
    dplyr::select(-RID, -VISCODE, -`Plate Bar Code`, -`Sample Type`) |>
    summarise(across(everything(), mean, na.rm = TRUE)) |>

    tidyr::pivot_longer(
        cols = everything(),
        names_to = "Analyte",
        values_to = "TargetValue"
    ) |>

    dplyr::mutate(QC_Level = "QC") |>
    dplyr::select(Analyte, QC_Level, TargetValue)

head(target_values)
head(transformed_data)

valid_data <- transformed_data # apply assertion if needed for TODO.

## CONFIGURATION ====
LOD_threshold <- 0.0001
use_median <- TRUE
normalization_type <- "target" # or "reference"
reference_sample <- "QC" # can be a sample name or sample type

## STEP 2: Calculate correction factors ====
if (normalization_type == "target") {

    correction_factors <- valid_data |>
        dplyr::filter(`Sample Type` == "QC", grepl(reference_sample,
                                                   `Sample Type`)) |>
        dplyr::group_by(Plate, Analyte) |>
        dplyr::summarize(
            Stat = if (use_median)
                median(Concentration)
            else
                mean(Concentration),
            .groups = "drop"
        ) |>
        dplyr::left_join(target_values, by = "Analyte") |>
        dplyr::mutate(CorrectionFactor = Stat / TargetValue)

} else if (normalization_type == "reference") {

    # step A: calculate A = mean/median per plate, per analyte/metabolite
    A <- valid_data |>
        dplyr::filter(SampleID == reference_sample | SampleType == reference_sample) |>
        dplyr::group_by(Plate, Analyte) |>
        dplyr::summarize(A = if (use_median)
            median(Concentration)
            else
                mean(Concentration),
            .groups = "drop")

    # step B: calculate B = overall mean/median per analyte/metabolite
    B <- A |>
        dplyr::group_by(Analyte) |>
        dplyr::summarize(B = if (use_median)
            median(A)
            else
                mean(A), .groups = "drop")

    # step C: correction factor = A / B
    correction_factors <- A |>
        dplyr::left_join(B, by = "Analyte") |>
        dplyr::mutate(CorrectionFactor = A / B)
}

## STEP 3: Normalize concentrations ====
normalized_data <- valid_data |>
    dplyr::left_join(correction_factors, by = c("Plate", "Analyte")) |>
    dplyr::mutate(NormalizedConcentration = Concentration / CorrectionFactor) |>
    dplyr::select(-c(Stat, QC_Level, TargetValue))
    # |>  filter(if_any(everything(), ~ is.na(.)))

# for Samples without their correspoonding plate for QCs
correction_factors_mean <- correction_factors |>
    dplyr::select(-c(Plate)) |>
    dplyr::group_by(Analyte) |>
    dplyr::summarize(
        mean = mean(CorrectionFactor, na.rm = TRUE),
        .groups = "drop"
    ) |>
    dplyr::rename(CorrectionFactor = mean)

normalized_data_mean <- normalized_data |> filter(if_any(everything(), ~ is.na(.))) |>
    dplyr::select(where(function(x) !all(is.na(x)))) |> unique() |>
    dplyr::left_join(correction_factors_mean, by = c("Analyte")) |>
    dplyr::mutate(NormalizedConcentration = Concentration / CorrectionFactor)

normalized_data <- normalized_data |> filter(if_all(everything(), ~ !is.na(.))) |>
    rbind(normalized_data_mean) |>
    dplyr::filter(`Sample Type` == "Sample") |> unique()

## Result ====
head(normalized_data)

reformatted_data <- normalized_data |>
    dplyr::select(RID, VISCODE, Plate, `Sample Type`, Analyte, NormalizedConcentration, Plate) |>
    tidyr::pivot_wider(
        id_cols = c(RID, VISCODE, Plate),
        names_from = Analyte,
        values_from = NormalizedConcentration,
        values_fn = ~ median(.x, na.rm = TRUE)
    )

reformatted_data_fia <- reformatted_data

# UPLC ====
## STEP 1: Preprocessing ====
# valid_data <- trans |>
#    dplyr::filter(Concentration > LOD_threshold, !is.na(Concentration))

tbbl <- UPLC_processed$qc_w_metabo_data
# todo: now FIA for test, then will create function
# to automatized the process with but FIA and UPLC

# check out the duplicated RID - VISCODE
tbbl |> dplyr::filter(`Sample Type` != "QC") |>
    dplyr::group_by(RID, VISCODE, `Plate Bar Code`) |>
    dplyr::filter(n() > 1 ) |>
    dplyr::summarize(n = n(),  .groups = "keep") -> duplicated_tbbl

duplicated_rows <- duplicated_tbbl$RID
tbbl <- tbbl |> dplyr::filter(!RID %in% c(duplicated_rows))
tbbl

orig_cols <- names(tbbl)

transformed_data <- tbbl |>
    tidyr::separate(`Plate Bar Code`,
                    into = c("Plate", "PlateID"),
                    sep = "-") |>
    dplyr::select(-PlateID) |>
    dplyr::mutate(SampleID = paste0(RID, "_", VISCODE)) |>
    tidyr::pivot_longer(
        cols = -c(RID, VISCODE, Plate, `Sample Type`, SampleID),
        names_to = "Analyte",
        values_to = "Concentration"
    )

id_cols <- names(transformed_data)

target_values <- tbbl |>
    dplyr::filter(`Sample Type` == "QC") |>
    dplyr::select(-RID, -VISCODE, -`Plate Bar Code`, -`Sample Type`) |>
    summarise(across(everything(), mean, na.rm = TRUE)) |>

    tidyr::pivot_longer(
        cols = everything(),
        names_to = "Analyte",
        values_to = "TargetValue"
    ) |>

    dplyr::mutate(QC_Level = "QC") |>
    dplyr::select(Analyte, QC_Level, TargetValue)

head(target_values)
head(transformed_data)

valid_data <- transformed_data # apply assertion if needed for TODO.

## CONFIGURATION ====
LOD_threshold <- 0.0001
use_median <- TRUE
normalization_type <- "target" # or "reference"
reference_sample <- "QC" # can be a sample name or sample type

## STEP 2: Calculate correction factors ====
if (normalization_type == "target") {

    correction_factors <- valid_data |>
        dplyr::filter(`Sample Type` == "QC", grepl(reference_sample,
                                                   `Sample Type`)) |>
        dplyr::group_by(Plate, Analyte) |>
        dplyr::summarize(
            Stat = if (use_median)
                median(Concentration)
            else
                mean(Concentration),
            .groups = "drop"
        ) |>
        dplyr::left_join(target_values, by = "Analyte") |>
        dplyr::mutate(CorrectionFactor = Stat / TargetValue)

} else if (normalization_type == "reference") {

    # step A: calculate A = mean/median per plate, per analyte/metabolite
    A <- valid_data |>
        dplyr::filter(SampleID == reference_sample | SampleType == reference_sample) |>
        dplyr::group_by(Plate, Analyte) |>
        dplyr::summarize(A = if (use_median)
            median(Concentration)
            else
                mean(Concentration),
            .groups = "drop")

    # step B: calculate B = overall mean/median per analyte/metabolite
    B <- A |>
        dplyr::group_by(Analyte) |>
        dplyr::summarize(B = if (use_median)
            median(A)
            else
                mean(A), .groups = "drop")

    # step C: correction factor = A / B
    correction_factors <- A |>
        dplyr::left_join(B, by = "Analyte") |>
        dplyr::mutate(CorrectionFactor = A / B)
}

## STEP 3: Normalize concentrations ====
normalized_data <- valid_data |>
    dplyr::left_join(correction_factors, by = c("Plate", "Analyte")) |>
    dplyr::mutate(NormalizedConcentration = Concentration / CorrectionFactor) |>
    dplyr::select(-c(Stat, QC_Level, TargetValue))
# |>  filter(if_any(everything(), ~ is.na(.)))

# for Samples without their correspoonding plate for QCs
correction_factors_mean <- correction_factors |>
    dplyr::select(-c(Plate)) |>
    dplyr::group_by(Analyte) |>
    dplyr::summarize(
        mean = mean(CorrectionFactor, na.rm = TRUE),
        .groups = "drop"
    ) |>
    dplyr::rename(CorrectionFactor = mean)

normalized_data_mean <- normalized_data |> filter(if_any(everything(), ~ is.na(.))) |>
    dplyr::select(where(function(x) !all(is.na(x)))) |> unique() |>
    dplyr::left_join(correction_factors_mean, by = c("Analyte")) |>
    dplyr::mutate(NormalizedConcentration = Concentration / CorrectionFactor)

normalized_data <- normalized_data |> filter(if_all(everything(), ~ !is.na(.))) |>
    rbind(normalized_data_mean) |>
    dplyr::filter(`Sample Type` == "Sample") |> unique()

## Result ====
head(normalized_data)

reformatted_data <- normalized_data |>
    dplyr::select(RID, VISCODE, Plate, `Sample Type`, Analyte, NormalizedConcentration, Plate) |>
    tidyr::pivot_wider(
        id_cols = c(RID, VISCODE, Plate),
        names_from = Analyte,
        values_from = NormalizedConcentration,
        values_fn = ~ median(.x, na.rm = TRUE)
    )

reformatted_data_uplc <- reformatted_data

reformatted_data_all <- reformatted_data_fia |>
    dplyr::left_join(reformatted_data_uplc, by = c("RID", "VISCODE")) |>
    dplyr::mutate(Plate = paste0(Plate.x, "_", Plate.y)) |>
    dplyr::relocate("Plate", .after = "VISCODE") |>
    dplyr::select(-c(Plate.x, Plate.y))

# Save the data sets ====
usethis::use_data(reformatted_data_all, overwrite = TRUE)

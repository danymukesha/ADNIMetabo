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

tbl <- FIA_processed$qc_w_metabo_data
# todo: now FIA for test, then will create function
# to automatized the process with but FIA and UPLC

data <- tbl |>
    tidyr::separate(`Plate Bar Code`, into = c("Plate", "SampleNumber"), sep = "-") %>%
    dplyr::rename(SampleType = `Sample Type`) %>%
    mutate(SampleID = paste0(SampleType, "_", SampleNumber)) %>%

    pivot_longer(
        cols = -c(RID, VISCODE, Plate, SampleType, SampleNumber, SampleID),
        names_to = "Analyte",
        values_to = "Concentration"
    ) %>%
    select(Plate, SampleID, SampleType, Analyte, Concentration)


target_values <- tbl %>%
    filter(`Sample Type` == "QC") %>%
    select(-RID, -VISCODE, -`Plate Bar Code`, -`Sample Type`) %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%

    pivot_longer(
        cols = everything(),
        names_to = "Analyte",
        values_to = "TargetValue"
    ) %>%

    mutate(QC_Level = "QC") %>%
    select(Analyte, QC_Level, TargetValue)

head(target_values)
head(transformed_data)

# CONFIGURATION ====
LOD_threshold <- 0.0001
use_median <- TRUE
normalization_type <- "target" # or "reference"
reference_sample <- "QC" # can be a sample name or sample type

# STEP 1: Preprocessing ====
valid_data <- data %>%
    dplyr::filter(Concentration > LOD_threshold, !is.na(Concentration))

# STEP 2: Calculate correction factors ====
if (normalization_type == "target") {

    correction_factors <- valid_data %>%
        filter(SampleType == "QC", grepl(reference_sample, SampleID)) %>%
        group_by(Plate, Analyte) %>%
        summarize(Stat = if (use_median) median(Concentration) else mean(Concentration), .groups = "drop") %>%
        left_join(target_values, by = "Analyte") %>%
        mutate(CorrectionFactor = Stat / TargetValue)

} else if (normalization_type == "reference") {

    # step A: calculate A = mean/median per plate, per analyte
    A <- valid_data %>%
        filter(SampleID == reference_sample | SampleType == reference_sample) %>%
        group_by(Plate, Analyte) %>%
        summarize(A = if (use_median) median(Concentration) else mean(Concentration), .groups = "drop")

    # step B: calculate B = overall mean/median per analyte
    B <- A %>%
        group_by(Analyte) %>%
        summarize(B = if (use_median) median(A) else mean(A), .groups = "drop")

    # step C: correction factor = A / B
    correction_factors <- A %>%
        left_join(B, by = "Analyte") %>%
        mutate(CorrectionFactor = A / B)
}

# STEP 3: Normalize concentrations ====
normalized_data <- valid_data %>%
    left_join(correction_factors, by = c("Plate", "Analyte")) %>%
    mutate(NormalizedConcentration = Concentration / CorrectionFactor)

# result ====
head(normalized_data)

## todo: reformate the normalized data into the original format `tbl`.

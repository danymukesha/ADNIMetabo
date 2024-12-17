library(data.table)
library(dplyr)
library(tidyr)

adnimerge <-  data.table::fread(
    input = "data-raw/ADNIMERGE_20Sep2024.csv",
    sep = ",")

adnimerge_dict <-  data.table::fread(
    input = "data-raw/ADNIMERGE_DICT_20Sep2024.csv",
    sep = ",") |> tibble::tibble()

transpose_df <- function(df) {
    t_df <- data.table::transpose(df)
    colnames(t_df) <- rownames(df)
    rownames(t_df) <- colnames(df)
    t_df <- t_df %>%
        tibble::rownames_to_column(.data = .) %>%
        tibble::as_tibble(.)
    return(t_df)
}

adnimerge_dict |> dplyr::select(FLDNAME, TEXT) %>%
    `row.names<-`(., NULL) %>%
    tibble::column_to_rownames(var = "FLDNAME") |>
    transpose_df() |>
    dplyr::rename(FLDNAME = rowname) ->  adnimerge_dict

# Demographic data ====
demographic_data <- adnimerge %>%
    select(RID, AGE, PTGENDER, PTEDUCAT, PTETHCAT, PTRACCAT, PTMARRY)


# Baseline clinical data ====
baseline_clinical_data <- adnimerge %>%
    select(RID, DX_bl, CDRSB_bl, ADAS11_bl, ADAS13_bl, MMSE_bl,
           RAVLT_immediate_bl, RAVLT_learning_bl, RAVLT_forgetting_bl,
           FAQ_bl, MOCA_bl)

# Follow-up clinical data ====
follow_up_clinical_data <- adnimerge %>%
    select(RID, VISCODE, EXAMDATE, DX, CDRSB, ADAS11, ADAS13, MMSE,
           RAVLT_immediate, RAVLT_learning, RAVLT_forgetting)

# Biomarker data ====
biomarker_data <- adnimerge %>%
    select(RID, ABETA, TAU, PTAU, FDG, PIB, AV45)

# Imaging data ====
imaging_data <- adnimerge %>%
    select(RID, FLDSTRENG, FSVERSION, Ventricles, Hippocampus, WholeBrain,
           Entorhinal, Fusiform, MidTemp, ICV)

# Cognitive composite scores ====
cognitive_composite_data <- adnimerge %>%
    select(RID, mPACCdigit, mPACCtrailsB)

# ECog  scores ====
ecog_scores_data <- adnimerge %>%
    select(RID, EcogPtMem, EcogPtLang, EcogPtVisspat, EcogPtPlan, EcogPtOrgan,
           EcogPtDivatt, EcogPtTotal, EcogSPMem, EcogSPLang, EcogSPVisspat,
           EcogSPPlan, EcogSPOrgan, EcogSPDivatt, EcogSPTotal)

# Save the data sets ====
usethis::use_data(demographic_data, overwrite = TRUE)
usethis::use_data(baseline_clinical_data, overwrite = TRUE)
usethis::use_data(follow_up_clinical_data, overwrite = TRUE)
usethis::use_data(biomarker_data, overwrite = TRUE)
usethis::use_data(imaging_data, overwrite = TRUE)
usethis::use_data(cognitive_composite_data, overwrite = TRUE)
usethis::use_data(ecog_scores_data, overwrite = TRUE)

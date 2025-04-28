library(data.table)
library(dplyr)
library(tidyr)
library(tidyverse)
library(tibble)

adnimerge <-  data.table::fread(
    input = "data-raw/ADNIMERGE_20Sep2024.csv",
    sep = ",") |> tibble::tibble()

visits <-  data.table::fread(
    input = "data-raw/VISITS_19Dec2024.csv",
    sep = ",") |> tibble::tibble()

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
    dplyr::select(RID, AGE, PTGENDER, PTEDUCAT, PTETHCAT, PTRACCAT, PTMARRY)


# Baseline clinical data ====
baseline_clinical_data <- adnimerge %>%
    dplyr::select(RID, DX_bl, CDRSB_bl, ADAS11_bl, ADAS13_bl, MMSE_bl,
                  RAVLT_immediate_bl, RAVLT_learning_bl, RAVLT_forgetting_bl,
                  FAQ_bl, MOCA_bl)

# Follow-up clinical data ====
follow_up_clinical_data <- adnimerge %>%
    dplyr::select(RID, VISCODE, EXAMDATE, DX, CDRSB, ADAS11, ADAS13, MMSE,
                  RAVLT_immediate, RAVLT_learning, RAVLT_forgetting)

# Biomarker data ====
biomarker_data <- adnimerge %>%
    dplyr::select(RID, ABETA, TAU, PTAU, FDG, PIB, AV45)

# Imaging data ====
imaging_data <- adnimerge %>%
    dplyr::select(RID, FLDSTRENG, FSVERSION, Ventricles, Hippocampus,
                  WholeBrain, Entorhinal, Fusiform, MidTemp, ICV)

# Cognitive composite scores ====
cognitive_composite_data <- adnimerge %>%
    dplyr::select(RID, mPACCdigit, mPACCtrailsB)

# ECog  scores ====
ecog_scores_data <- adnimerge %>%
    dplyr::select(RID, EcogPtMem, EcogPtLang, EcogPtVisspat, EcogPtPlan,
                  EcogPtOrgan, EcogPtDivatt, EcogPtTotal, EcogSPMem,
                  EcogSPLang, EcogSPVisspat, EcogSPPlan, EcogSPOrgan,
                  EcogSPDivatt, EcogSPTotal)

# link RID to PTID
RID_to_PTID <- adnimerge %>%
    dplyr::select(RID, PTID)

# re-merge all the datasets
ADNIMERGE <- new(
    "ADNIMERGE",
    complete_metadata = adnimerge,
    demographic_data = demographic_data,
    baseline_clinical_data = baseline_clinical_data,
    follow_up_clinical_data = follow_up_clinical_data,
    biomarker_data = biomarker_data,
    imaging_data = imaging_data,
    cognitive_composite_data = cognitive_composite_data,
    ecog_scores_data = ecog_scores_data,
    RID_to_PTID = RID_to_PTID,
    logs = character(0),  # optional logs if needed
    misc_data = list()     # optional miscellaneous data
)

# Save the data sets ====
usethis::use_data(demographic_data, overwrite = TRUE)
usethis::use_data(baseline_clinical_data, overwrite = TRUE)
usethis::use_data(follow_up_clinical_data, overwrite = TRUE)
usethis::use_data(biomarker_data, overwrite = TRUE)
usethis::use_data(imaging_data, overwrite = TRUE)
usethis::use_data(cognitive_composite_data, overwrite = TRUE)
usethis::use_data(ecog_scores_data, overwrite = TRUE)
usethis::use_data(RID_to_PTID, overwrite = TRUE)
usethis::use_data(ADNIMERGE, overwrite = TRUE)
usethis::use_data(visits, overwrite = TRUE)
usethis::use_data(visits, overwrite = TRUE)

## Explore the ADNI metadata ====

library(DescrTab2)
adnimerge$APOE4 <- adnimerge$APOE4 |> as.factor()
adnimerge$PTAU <- adnimerge$PTAU |> as.numeric()
adnimerge$TAU <- adnimerge$TAU |> as.numeric()
adnimerge |> dplyr::filter(VISCODE == "bl") |>
    dplyr::select(DX_bl, MMSE, AGE, PTGENDER,
                  PTEDUCAT, APOE4, ADAS13,
                  RAVLT_perc_forgetting,
                  PTAU, TAU, Hippocampus, MOCA,
                  EcogPtTotal) |>
    dplyr::filter(DX_bl != "") |>
    DescrTab2::descr(group = "DX_bl") -> ADNI_DX_bl


## Prepare `adnimerge` dataset ====

library(anscombiser)
library(MOFA2)

adnimerge <-  data.table::fread(
    input = "data-raw/ADNIMERGE_20Sep2024.csv",
    sep = ",") |> tibble::tibble()

## select variables of interest ##
# csf: ABETA_bl, PTAU_bl, TAU_bl
# pet: AV45_bl
# cognition: CDRSB, CDRSB_bl, ADAS13, ADAS13_bl, MMSE, MMSE_bl
# mri: Hippocampus_bl
# demo: AGE, PTGENDER, PTEDUCAT, APOE4
# subj: RID, VISCODE, Years_bl, DX_bl
adnimerge <- adnimerge %>%
    select(
        RID, VISCODE, Years_bl, DX_bl,
        AGE, PTGENDER, PTEDUCAT, APOE4,
        CDRSB, CDRSB_bl, ADAS13, ADAS13_bl, MMSE, MMSE_bl,
        ABETA_bl, PTAU_bl, TAU_bl,
        Hippocampus_bl,
        AV45_bl
    )

## process CSF data ##
adnimerge <- adnimerge %>%
    mutate(
        across(
            c(ABETA_bl, PTAU_bl, TAU_bl),
            ~stringr::str_replace_all(., c('<'='','>'='')) %>% as.numeric()
        )
    )

# filter
adnimerge <- adnimerge %>%
    #filter(complete.cases(.)) %>% # here I tried to remove rows with at least 1 missing value
    arrange(RID, Years_bl) %>%
    mutate(
        DX_bl = fct_recode(
            DX_bl,
            'CU' = 'CN',
            'CU' = 'SMC',
            'MCI' = 'EMCI',
            'MCI' = 'LMCI'
        )
    )

adnimerge <- adnimerge %>%
    mutate(
        PET_ABETA_STATUS_bl = as.integer(AV45_bl > 1.11),
        APOE4 = as.integer(APOE4 >= 1),
        PTGENDER = as.integer(PTGENDER == 'Male')
    )

adnimerge <- adnimerge %>%
    rename(
        PET_ABETA_bl = AV45_bl,
        CSF_ABETA_bl = ABETA_bl,
        CSF_PTAU_bl = PTAU_bl,
        CSF_TAU_bl = TAU_bl,
        GENDER = PTGENDER,
        EDUCATION = PTEDUCAT,
        YEARS_bl = Years_bl,
        MRI_HIPP_bl = Hippocampus_bl
    )


data <- adnimerge %>% dplyr::filter(VISCODE == 'bl')

# build and fit an aba model with multiple outcomes and predictors
model <- data %>% aba_model() %>%
    set_outcomes(ConvertedToAlzheimers, CSF_ABETA_STATUS_bl) %>%
    set_predictors(
        PLASMA_ABETA_bl, PLASMA_PTAU181_bl, PLASMA_NFL_bl,
        c(PLASMA_ABETA_bl, PLASMA_PTAU181_bl, PLASMA_NFL_bl)
    ) %>%
    set_stats('glm') %>%
    fit()

# summarise the model results (coefficients and metrics)
model_summary <- model_fit %>% summary()


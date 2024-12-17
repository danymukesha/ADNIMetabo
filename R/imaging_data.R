#' Imaging data
#'
#' This dataset contains imaging measurements from the Alzheimer's Disease
#' Neuroimaging Initiative (ADNI).
#'
#' @format A data frame with rows and variables as described in the ADNI
#' documentation.
#' \describe{
#'   \item{RID}{Participant roster ID}
#'   \item{FLDSTRENG}{MRI Field Strength}
#'   \item{FSVERSION}{FreeSurfer Software Version}
#'   \item{Ventricles}{UCSF Ventricles}
#'   \item{Hippocampus}{UCSF Hippocampus}
#'   \item{WholeBrain}{UCSF WholeBrain}
#'   \item{Entorhinal}{UCSF Entorhinal}
#'   \item{Fusiform}{UCSF Fusiform}
#'   \item{MidTemp}{UCSF Med Temp}
#'   \item{ICV}{UCSF ICV}
#' }
#' @source Alzheimer's Disease Neuroimaging Initiative (ADNI)
"imaging_data"

#' Load the Imaging Data
#' data("imaging_data", package = "ADNIData", envir = environment())

#' Biomarker data
#'
#' This dataset contains biomarker measurements from the Alzheimer's Disease
#' Neuroimaging Initiative (ADNI).
#'
#' @format A data frame with rows and variables as described in the ADNI
#' documentation.
#' \describe{
#'   \item{RID}{Participant roster ID}
#'   \item{ABETA}{CSF ABETA}
#'   \item{TAU}{CSF TAU}
#'   \item{PTAU}{CSF PTAU}
#'   \item{FDG}{Average FDG-PET of angular, temporal, and posterior cingulate}
#'   \item{PIB}{Average PIB SUVR of frontal cortex, anterior cingulate,
#'   precuneus cortex, and parietal cortex}
#'   \item{AV45}{Reference region - florbetapir mean of whole cerebellum}
#' }
#' @source Alzheimer's Disease Neuroimaging Initiative (ADNI)
"biomarker_data"

#' Load the Biomarker
#' data("biomarker_data", package = "ADNIData", envir = environment())

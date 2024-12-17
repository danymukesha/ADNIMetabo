#' Cognitive composite data
#'
#' This dataset contains cognitive composite scores from the Alzheimer's
#' Disease Neuroimaging Initiative (ADNI).
#'
#' @format A data frame with rows and variables as described in the ADNI
#' documentation.
#' \describe{
#'   \item{RID}{Participant roster ID}
#'   \item{mPACCdigit}{ADNI modified Preclinical Alzheimer's Cognitive
#'   Composite (PACC) with Digit Symbol Substitution}
#'   \item{mPACCtrailsB}{ADNI modified Preclinical Alzheimer's Cognitive
#'   Composite (PACC) with Trails B}
#' }
#' @source Alzheimer's Disease Neuroimaging Initiative (ADNI)
"cognitive_composite_data"

#' Load the Cognitive Composite Data
#' data("cognitive_composite_data", package = "ADNIData",
#' envir = environment())

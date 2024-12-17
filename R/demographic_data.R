#' Demographic data
#'
#' This dataset contains demographic information from the Alzheimer's Disease
#' Neuroimaging Initiative (ADNI).
#'
#' @format A data frame with rows and variables as described in the ADNI
#' documentation.
#' \describe{
#'   \item{RID}{Participant roster ID}
#'   \item{AGE}{Age}
#'   \item{PTGENDER}{Sex}
#'   \item{PTEDUCAT}{Education}
#'   \item{PTETHCAT}{Ethnicity}
#'   \item{PTRACCAT}{Race}
#'   \item{PTMARRY}{Marital status}
#' }
#' @source Alzheimer's Disease Neuroimaging Initiative (ADNI)
"demographic_data"

#' Load the Demographic Data
#' data("demographic_data", package = "ADNIData", envir = environment())

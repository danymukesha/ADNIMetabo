#' Follow-Up clinical data
#'
#' This dataset contains follow-up clinical assessments from the Alzheimer's
#' Disease Neuroimaging Initiative (ADNI).
#'
#' @format A data frame with rows and variables as described in the ADNI
#' documentation.
#' \describe{
#'   \item{RID}{Participant roster ID}
#'   \item{VISCODE}{Visit code}
#'   \item{EXAMDATE}{Date}
#'   \item{DX}{Diagnosis}
#'   \item{CDRSB}{CDR-SB}
#'   \item{ADAS11}{ADAS 11}
#'   \item{ADAS13}{ADAS 13}
#'   \item{MMSE}{MMSE}
#'   \item{RAVLT_immediate}{RAVLT Immediate}
#'   \item{RAVLT_learning}{RAVLT Learning}
#'   \item{RAVLT_forgetting}{RAVLT Forgetting}
#' }
#' @source Alzheimer's Disease Neuroimaging Initiative (ADNI)
"follow_up_clinical_data"

#' Load the Follow-Up Clinical Data
#' data("follow_up_clinical_data", package = "ADNIData",
#' envir = environment())

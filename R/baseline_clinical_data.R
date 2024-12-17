#' Baseline clinical data
#'
#' This dataset contains baseline clinical assessments from the Alzheimer's
#' Disease Neuroimaging Initiative (ADNI).
#'
#' @format A data frame with rows and variables as described in the ADNI
#' documentation.
#' \describe{
#'   \item{RID}{Participant roster ID}
#'   \item{DX_bl}{Baseline diagnosis}
#'   \item{CDRSB_bl}{CDR-SB at baseline}
#'   \item{ADAS11_bl}{ADAS 11 at baseline}
#'   \item{ADAS13_bl}{ADAS 13 at baseline}
#'   \item{MMSE_bl}{MMSE at baseline}
#'   \item{RAVLT_immediate_bl}{RAVLT Immediate at baseline}
#'   \item{RAVLT_learning_bl}{RAVLT Learning at baseline}
#'   \item{RAVLT_forgetting_bl}{RAVLT Forgetting at baseline}
#'   \item{FAQ_bl}{FAQ at baseline}
#'   \item{MOCA_bl}{MOCA at baseline}
#' }
#' @source Alzheimer's Disease Neuroimaging Initiative (ADNI)
"baseline_clinical_data"

#' Load the Baseline Clinical Data
#' data("baseline_clinical_data", package = "ADNIData", envir = environment())

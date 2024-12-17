#' ECog scores data
#'
#' This dataset contains ECog scores for both participants and study
#' partners from the Alzheimer's Disease Neuroimaging Initiative (ADNI).
#'
#' @format A data frame with rows and variables as described in the ADNI
#' documentation.
#' \describe{
#'   \item{RID}{Participant roster ID}
#'   \item{EcogPtMem}{Pt ECog - Mem}
#'   \item{EcogPtLang}{Pt ECog - Lang}
#'   \item{EcogPtVisspat}{Pt ECog - Vis/Spat}
#'   \item{EcogPtPlan}{Pt ECog - Plan}
#'   \item{EcogPtOrgan}{Pt ECog - Organ}
#'   \item{EcogPtDivatt}{Pt ECog - Div atten}
#'   \item{EcogPtTotal}{Pt ECog - Total}
#'   \item{EcogSPMem}{SP ECog - Mem}
#'   \item{EcogSPLang}{SP ECog - Lang}
#'   \item{EcogSPVisspat}{SP ECog - Vis/Spat}
#'   \item{EcogSPPlan}{SP ECog - Plan}
#'   \item{EcogSPOrgan}{SP ECog - Organ}
#'   \item{EcogSPDivatt}{SP ECog - Div atten}
#'   \item{EcogSPTotal}{SP ECog - Total}
#' }
#' @source Alzheimer's Disease Neuroimaging Initiative (ADNI)
"ecog_scores_data"

#' Load the Cognitive Composite Data
#' data("ecog_scores_data", package = "ADNIData", envir = environment())

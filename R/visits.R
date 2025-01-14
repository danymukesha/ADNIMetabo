#' Baseline clinical data
#'
#' This dataset contains visits details from the Alzheimer's
#' Disease Neuroimaging Initiative (ADNI).
#'
#' @format A data frame with rows and variables as described in the ADNI
#' documentation.
#' \describe{
#'   \item{Phase}{ADNI stady phase}
#'   \item{ID}{Unique identifier}
#'   \item{VISCODE}{Visit code}
#'   \item{VISNAME}{Visit name}
#'   \item{VISORDER}{Visit order}
#' }
#' @source Alzheimer's Disease Neuroimaging Initiative (ADNI)
"visits"

#' Load the Baseline Clinical Data
#' data("visits", package = "ADNIData", envir = environment())

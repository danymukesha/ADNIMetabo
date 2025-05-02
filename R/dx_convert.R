#' Generate diagnostic conversion indicators for ADNI Data
#'
#' @description With this function you will be able to create binary flags
#' tracking transitions between baseline and follow-up diagnoses
#' in the ADNI dataset. Optimized for performance and compatibility
#' with tidyverse workflows.
#'
#' @param dt is a data frame containing ADNI clinical data
#' @param dx_only is logical indicating whether to return
#' only diagnosis-related columns (default: TRUE)
#' @param check_validity is a logical indicating whether to validate
#' diagnostic codes (default: TRUE)
#'
#' @return The function should return a tibble with original columns
#' (if dx_only = FALSE) and added conversion indicators:
#' \itemize{
#'   \item{`ad_conv`: conversion from any non-AD baseline
#'   to dementia (1 = conversion)}
#'   \item{`cn_to_ad`: conversion from cognitively normal (CN) to dementia}
#'   \item{`cn_to_mci`: conversion from CN to mild cognitive impairment (MCI)}
#'   \item{`emci_to_ad`: conversion from early MCI to dementia}
#'   \item{`lmci_to_ad`: conversion from late MCI to dementia}
#' }
#'
#' @section Input requirements:
#' - should contain columns: RID (participant ID), VISCODE (visit code)
#' - should contain diagnostic columns: DX_bl (baseline diagnosis),
#' DX (current diagnosis)
#' - expected diagnostic values: "CN", "MCI", "EMCI", "LMCI", "Dementia"
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#'
#' # eg.
#' adni_data <- tibble(
#'   RID = c(1, 2, 3),
#'   VISCODE = c("bl", "m12", "m24"),
#'   DX_bl = c("CN", "EMCI", "LMCI"),
#'   DX = c("MCI", "Dementia", "Dementia")
#' )
#'
#' # generate conversion indicators
#' dx_convert(adni_data, dx_only = TRUE)
#' }
#'
#' @export
dx_convert <- function(data, dx_only = TRUE, check_validity = TRUE) {

    required_cols <- c("RID", "VISCODE", "DX_bl", "DX")
    if (!all(required_cols %in% colnames(data))) {
        stop("Input data missing required columns: ",
             paste(setdiff(required_cols, colnames(data)), collapse = ", "))
    }

    if (check_validity) {
        valid_dx <- c("CN", "MCI", "EMCI", "LMCI", "AD", "SMC", "Dementia", "")
        invalid_dx <- unique(c(data$DX_bl, data$DX)) %>%
            setdiff(valid_dx) %>%
            na.omit()

        if (length(invalid_dx) > 0) {
            stop(
                "Invalid diagnostic codes detected: ",
                paste(invalid_dx, collapse = ", "),
                "\nExpected values: ",
                paste(valid_dx, collapse = ", ")
            )
        }
    }

    # here i'm creating conversion indicators using vectorized operations.
    data <- data %>%
        dplyr::as_tibble() %>%
        dplyr::mutate(
            ad_conv = dplyr::if_else(
                .data$DX_bl != "AD" & .data$DX == "Dementia", 1L, 0L),
            cn_to_ad = dplyr::if_else(
                .data$DX_bl == "CN" & .data$DX == "Dementia", 1L, 0L),
            cn_to_mci = dplyr::if_else(
                .data$DX_bl == "CN" & .data$DX == "MCI", 1L, 0L),
            emci_to_ad = dplyr::if_else(
                .data$DX_bl == "EMCI" & .data$DX == "Dementia", 1L, 0L),
            lmci_to_ad = dplyr::if_else(
                .data$DX_bl == "LMCI" & .data$DX == "Dementia", 1L, 0L),
            .after = "DX"
        )

    if (dx_only) {
        data <- data %>%
            dplyr::select(
                "RID", "VISCODE", "DX_bl", "DX",
                dplyr::starts_with("ad_conv"),
                dplyr::matches("^cn_|mci_to_ad$")
            )
    }

    data
}

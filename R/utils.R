#' Process and impute data below Limits of Detection (LOD) with reporting
#'
#' this `f()` processes a dataset by identifying analytes with values
#' below the Limit of Detection (LOD), imputing these values using a logspline
#' distribution, and removing columns with specific conditions such as all
#' values being below LOD or missing. It provides a summary report
#' of the processing steps.
#'
#' @param raw_data is a data frame containing the raw data to be processed.
#' Columns should represent analytes, and rows should represent samples.
#' @param LODs is a data frame containing information about analytes and
#' their respective LOD values. Must include columns:
#'   - `Analyte (short name)`: Short names of analytes matching
#'   the column names in `raw_data`.
#'   - `LOD (µM)`: The limit of detection for each analyte.
#' @param inj_type is the injection type. e.g.: FIA or UPLC.
#'  Please provide a corrected name of injection method in quotes.
#' @param qc+samples If TRUE, you intend to process QC samples, otherwise
#'  study samples. default is FALSE.
#'
#' @return a list with two components:
#'   - `processed_data`: is a processed version of the input `raw_data`.
#'   - `report`: is a list containing details of removed and imputed columns.
#'
#' @examples
#' \dontrun{
#' # example raw data
#' raw_data <- data.frame(
#'   Analyte1 = c(0.1, 0.05, NA, 0.2),
#'   Analyte2 = c(0, 0, 0, 0),
#'   Analyte3 = c(0.3, 0.4, 0.2, 0.1)
#' )
#'
#' # example LODs table
#' LODs <- data.frame(
#'   `Analyte (short name)` = c("Analyte1", "Analyte2", "Analyte3"),
#'   `LOD (µM)` = c(0.15, 0.1, 0.25)
#' )
#'
#' # process the raw data
#' result <- process_data_with_report(raw_data, LODs)
#' }
#' @note
#' The `logspline` package doesn't support setting minimum and maximum values
#' directly in the `rlogspline` function, so instead, we needed to fit
#' a logspline model to a set of values between `LOD/2` and `LOD`.
#' Then, we can sample from this fitted distribution and clip the values
#' to stay within the desired range.
#' @importFrom logspline logspline rlogspline
process_data_with_report <- function(raw_data, LODs, inj_type = "FIA",
                                     qc_samples = FALSE) {
    library(logspline)
    col_to_remove <- c()
    removed_columns <- list()
    imputed_columns <- list()

    impute_below_lod <- function(values, lod) {
        sample_range <- seq(lod / 2, lod, length.out = 100)
        logspline_fit <- logspline(sample_range)
        sapply(values, function(x) logspline::rlogspline(1, logspline_fit))
    }

    if (inj_type == "FIA") {
        for (i in LODs$`Short.Name/Injection.Number`) {
            if (i %in% colnames(raw_data)) {
                lod <- LODs$`LOD (calc.) 1034741038/1 [µM]`[which(LODs$`Short.Name/Injection.Number` == i)]
                lod <- lod |> as.numeric()
                col_pos <- which(colnames(raw_data) == i)
                raw_data[, col_pos] <- ifelse(
                    qc_samples,
                    raw_data[, col_pos] |>
                        dplyr::mutate_if(is.character, as.numeric),
                    raw_data[, col_pos] |>
                        dplyr::mutate_if(is.integer, as.numeric)
                )

                row_pos_lowLOD <- which(raw_data[, col_pos] <= lod)
                if (length(row_pos_lowLOD) > 0) {
                    ki <- length(row_pos_lowLOD) + sum(is.na(raw_data[, col_pos]))

                    if (ki == nrow(raw_data) && !qc_samples) {
                        col_to_remove <- c(col_to_remove, col_pos)
                        removed_columns[[i]] <- "All values below LOD or missing"
                    } else {
                        num_imputed <- length(row_pos_lowLOD)
                        raw_data[row_pos_lowLOD, col_pos] <-
                            impute_below_lod(raw_data[row_pos_lowLOD, col_pos], lod)
                        imputed_columns[[i]] <- num_imputed
                    }
                }
            } else {
                # handle columns with numeric data containing only zeros
                col_pos <- which(colnames(raw_data) == i)
                if (is.numeric(unlist(raw_data[, col_pos])) && sum(raw_data[, col_pos],
                                                                   na.rm = TRUE) == 0) {
                    col_to_remove <- c(col_to_remove, col_pos)
                    removed_columns[[i]] <- "All values are zero"
                }
            }
        }
    } else if (inj_type == "UPLC") {
        for (i in LODs$`Short.Name/Injection.Number`) {
            if (i %in% colnames(raw_data)) {
                lod <- LODs$`LOD (from OP) [?M]`[which( LODs$`Short.Name/Injection.Number` == i)]
                lod <- lod |> as.numeric()
                col_pos <- which(colnames(raw_data) == i)

                row_pos_lowLOD <- which(raw_data[, col_pos] == "<LOD")
                if (length(row_pos_lowLOD) > 0) {
                    ki <- length(row_pos_lowLOD) + sum(is.na(raw_data[, col_pos]))

                    if (ki == nrow(raw_data) && !qc_samples) {
                        col_to_remove <- c(col_to_remove, col_pos)
                        removed_columns[[i]] <- "All values below LOD or missing"
                    } else {
                        num_imputed <- length(row_pos_lowLOD)
                        raw_data[, col_pos] <- raw_data[, col_pos] |>
                            dplyr::mutate(across(everything(),
                                                 ~ replace(., . ==  "<LOD" , lod))) |>
                            dplyr::mutate_if(is.character, as.numeric)
                        imputed_columns[[i]] <- num_imputed
                    }
                }
            } else {
                # handle columns with numeric data containing only zeros
                col_pos <- which(colnames(raw_data) == i)
                if (is.numeric(unlist(raw_data[, col_pos])) && sum(raw_data[, col_pos],
                                                                   na.rm = TRUE) == 0) {
                    col_to_remove <- c(col_to_remove, col_pos)
                    removed_columns[[i]] <- "All values are zero"
                }
            }
        }

    } else {
        warning("Please provide a valide injuctation method name.
                eg: 'FIA' or 'UPLC'.")
    }

    if (length(col_to_remove) > 0) {
        raw_data <- raw_data[, -col_to_remove, drop = FALSE]
    }

    report <- list(
        removed_columns = removed_columns,
        imputed_columns = imputed_columns,
        original_data = raw_data
    )

    cat("Processing Summary:\n")
    if (length(removed_columns) > 0) {
        cat("Removed Columns:\n")
        for (col_name in names(removed_columns)) {
            cat(sprintf("- %s: %s\n", col_name, removed_columns[[col_name]]))
        }
    } else {
        cat("No columns were removed.\n")
    }
    message("Total removed: ", length(removed_columns))

    if (length(imputed_columns) > 0) {
        cat("\nImputed Columns:\n")
        for (col_name in names(imputed_columns)) {
            cat(sprintf("- %s: %d values imputed\n", col_name,
                        imputed_columns[[col_name]]))
        }
    } else {
        cat("\nNo values were imputed.\n")
    }
    message("Total imputed: ", length(imputed_columns))

    return(list(processed_data = raw_data, report = report))
}


#' Add missing columns to data frame
#'
#' this `f()` ensures a data frame contains specified columns,
#' adding them as `NA` if missing.
#'
#' @param data Data frame to check.
#' @param cname Vector of column names to add if missing.
#' @return A data frame with added columns (if missing).
#' @examples
#' \dontrun{
#' data <- fncols(data, c("new_col1", "new_col2"))
#' }
fncols <- function(data, cname) {
    add <- cname[!cname %in% names(data)]

    if (length(add) != 0) data[add] <- NA
    data
}



#' Add DX_bl, AGE, SEX, ... data
#'
#' this `f()` adds DX_bl, AGE, SEX, ... information to the dataset
#' based on sample identifiers.
#'
#' @param mt Main dataset.
#' @param other_info Additional information containing DX_bl, AGE, SEX, ... data.
#' @return Data frame with added DX_bl, AGE, SEX, ... information.
#' @examples
#' \dontrun{
#' updated_data <- add_metadata(mt, other_info)
#' }
add_metadata <- function(mt,
                     other_info) {
    samples <- Filter(nzchar, mt$RID)
    ## Disease baseline
    mt$DX_bl <- NA
    ## Age
    mt$AGE <- NA
    ## Patient gender
    mt$SEX <- NA
    ## APOE4
    mt$APOE4 <- NA
    ## ABETA baseline
    mt$ABETA <- NA
    ## TAU baseline
    mt$TAU <- NA
    ## PTAU baseline
    mt$PTAU <- NA
    ## MMSE baseline
    ##
    ## MOCA baseline
    ##
    ## ADAS13 baseline
    ##
    ## Intracerebroventricular injection (ICV)
    ##
    ## Average FDG-PET of angular, temporal, and posterior cingulate
    for (sample in samples) {
        pattern <- sample
        k <- which(other_info$RID %in% pattern)
        j <- which(grepl(paste0("\\b", sample, "\\b"), mt$RID))
        visits <- mt |> dplyr::slice(j) |> select(VISCODE2)
        allSame <- function(x) length(unique(x)) == 1
        if (length(j) > 0) {
            # Disease
            DX_bl <- other_info[k, ]$DX_bl
            if (length(DX_bl) == 0) {
                ms <- mt[j, ]$RID
            }
            if (allSame(DX_bl)) {
                mt[j, ]$DX_bl <- DX_bl[1]
            }
            # AGE
            AGE <- other_info[k, ]$AGE
            if (length(AGE) == 0) {
                ms <- mt[j, ]$RID
            }
            if (allSame(AGE)) {
                mt[j, ]$AGE <- AGE[1]
            }
            # SEX
            SEX <- other_info[k, ]$PTGENDER
            if (length(SEX) == 0) {
                ms <- mt[j, ]$RID
            }
            if (allSame(SEX)) {
                mt[j, ]$SEX <- SEX[1]
            }
            # ABETA baseline, and other timepoints
            ABETA <- other_info[k, ]$ABETA
            if (length(ABETA) == 0) {
                ms <- mt[j, ]$RID
            }
            ABETA <- other_info |> dplyr::slice(k) |>
                dplyr::select(RID, VISCODE, ABETA) |>
                dplyr::mutate(across(ABETA, gsub,
                              pattern = "[^0-9.-]",
                              replacement = "")) |>
                dplyr::mutate_at("ABETA", as.numeric)

            mt[j, "ABETA"] <- ABETA |>
                dplyr::filter(VISCODE %in% list(visits$VISCODE2)[[1]]) |>
                dplyr::select(ABETA)

            # TAU baseline and other timepoints
            TAU <- other_info[k, ]$TAU
            if (length(TAU) == 0) {
                ms <- mt[j, ]$RID
            }
            TAU <- other_info |> dplyr::slice(k) |>
                dplyr::select(RID, VISCODE, TAU) |>
                dplyr::mutate(across(TAU, gsub,
                                     pattern = "[^0-9.-]",
                                     replacement = "")) |>
                dplyr::mutate_at("TAU", as.numeric)

            mt[j, "TAU"] <- TAU |>
            dplyr::filter(VISCODE %in% list(visits$VISCODE2)[[1]]) |>
            dplyr::select(TAU)
            # PTAU baseline
            PTAU <- other_info[k, ]$TAU
            if (length(PTAU) == 0) {
                ms <- mt[j, ]$RID
            }
            PTAU <- other_info |> dplyr::slice(k) |>
                dplyr::select(RID, VISCODE, PTAU) |>
                dplyr::mutate(across(PTAU, gsub,
                                     pattern = "[^0-9.-]",
                                     replacement = "")) |>
                dplyr::mutate_at("PTAU", as.numeric)

            mt[j, "PTAU"] <- PTAU |>
            dplyr::filter(VISCODE %in% list(visits$VISCODE2)[[1]]) |>
            dplyr::select(PTAU)
            # APOE4
            APOE4 <- other_info[k, ]$APOE4
            if (length(APOE4) == 0) {
                ms <- mt[j, ]$RID
            }
            if (allSame(APOE4)) {
                mt[j, ]$APOE4 <- APOE4[1]
            }

        }
    }
    mt <- mt |>
        dplyr::relocate(DX_bl, .after = VISCODE2) |>
        dplyr::relocate(SEX, .after = DX_bl)  |>
        dplyr::relocate(AGE, .after = SEX) |>
        dplyr::relocate(ABETA, .after = AGE) |>
        dplyr::relocate(TAU, .after = ABETA) |>
        dplyr::relocate(PTAU, .after = TAU) |>
        dplyr::relocate(APOE4, .after = PTAU)

    return(mt)
    message(sum(mt$DX_bl != "N/A"),
            " out of ",
            nrow(mt),
            " samples have DX_baseline")
}

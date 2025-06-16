#' batch normalization for a tibble
#'
#' Normalization data by the median value of each batch in a tibble
#'
#' @param tbl is a tibble object. The first two columns are "ID" and "Batch", the remaining columns for batch normalization
#' @param try test the function for the first 20 columns.
#' @param verbose print log information.
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @return A data.frame object after normalization.
#' @export
batch_normalization <- function(tbl, try = FALSE,
                                batch_name = "Plate Bar Code",
                                verbose = TRUE) {
    stopifnot(inherits(tbl, "list"))
    output_data <- copy(tbl)
    input_data <- copy(tbl)

    platforms <- tbl[, batch_name]

    if(verbose) {
        cat("\nBatch information:\n")
        print(table(platforms))
        cat("\n")
    }

    v_n_col <- dim(tbl)[2]

    if(isTRUE(test)) {
        v_n_col <- 10
    } else if(test >=2 ) {
        v_n_col <- test
    }

    pb <- txtProgressBar(min = 0, max = v_n_col, style = 3, file = stderr())

    sample_IDs <- tbl[, "RID"]

    v_batch <- platforms

    for(j in 5L:v_n_col) {
        setTxtProgressBar(pb = pb, value = j)

        v_metab <- names(tbl)[j]
        InputData_each <- unlist(input_data[, j, with = FALSE])


        for (id_batch in unique(v_batch[! is.na(v_batch)])) {

            Index_each <- (v_batch  == id_batch)

            v_median <- median(InputData_each[Index_each], na.rm = TRUE)

            set(OutputData, which(Index_each), j, InputData_each[which(Index_each)]/v_median)
        }
    }

    close(con = pb)

    return(OutputData)
}

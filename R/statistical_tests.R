#' @title student_fdr
#' @description
#' Calculate a Student test and the foldchange
#'
#' @param start Numeric, the starting column index of your matrix
#' or data frame.
#' @param end Numeric, the ending column index of your matrix or data frame.
#' @param group_1 Data frame or matrix for the first group.
#' @param group_2 Data frame or matrix for the second group.
#' @param log Logical, TRUE if your data are log-transformed.
#' @param base Numeric, base of the log transformation
#' @param group_1_original first group group original data set
#' @param group_2_original second group original data set
#' (ignored if log = FALSE).
#' @param comparison Character, name of the comparison for the output file.
#' @return the function returns a data frame with mean, standard deviation,
#' p-value, fold-change, and FDR-adjusted p-value for each variable.
#' @export
student_fdr <- function(start, end, group_1, group_2, log, base, comparison,
                        group_1_original = NULL, group_2_original = NULL) {

    # argument validation
    if (!is.numeric(start) || start < 1)
        stop("'start' must be a positive numeric value.")
    if (!is.numeric(end) || end <= start)
        stop("'end' must be greater than 'start'.")
    if (!is.matrix(group_1) && !is.data.frame(group_1))
        stop("'group_1' must be a matrix or data frame.")
    if (!is.matrix(group_2) && !is.data.frame(group_2))
        stop("'group_2' must be a matrix or data frame.")
    if (!is.logical(log)) stop("'log' must be a logical value.")
    if (!is.numeric(base) || base <= 0)
        stop("'base' must be a positive numeric value.")
    if (!is.character(comparison))
        stop("'comparison' must be a character string.")

    # identify and filter out columns with constant data
    constant_columns <- which(
        apply(group_1, 2, sd, na.rm = TRUE) == 0 |
            apply(group_2, 2, sd, na.rm = TRUE) == 0
    )
    if (length(constant_columns) > 0) {
        message("Skipping columns with constant data at positions: ",
                paste(constant_columns, collapse = ", "))
        group_1 <- group_1[, -constant_columns]
        group_2 <- group_2[, -constant_columns]
        # adjust the end index after removing columns
        end <- end - length(constant_columns)
    }

    # start of the function
    test <- NULL
    for (i in start:end) {

        # calculate the P-value form a Student test and the mean and SD for each
        pvalue <- NULL
        pvalue <- c(
            variable <- colnames(group_2)[i],
            moygroup2 <- mean(group_2[, i], na.rm = TRUE),
            sdgroup2 <- sd(group_2[, i], na.rm = TRUE),
            moygroup1 <- mean(group_1[, i], na.rm = TRUE),
            sdgroup1 <- sd(group_1[, i], na.rm = TRUE),
            pvalue <- t.test(group_1[, i], group_2[, i], na.rm = TRUE)$p.value,
            statistic <- t.test(group_1[, i], group_2[, i], na.rm = TRUE)$statistic |>
                as.double()
        )
        # calculate a fold-change depending of we do a log transformation or not
        foldchange <- NULL
        if (log == TRUE) {
            foldchange <- base^(mean(group_1[,i],na.rm=T)-mean(group_2[,i],na.rm=T))
        } else {
            foldchange <- mean(group_1[, i], na.rm = TRUE) / mean(group_2[, i],
                                                                  na.rm = TRUE)
        }
        pvalue <- c(pvalue, foldchange <- foldchange)
        test <- rbind(test, pvalue)
    }
    test=cbind(test,p.adjust(test[,6],"fdr"))
    if (is.null(group_2_original) && is.null(group_1_original)) {
        test <- cbind(test, NA)
        test <- cbind(test, NA)
    } else {
        test <- cbind(test, sapply(group_2_original[,start:end], mean) |>
                          as.data.frame())
        test <- cbind(test, sapply(group_1_original[,start:end], mean) |>
                          as.data.frame())
    }
    colnames(test) <- c(
        "variable",
        "mean group2",
        "sd group2",
        "mean group1",
        "sd group1",
        "pvalue",
        "statistic",
        "foldchange",
        "pvalue adjusted",
        "mean group2 original",
        "mean group1 original"
    )
    folder_name <- "student_test"
    file_name <- paste0("Summary_", comparison, ".csv")
    file_path <- file.path(folder_name, file_name)

    if (!dir.exists(folder_name)) {
        dir.create(folder_name)
    }
    write.csv(test, file_path, row.names = FALSE)
    return(test)
}

#' @title wilcox_fdr
#' @description
#' this `f()` calculates a Wilcoxon test and the foldchange
#'
#' @param start the start of your matrix or data frame
#' @param end the end of the matrix or data frame
#' @param group_1 data frame or matrix for  the first group of your comparison
#' for the Wilcoxon test
#' @param group_2 data frame or matrix for  the second group of your comparison
#' for the Wilcoxon test
#' @param log logical TRUE if your data are log transformed
#' @param base base for your log transformation
#' @param comparison name of your comparison that it will the title
#' of your output
#' @return the function give a csv file with a mean, the standard deviation,
#' the p-value of a Wilcoxon test and the foldchange for all variable
#' in your matrix or data frame
#' @export


wilcox_fdr <- function(start, end, group_1, group_2, log, base, comparison) {

    # test for the argument
    if (is.numeric(start) == FALSE) {
        stop("The argument start need to be a numeric argument")}
    if (is.numeric(end) == FALSE) {
        stop("The argument end need to be a numeric argument")}
    if (end == start) {
        stop("You need to have at least 2 variable in your data frame or matrix")}
    if (is.matrix(group_1) == FALSE & is.data.frame(group_1) == FALSE) {
        stop("The argument group_1 need to be a data frame or a matrix")}
    if (is.matrix(group_2) == FALSE & is.data.frame(group_2) == FALSE) {
        stop("The argument group_2 need to be a data frame or a matrix")}
    if (is.logical(log) == FALSE) {
        stop("The argument log need to be a logical argument")}
    if (is.numeric(base) == FALSE) {
        stop("The argument base need to be a numeric argument")}
    if (is.character(comparison) == FALSE) {
        stop("The argument comparison need to be a characther")}

    # start of the function
    test <- NULL
    for (i in start:end) {
        # Calculation of the P-value form a Wilcoxon test and the mean and SD for each
        pvalue <- NULL
        pvalue <- c(
            variable <- colnames(group_2)[i],
            moygroup2 <- mean(group_2[, i], na.rm = TRUE),
            sdgroup2 <- sd(group_2[, i], na.rm = TRUE),
            moygroup1 <- mean(group_1[, i], na.rm = TRUE),
            sdgroup1 <- sd(group_1[, i], na.rm = TRUE),
            pvalue <- wilcox.test(group_1[, i], group_2[, i], na.rm = TRUE)$p.value
        )
        # Calculation of a foldchange depending f we do a log transformation or not
        foldchange <- NULL
        if (log == TRUE) {
            foldchange <- base ^ (mean(group_1[, i], na.rm = T) - mean(group_2[, i], na.rm =
                                                                           T))
        } else {
            foldchange <- mean(group_1[, i], na.rm = TRUE) / mean(group_2[, i], na.rm = TRUE)
        }
        pvalue <- c(pvalue, foldchange <- foldchange)
        test <- rbind(test, pvalue)
    }
    test = cbind(test, p.adjust(test[, 6], "fdr"))
    colnames(test) <- c(
        "variable",
        "mean group2",
        "sd group2",
        "mean group1",
        "sd group1",
        "pvalue",
        "foldchange",
        "pvalue adjusted"
    )
    folder_name <- "wilcox_fdr"
    file_name <- paste0("Summary_", comparison, ".csv")
    file_path <- file.path(folder_name, file_name)

    if (!dir.exists(folder_name)) {
        dir.create(folder_name)
    }
    write.csv(test, file_path, row.names = FALSE)
    return(test)
}

#' @title welch_fdr
#' @description
#' Calculate a Welch t-test and the fold-change
#'
#' @param start Numeric, the starting column index of your matrix
#' or data frame.
#' @param end Numeric, the ending column index of your matrix or data frame.
#' @param group_1 Data frame or matrix for the first group.
#' @param group_2 Data frame or matrix for the second group.
#' @param log Logical, TRUE if your data are log-transformed.
#' @param base Numeric, base of the log transformation
#' (ignored if log = FALSE).
#' @param comparison Character, name of the comparison for the output file.
#' @return The function returns a data frame with mean, standard deviation,
#' p-value, fold-change, and FDR-adjusted p-value for each variable.
#' @export
welch_fdr <- function(start, end, group_1, group_2, log, base, comparison) {

    # argument validation
    if (!is.numeric(start) || start < 1)
        stop("'start' must be a positive numeric value.")
    if (!is.numeric(end) || end <= start)
        stop("'end' must be greater than 'start'.")
    if (!is.matrix(group_1) && !is.data.frame(group_1))
        stop("'group_1' must be a matrix or data frame.")
    if (!is.matrix(group_2) && !is.data.frame(group_2))
        stop("'group_2' must be a matrix or data frame.")
    if (!is.logical(log)) stop("'log' must be a logical value.")
    if (!is.numeric(base) || base <= 0)
        stop("'base' must be a positive numeric value.")
    if (!is.character(comparison))
        stop("'comparison' must be a character string.")

    # identify and filter out columns with constant data
    constant_columns <- which(
        apply(group_1, 2, sd, na.rm = TRUE) == 0 |
            apply(group_2, 2, sd, na.rm = TRUE) == 0
    )
    if (length(constant_columns) > 0) {
        message("Skipping columns with constant data at positions: ",
                paste(constant_columns, collapse = ", "))
        group_1 <- group_1[, -constant_columns]
        group_2 <- group_2[, -constant_columns]
        # adjust the end index after removing columns
        end <- end - length(constant_columns)
    }

    # start of the function
    test <- NULL
    for (i in start:end) {

        # calculate the p-value from a Welch t-test and the mean and SD for each
        pvalue <- NULL
        pvalue <- c(
            variable <- colnames(group_2)[i],
            moygroup2 <- mean(group_2[, i], na.rm = TRUE),
            sdgroup2 <- sd(group_2[, i], na.rm = TRUE),
            moygroup1 <- mean(group_1[, i], na.rm = TRUE),
            sdgroup1 <- sd(group_1[, i], na.rm = TRUE),
            pvalue <- t.test(group_1[, i], group_2[, i], var.equal = FALSE, na.rm = TRUE)$p.value
        )

        # calculate a fold-change depending on whether data are log-transformed
        foldchange <- NULL
        if (log == TRUE) {
            foldchange <- base^(mean(group_1[, i], na.rm = TRUE) - mean(group_2[, i], na.rm = TRUE))
        } else {
            foldchange <- mean(group_1[, i], na.rm = TRUE) / mean(group_2[, i], na.rm = TRUE)
        }
        pvalue <- c(pvalue, foldchange <- foldchange)
        test <- rbind(test, pvalue)
    }
    test <- cbind(test, p.adjust(test[, 6], "fdr"))
    colnames(test) <- c(
        "variable",
        "mean group2",
        "sd group2",
        "mean group1",
        "sd group1",
        "pvalue",
        "foldchange",
        "pvalue adjusted"
    )
    folder_name <- "welch_test"
    file_name <- paste0("Summary_", comparison, ".csv")
    file_path <- file.path(folder_name, file_name)

    if (!dir.exists(folder_name)) {
        dir.create(folder_name)
    }
    write.csv(test, file_path, row.names = FALSE)
    return(test)
}


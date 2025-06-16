#' @title AUC_indiv
#' @description Calculates and visualizes the Area Under the Curve (AUC)
#' for individual biomarkers.
#'
#' @details This function computes the AUC for each biomarker in a dataset
#' using the ROCR package, generates sensitivity-specificity plots for each
#' biomarker, and saves the results to a PDF and CSV file. If the ROCR package
#' is not installed, it will be automatically installed and loaded.
#'
#' @param dataX is a data frame or matrix containing the biomarker values.
#' Each column corresponds to a biomarker.
#' @param dataY is a vector containing the binary class labels (e.g., 0 and 1).
#' @param n Numeric, the number of biomarkers (columns in `dataX`) to analyze.
#' @param name Character, a name to identify the output files (PDF and CSV).
#' @return this `f()` returns data frame with the AUC values for each biomarker
#' @export
#' @importFrom ROCR prediction performance
#' @examples
#' \dontrun{
#' # Example usage:
#' # Assuming `biomarkers` is a matrix of biomarker data and `labels` is a binary vector:
#' # AUC_indiv(biomarkers, labels, n = ncol(biomarkers), name = "Biomarker_AUC")
#' }
AUC_indiv = function(dataX, dataY, n, name) {
  auc = NULL
  nom = paste("AUC indiv ", name, ".pdf", sep = "")
  pdf(file = nom)
  for (i in 1:n) {
    dat = na.omit(cbind(dataX[, i], dataY))
    pred = ROCR::prediction(dat[, 1], dat[, 2])
    perf = performance(pred, "sens", "spec")
    perf2 = performance(pred, "auc")
    auc = rbind(auc, c(biomarker = colnames(dataX)[i],
                       auc = perf2@y.values[[1]]))
    plot(1 - perf@x.values[[1]],
         perf@y.values[[1]],
         main = colnames(dataX)[i],
         type = 'l',
         xlab = "1-specificity",
         ylab = "sensitivity"
    )
    legend("bottomright",
           legend = paste("AUC=", c(round(
             perf2@y.values[[1]], 3
           )), sep = ""),
           bty = "n")
  }
  folder_name <- "AUC_ROC"
  file_name <- paste("AUC indiv ", name, ".csv", sep = "")
  file_path <- file.path(folder_name, file_name)

  if (!dir.exists(folder_name)) {
    dir.create(folder_name)
  }
  write.csv(auc, file_path)
  dev.off()
  return(auc)
}

# function to calcuAUC for the Precision-Recall Curve (PR-AUC) instead of the ROC curve
AUC_indiv_PR <- function(dataX, dataY, n, name) {
  auc <- NULL
  nom <- paste("AUC_PR indiv ", name, ".pdf", sep = "")
  pdf(file = nom)

  for (i in 1:n) {
    dat <- na.omit(cbind(dataX[, i], dataY))
    pred <- ROCR::prediction(dat[, 1], dat[, 2])

    # Calculate Precision-Recall performance
    perf <- ROCR::performance(pred, "prec", "rec")
    perf2 <- ROCR::performance(pred, "aucpr")  # AUC for Precision-Recall curve
    auc <- rbind(auc, c(biomarker = colnames(dataX)[i], auc = perf2@y.values[[1]]))
    plot(perf@x.values[[1]], perf@y.values[[1]],
         main = colnames(dataX)[i], type = 'l',
         xlab = "Recall", ylab = "Precision")
    legend("bottomright", legend = paste("AUC-PR =", round(perf2@y.values[[1]], 3)), bty = "n")
  }
  nom1 <- paste("AUC_PR indiv ", name, ".csv", sep = "")
  write.csv(auc, nom1)

  dev.off()
  return(auc)
}

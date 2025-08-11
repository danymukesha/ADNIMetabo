#' Plot and Save a ROC Curve with AUC Annotation
#'
#' Generates and saves a ggplot2-based ROC curve from a `pROC::roc` object,
#' with the AUC displayed on the plot.
#'
#' @param roc_obj A `pROC::roc` object.
#' @param filename Output filename for the saved plot (e.g., "roc_curve.png").
#' @param title Title of the plot.
#' @param width Width of the saved image in inches.
#' @param height Height of the saved image in inches.
#' @param line_color Color of the ROC curve line.
#' @param auc_digits Number of digits to round the AUC value.
#' @param auc_position A numeric vector of length 2 indicating the (x, y)
#' position of the AUC label (in plot coordinates).
#' @param smooth Logical; if `TRUE`, the ROC curve will be smoothed (unless
#' @param ci Logical; if TRUE, plots confidence interval around the ROC curve.
#' @param ci_level Confidence level for interval (default 0.95).
#' @param ci_size Line size for confidence interval ribbon.
#' @param ci_alpha Transparency for confidence interval ribbon.
#' @param ci_auc if TRUE, add confidence interval for AUC
#' already).
#'
#' @return The ggplot object (invisible).
#' @export
#'
#' @examples
#' \dontrun{
#' roc_obj <- pROC::roc(
#'     response = test$Allgr, predictor = top_probs, levels = c("CN", "AD")
#' )
#' plot_and_save_roc(roc_obj, "roc_curve_top.png", title = "Top 5 Metabolites")
#' }
plot_and_save_roc <- function(roc_obj,
                              filename = "roc_curve.png",
                              title = "ROC Curve",
                              width = 5,
                              height = 5,
                              line_color = "blue",
                              auc_digits = 2,
                              auc_position = c(0.6, 0.1),
                              smooth = FALSE,
                              ci = FALSE,
                              ci_level = 0.95,
                              ci_size = 0.5,
                              ci_alpha = 0.3,
                              ci_auc = TRUE) {
    # if (!inherits(roc_obj, "roc")) {
    #     stop("roc_obj must be a pROC::roc object.")
    # }
    if (smooth && !inherits(roc_obj, "smooth.roc")) {
        roc_obj <- pROC::smooth(roc_obj)
    }
    if (is.null(roc_obj$specificities) || is.null(roc_obj$sensitivities)) {
        stop("Invalid ROC object: specificities or sensitivities not found.")
    }

    if (ci) {
        ci_obj <- pROC::ci.se(roc_obj,
            specificities = roc_obj$specificities,
            conf.level = ci_level
        )
        ci_df <- data.frame(
            fpr = rev(roc_obj$specificities),
            lower = rev(ci_obj[, 1]),
            upper = rev(ci_obj[, 3])
        )

        p <- p + ggplot2::geom_ribbon(
            data = ci_df,
            aes(x = 1 - fpr, ymin = lower, ymax = upper),
            fill = line_color, alpha = ci_alpha
        )
    }

    auc_value <- pROC::auc(roc_obj)
    if (ci_auc) {
        auc_ci <- pROC::ci.auc(roc_obj)
        auc_label <- paste0(
            "AUC = ",
            formatC(auc_value, format = "f", digits = auc_digits),
            " [",
            formatC(auc_ci[1], format = "f", digits = auc_digits),
            "–",
            formatC(auc_ci[3], format = "f", digits = auc_digits),
            "]"
        )
    }

    roc_df <- data.frame(
        fpr = rev(roc_obj$specificities),
        tpr = rev(roc_obj$sensitivities)
    )

    p <- ggplot2::ggplot(roc_df, ggplot2::aes(x = 1 - fpr, y = tpr)) +
        ggplot2::geom_line(color = line_color, size = 1) +
        ggplot2::geom_abline(
            intercept = 0, slope = 1, color = "darkgrey",
            linetype = "dashed"
        ) +
        ggplot2::annotate("text",
            x = auc_position[1],
            y = auc_position[2],
            label = auc_label, size = 5, hjust = 0
        ) +
        ggplot2::labs(
            title = title, x = "False Positive Rate",
            y = "True Positive Rate"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5),
            axis.text.y = element_text(size = 13, angle = 90, hjust = 0.5),
            axis.text.x = element_text(size = 13),
            axis.title = element_text(size = 13)
        )

    ggplot2::ggsave(filename, plot = p, width = width, height = height)

    invisible(p)
}

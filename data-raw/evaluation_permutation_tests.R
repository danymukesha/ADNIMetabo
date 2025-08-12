# Load required libraries
library(caret)
library(pROC)
library(dplyr)
library(doParallel)
library(data.table)
library(readr)

library(progressr)
handlers(global = TRUE)
handlers("txtprogressbar")

here::dr_here()
unregister_dopar <- function() {
    env <- foreach:::.foreachGlobals
    rm(list = ls(name = env), pos = env)
}

data0 <- read_csv("onlyselectedADvsCN.csv")
final <- read_csv("../lasso_varimpAD vs CN.csv")
selected_features <- final |>
    filter(numberoftimes >= 70) |>
    pull(1) |>
    as.character()

set.seed(123)
CN <- data0 |> filter(Allgr == "CN")
AD <- data0 |> filter(Allgr == "AD")

train_nCN <- sample(nrow(CN), nrow(CN) * 0.7)
train_nAD <- sample(nrow(AD), nrow(AD) * 0.7)

train_CN <- CN[train_nCN, ]
train_AD <- AD[train_nAD, ]

train <- rbind(train_CN, train_AD)

test_CN <- CN[-train_nCN, ]
test_AD <- AD[-train_nAD, ]
test <- rbind(test_CN, test_AD)

train$Allgr <- as.factor(train$Allgr)
test$Allgr <- as.factor(test$Allgr)

set.seed(123)
n_cv <- 5
n_repeat <- 20
n_param_max <- 500
my_folds <- createMultiFolds(train$Allgr, k = n_cv, times = n_repeat)
seeds <- list()

for (i in seq_len(n_cv * n_repeat)) {
    seeds[[i]] <- sample.int(n = 1000, n_param_max)
}
seeds[[n_cv * n_repeat + 1]] <- sample.int(1000, 1)

control <- trainControl("repeatedcv",
    index = my_folds,
    selectionFunction = "best",
    classProbs = TRUE, savePredictions = "final",
    summaryFunction = twoClassSummary, seeds = seeds
)
lasso_path <- "cache/lasso_model.rds"
if (file.exists(lasso_path)) {
    cat("Loading saved LASSO model...\n")
    lasso_model <- readRDS(lasso_path)
} else {
    cat("Training LASSO model...\n")
    cl <- makeCluster(detectCores() - 5)
    registerDoParallel(cl)
    lasso_model <- train(train[, selected_features], train$Allgr,
        method = "glmnet", metric = "ROC",
        trControl = control,
        preProcess = c("zv", "center", "scale", "medianImpute")
    )
    stopCluster(cl)
    unregister_dopar()
    dir.create("cache", showWarnings = FALSE)
    saveRDS(lasso_model, lasso_path)
}

# Section: AUC Evaluation with Increasing Metabolites

# Get feature importance from LASSO model
lasso_imp <- varImp(lasso_model)$importance
ranked_features <- rownames(lasso_imp)[order(lasso_imp$Overall, decreasing = TRUE)]

# Limit to 57 metabolites
ranked_features <- ranked_features[1:min(196, length(ranked_features))]

auc_path <- "cache/auc_increase.rds"
if (file.exists(auc_path)) {
    cat("Loading saved AUC progression...\n")
    auc_increase <- readRDS(auc_path)
} else {
    cat("Training models with increasing metabolites...\n")
    auc_increase <- numeric(length(ranked_features))
    cl <- makeCluster(detectCores() - 5)
    registerDoParallel(cl)
    pb <- txtProgressBar(min = 0, max = length(ranked_features), style = 3)

    for (i in seq_along(ranked_features)) {
        current_features <- ranked_features[1:i]
        x_train <- train[, current_features, drop = FALSE]
        x_test <- test[, current_features, drop = FALSE]
        x_train$dummy <- "cdsa"
        x_test$dummy <- "dasdas"
        model <- train(
            x = x_train,
            y = train$Allgr,
            method = "glmnet",
            trControl = control,
            preProcess = c("zv", "center", "scale")
        )
        probs <- predict(model, x_test, type = "prob")[, "AD"]
        roc_obj <- roc(test$Allgr, probs, levels = c("CN", "AD"), direction = "<")
        auc_increase[i] <- pROC::auc(roc_obj)
        setTxtProgressBar(pb, i)
    }
    close(pb)
    stopCluster(cl)
    unregister_dopar()
    saveRDS(auc_increase, auc_path)
}

cat("....
    Finished increasing number of metabolites...\n\n")

cat("Started permutations...\n")

dir.create("results", showWarnings = FALSE)
png("results/auc_vs_metabolites.png", width = 7.5, height = 7.5, units = "in", res = 300)
plot(1:length(auc_increase), auc_increase,
    type = "l",
    xlab = "Number of Metabolites", ylab = "AUC",
    main = "AUC vs. Number of Metabolites",
    cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.5
)
optimal_index <- which.max(auc_increase)
optimal_auc <- auc_increase[optimal_index]
abline(v = optimal_index, col = "red", lwd = 2, lty = 2)
text(optimal_index + 0.1, max(auc_increase) - 0.036,
    labels = paste0("Optimal \n(", optimal_index, " Metabolites) \nAUC:", round(optimal_auc, 2)),
    col = "red", pos = 4, cex = 1.1
)
dev.off()

# Section: Permutation Test

# Number of permutations (optimized for reasonable runtime)
perm_path <- "cache/permutation_aucs.rds"
if (file.exists(perm_path)) {
    cat("Loading saved permutation AUCs...\n")
    perm_aucs <- readRDS(perm_path)
} else {
    cat("Running permutation test...\n")
    n_perm <- 100
    log_progress <- function(step, total_steps) {
        cat(sprintf("Step %d of %d at %s\n", step, total_steps, Sys.time()),
            file = "progress_log.txt", append = TRUE
        )
    }

    cl <- makeCluster(detectCores() - 5)
    registerDoParallel(cl)
    with_progress({
        p <- progressor(steps = n_perm)
        perm_aucs <- foreach(
            i = 1:n_perm, .combine = "c",
            .packages = c("caret", "pROC", "glmnet")
        ) %dopar% {
            p(sprintf("Permutation %d", i))
            log_progress(i, n_perm)
            set.seed(i)
            perm_labels <- sample(train$Allgr)
            model_perm <- train(
                x = train[, ranked_features], y = perm_labels,
                method = "glmnet", trControl = control,
                preProcess = c("zv", "center", "scale")
            )
            probs_perm <- predict(model_perm, test[, ranked_features], type = "prob")[, "AD"]
            pROC::auc(roc(test$Allgr, probs_perm, levels = c("CN", "AD"), direction = "<"))
        }
    })
    stopCluster(cl)
    unregister_dopar()
    saveRDS(perm_aucs, perm_path)
}
# Calculate p-value
actual_auc <- auc_increase[length(auc_increase)]
p_value <- mean(perm_aucs >= actual_auc)
cat("Permutation test p-value:", p_value, "\n")

# Plot permutation test results
png("results/permutation_test.png", width = 7.5, height = 7.5, units = "in", res = 300)
hist(perm_aucs,
    main = "Permutation Test AUC Distribution",
    xlab = "AUC", breaks = 20, col = "lightblue",
    xlim = c(0, 1), ylim = c(0, 30),
    cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.5
)
abline(v = actual_auc, col = "red", lwd = 2)
text(actual_auc - 0.14, 25, paste("Actual AUC:", round(actual_auc, 2)),
    cex = 1, font = 3, col = "red"
)
dev.off()

## best model
best_model <- train(
    x = train[, ranked_features[1:151]], y = train$Allgr, method = "glmnet", metric = "ROC",
    trControl = control, preProcess = c("zv", "center", "scale", "medianImpute")
)

Imp <- varImp(best_model)$importance

Imp$Metabolite <- rownames(Imp)

top_metabolites <- Imp |>
    dplyr::arrange(desc(Overall)) |>
    dplyr::rename(Importance = Overall) |>
    tibble::remove_rownames() |>
    dplyr::slice_head(n = 25)

ggplot(top_metabolites, aes(x = Importance, y = forcats::fct_reorder(Metabolite, Importance))) +
    geom_bar(stat = "identity", fill = "steelblue", width = 0.6) +
    # geom_text(aes(label = round(Importance, 2)), hjust = -0.2, size = 5, fontface = "bold") +
    labs(title = "Top 25 Metabolites (LASSO Importance)", x = "Importance Score", y = "Metabolite") +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

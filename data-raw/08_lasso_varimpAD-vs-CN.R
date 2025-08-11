## code to prepare `lasso_varimpAD vs CN` dataset goes here
####### Modelling for AD vs CN
library(caret)
library(readr)
data3_1 <- read_csv("./vignettes/03_EffectsAdjusted.csv")
data3_1 <- as.data.frame(data3_1) ## data involving only significant biomarkers based on q values of group

CN <- data3_1[which(data3_1$Allgr == "CN"), ]
AD <- data3_1[which(data3_1$Allgr == "AD"), ]

train_nCN <- sample(nrow(CN), nrow(CN) * 0.7)
train_nAD <- sample(nrow(AD), nrow(AD) * 0.7)

train_CN <- CN[train_nCN, ]
train_AD <- AD[train_nAD, ]
train <- rbind(train_CN, train_AD)

test_verif_CN <- CN[-train_nCN, ]
test_verif_AD <- AD[-train_nAD, ]

test <- rbind(test_verif_CN, test_verif_AD)

dat <- train

dat$Allgr <- as.factor(dat$Allgr)
n <- ncol(CN)
n_cv <- 5
n_repeat <- 20
n_param_max <- 400

my_folds <- createMultiFolds(dat$Allgr, k = n_cv, times = n_repeat)
seeds <- list()
# one seed for each tuning parameter at each resampling (can be longer than the number of parameters to tune)
for (i in seq_len(n_cv * n_repeat)) {
    seeds[[i]] <- sample.int(n = 1000, n_param_max)
}
# for the last model
seeds[[n_cv * n_repeat + 1]] <- sample.int(1000, 1)

control <- trainControl(
    "repeatedcv",
    index = my_folds,
    selectionFunction = "best",
    classProbs = TRUE,
    savePredictions = "final",
    summaryFunction = twoClassSummary,
    seeds = seeds, verboseIter = T
)
# control <- trainControl(
#     method = "LOOCV",
#     classProbs = TRUE,
#     savePredictions = "final",
#     summaryFunction = twoClassSummary,
#     verboseIter = TRUE
# )
dir.create(paste0(getwd(), "/export"), showWarnings = F)
n <- ncol(dat)

library(doParallel)
cores <- detectCores() - 5
cl <- makeCluster(cores)
registerDoParallel(cl)
lasso_grid <- train(
    x = dat[, c(8:n)],
    y = dat$Allgr,
    method = "glmnet",
    metric = "ROC",
    trControl = control,
    tuneLength = 40,
    preProcess = c("zv", "center", "scale", "medianImpute")
)

stopCluster(cl)
# unregister_dopar()

saveRDS(lasso_grid, file = paste0(getwd(), "/export/lasso_grid.Rds"), compress = F)

unlink(paste0(getwd(), "/export/lasso_models"), recursive = T)
dir.create(paste0(getwd(), "/export/lasso_models"))

lasso_funcs <- getModelInfo("glmnet", regex = FALSE)[[1]]

lasso_funcs$fit <- eval(str2lang(sub(
    "(.*?)\\{(.*)",
    '\\1 {
    on.exit({
      if (!last) {
        saveRDS(out,
          file = paste0("export/lasso_models/lasso_", gsub(" |:", "_",
            Sys.time()), "_", sample(10**6, 1), ".Rds"),
          compress = F)
      }
    }, add = T)
    \\2', deparse1(lasso_funcs$fit, collapse = "\n")
)))

n <- ncol(dat)
cores <- detectCores() - 5
cl <- makeCluster(cores)
registerDoParallel(cl)
lasso <- train(
    x = dat[, c(8:ncol(dat))],
    y = dat$Allgr,
    method = lasso_funcs, tuneGrid = lasso_grid$bestTune,
    metric = "ROC",
    tuneLength = 1,
    trControl = control,
    preProcess = c("zv", "center", "scale", "medianImpute")
)
stopCluster(cl)

cv_models <- list.files("export/lasso_models/", full.names = T)
# if (length(cv_models) != n_cv * n_repeat) {
#     stop("Different number of lasso models than folds.")
# }

n <- length(cv_models)
# lasso_varimp <- colnames(CN)[8:ncol(CN)]
# for (i in 1:n) {
#     my_lasso <- readRDS(cv_models[[i]])
#     lasso_varimp <- cbind(lasso_varimp,varImp(my_lasso))
# }

varimp_list <- list()

for (i in 1:n) {
    my_lasso <- readRDS(cv_models[[i]])
    lasso_imp_df <- varImp(my_lasso) %>%
        as.data.frame() %>%
        rownames_to_column(var = "variable") %>%
        rename(!!paste0("Model_", i) := Overall)

    varimp_list[[i]] <- lasso_imp_df
}

# Merge all importance data frames by variable name
lasso_varimp <- purrr::reduce(varimp_list, full_join, by = "variable")


varimp_binary <- apply(lasso_varimp[, -1], 2, function(x) ifelse(x > 0, 1, 0))
lasso_VIS_ADvsCN <- NULL
lasso_VIS_ADvsCN <- cbind(
    lasso_VIS_ADvsCN,
    cbind(
        Metabolites = lasso_varimp[, 1],
        mean_VIS = rowMeans(lasso_varimp[, -1]),
        numberoftimes = rowSums(varimp_binary)
    )
)
write.csv(lasso_VIS_ADvsCN, "lasso_varimpAD vs CN.csv", row.names = F)

usethis::use_data(lasso_VIS_ADvsCN, overwrite = TRUE)


lasso_coef <- coef(lasso_grid$finalModel, s = lasso_grid$bestTune$lambda)
coef_df <- as.data.frame(as.matrix(lasso_coef))
coef_df$Feature <- rownames(coef_df)
colnames(coef_df) <- c("Coefficient", "Feature")

coef_df <- coef_df[coef_df$Coefficient != 0 & coef_df$Feature != "(Intercept)", ]

top_features <- coef_df[order(abs(coef_df$Coefficient), decreasing = TRUE), ][1:50, ]

ggplot(top_features, aes(
    x = reorder(Feature, abs(Coefficient)),
    y = Coefficient, fill = Coefficient > 0
)) +
    coord_flip() +
    geom_bar(stat = "identity", show.legend = FALSE) +
    scale_fill_manual(values = c("#D55E00", "#0072B2")) +
    labs(
        title = "Top 50 LASSO features (AD vs CN)",
        subtitle = "Ranked by absolute coefficient magnitude",
        x = "Metabolites features", y = "Coefficient value",
        panel.grid.major.y = element_blank()
    ) +
    theme_hc()

ggsave("top50_lasso_features_ADvsCN.png", width = 6, height = 7.5, dpi = 300)

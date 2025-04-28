library(tidyverse)
library(caret)
library(randomForest)
# library(warnings)
library(VIM)
library(class) #for kNN
library(FNN)
library(caret) #for train-test split,  logistic regression and confusion matrix
library(ggplot2)    #for visualization
library(gridExtra)  #for arranging multiple plots

options(warn = -1)

set.seed(1234)

adni_df <- read.csv('data-raw/ADNIMERGE_20Sep2024.csv')

## data processing ====

# "Dementia" -> "AD"
# With the exception of 6 cases, all observations for which DX=="Dementia" had
# baseline diagnoses of "AD", suggesting that Dementia can be recoded to AD
erronious_obs <- adni_df |>
    dplyr::filter(DX == "Dementia", VISCODE == "bl",
                DX_bl != "AD")

adni_df <- adni_df |>
    dplyr::filter(!(row_number() %in% row_number(erronious_obs))) %>%
    dplyr::mutate(DX = ifelse(DX == "Dementia", "AD", DX))

# collapse DX_bl groups: "EMCI"/"LMCI" --> "MCI"; "SMC" --> "CN"
adni_df <- adni_df |>
    dplyr::mutate(DX_bl = case_when(
        DX_bl %in% c("EMCI", "LMCI") ~ "MCI",
        DX_bl == "SMC" ~ "CN",
        TRUE ~ DX_bl
    ))

# identifying erroneous observations for baseline groups
erronious_obs_AD_MCI <- adni_df |>
    dplyr::filter(DX_bl == "AD", DX == "MCI", VISCODE == "bl")

erronious_obs_MCI_CN <- adni_df |>
    dplyr::filter(DX_bl == "MCI", DX == "CN", VISCODE == "bl")

erronious_obs <- bind_rows(erronious_obs_AD_MCI,
                           erronious_obs_MCI_CN,
                           erronious_obs)

adni_df <- adni_df |>
    dplyr::filter(!(row_number() %in% row_number(erronious_obs_AD_MCI))) %>%
    dplyr::filter(!(row_number() %in% row_number(erronious_obs_MCI_CN)))

# splitting dataset by original protocol
adni1_df <- adni_df %>% filter(ORIGPROT == "ADNI1")
adni2_df <- adni_df %>% filter(ORIGPROT == "ADNI2")
adni3_df <- adni_df %>% filter(ORIGPROT == "ADNI3")
adnigo_df <- adni_df %>% filter(ORIGPROT == "ADNIGO")


# description of the dataset ====
description <- "
Subsetting our Data
We use the ADNI 1 phase for our cross-sectional models. The justification
for using one of the four phases is that we've noticed that there are certain
predictors that are entirely missing for certain phases of the study.
For example, Everyday Cognition (Ecog) is missing in general for the ADNI 1 phase.
Furthermore, there are differences in the ADNI phases with regard to the amount
of information collected and the study procedure followed.
This means that patient records in ADNI 1 are not directly comparable
with patient records in ADNI GO. In addition, our baseline model aims
to predict baseline diagnosis, so we only choose records where the patient
was visiting for a baseline test. In addition, we drop predictors
where the number of missing values is greater than or equal to 90% of
the total number of records. In addition, we drop obviously collinear
features (such as Month and Year after baseline).
Finally, we create dummies (one-hot-encoding) for our categorical variables
(such as marital status, Race etc).

Dealing with Missing Data
Our EDA also showedd us that certain predictors, such as the ones taken from
the Cerebrospinal Fluid (CSF) such as Amyloid Beta (ABeta), Tau and Ptau, are
missing in large chunks for participants who did not have this procedure done.
This means that we could not simply follow a mean or median value imputation
since the patients who displayed only minor symptoms of cognitive decline are
systemically less likely to have their CSF taken than patients who display a
greater degree of cognitive decline. In order to deal with these concerns, we
use k-Nearest Neighbors (kNN) with 5 neighbors to impute values that are
missing. This appears to be the fairest system of imputing missing values.
"

## model 1: predicting baseline diagnosis using all baseline features ====
description <- "We began by developing a basic multiclass logistic regression
model to predict baseline diagnosis from data available at that initial visit
to get a rough understanding of the strength of the available predictors,
as well as to compare how well they perform in predicting cross-sectional
diagnosis versus longitudinal change in diagnosis.

The one-versus-all logistic regression classifier achieves a training score
of 80.6% and a test score of 78.3% (based on a 25% train-test split stratified
by diagnosis). The confusion matrices for the training and test sets show that
the model does a good job discerning between the three classes, without always
defaulting to -- though clearly biased toward -- the dominant group. (48.8% of
observations were classified as MCI, 28.2% as CN, and 23.0% as AD).
"

# ADNI1 phase for cross-sectional models
adni1_df <- adni1_df |>
    dplyr::filter(M == 0) # only keeping baseline measures

# dealing with '<' and '>' in ABETA, TAU, and PTAU
genes <- c("TAU", "ABETA", "PTAU")
for (val in genes) {
    all_vals <- adni1_df[[val]]
    adni1_df <- adni1_df |>
        dplyr:::select(-all_of(val)) |>
        dplyr::mutate(!!val := as.numeric(gsub("[<>]", "", all_vals)))
}

# drop collinear, unrelated, or missing variables
collinear <- c("DX_bl", "M", "VISCODE", "update_stamp", "Month")
unrelated <- c("RID", "PTID", "SITE", "ORIGPROT", "COLPROT",
               "FSVERSION", "IMAGEUID", "FLDSTRENG", "EXAMDATE")

# identify baseline variables and those with high missingness
nulls <- colSums(is.na(adni1_df))
baseline_vars <- names(adni1_df)[grepl("_bl|_BL", names(adni1_df))]
null_cols <- names(adni1_df)[nulls > nrow(adni1_df) / 1.1]

# drop unnecessary variables
adni1_df <- adni1_df |>
    dplyr::select(-all_of(c(collinear, unrelated,
                            baseline_vars, null_cols)))

# dealing with categorical variables
adni1_df <- adni1_df |>
    dplyr::mutate(PTGENDER_Female = ifelse(PTGENDER == "Female", 1, 0)) |>
    dplyr::select(-PTGENDER)

categoricals <- c("PTETHCAT", "PTRACCAT", "PTMARRY")
dummies <- model.matrix(~ . - 1, data = adni1_df[categoricals]) |>
    as.data.frame()

adni1_df <- adni1_df |>
    dplyr::bind_cols(dummies) |>
    dplyr::select(-all_of(categoricals))

# drop rows missing the response variable (DX)
adni1_df <- adni1_df |>
    drop_na(DX)

# make the 'DX' variable numeric
adni1_df <- adni1_df |>
    mutate(DX = case_when(
        DX == "CN" ~ 1,
        DX == "MCI" ~ 2,
        DX == "AD" ~ 3,
        TRUE ~ NA_real_
    ))

df <- adni1_df

# kNN Imputation with k = 5
weird_value <- -247
df[is.na(df)] <- weird_value  # Replace NAs with a placeholder value

nulls <- colnames(df)[colSums(df == weird_value) > 0]
ind_dict <- list()

for (col in nulls) {
    pseudo_train <- df[df[[col]] != weird_value, ]
    pseudo_x <- pseudo_train[, setdiff(names(df), c(col, "DX"))]
    pseudo_y <- pseudo_train[[col]]

    pseudo_test <- df[df[[col]] == weird_value, ]
    pseudo_x_tes <- pseudo_test[, setdiff(names(df), c(col, "DX"))]

    # kNN regression
    pred <- FNN::knn.reg(train = pseudo_x, y = pseudo_y, test = pseudo_x_tes, k = 5)

    # predictions for imputation
    ind_dict[[col]] <- pred$pred
}

# actual imputation
for (col in nulls) {
    indices <- which(df[[col]] == weird_value)
    df[indices, col] <- ind_dict[[col]]
}

# categorizing variables
demographic <- c('AGE', 'PTGENDER_Female', 'PTEDUCAT')
ravlt <- c('RAVLT_immediate', 'RAVLT_learning', 'RAVLT_forgetting', 'RAVLT_perc_forgetting')
adas <- c('ADAS11', 'ADAS13', 'ADASQ4')
exams <- c('MMSE', 'CDRSB', 'TRABSCOR', 'FAQ', 'mPACCdigit', 'mPACCtrailsB', 'LDELTOTAL')
genetic <- c('APOE4', 'TAU', 'ABETA', 'PTAU', 'FDG')
brain <- c('Ventricles', 'Hippocampus', 'WholeBrain', 'Entorhinal', 'Fusiform', 'MidTemp', 'ICV')
adni1 <- c('DIGITSCOR')

adni2 <- c('AV45', 'MOCA',
           'EcogPtMem', 'EcogPtLang', 'EcogPtVisspat', 'EcogPtPlan', 'EcogPtOrgan', 'EcogPtDivatt',
           'EcogPtTotal', 'EcogSPMem', 'EcogSPLang', 'EcogSPVisspat', 'EcogSPPlan', 'EcogSPOrgan', 'EcogSPTotal',
           'PTRACCAT_Unknown', 'PTRACCAT_Hawaiian/Other PI')

response <- c('DX')

# Train-test split
set.seed(209)
train_index <- caret::createDataPartition(df$DX, p = 0.75, list = FALSE)
train <- df[train_index, ]
test <- df[-train_index, ]

X_train <- train[, setdiff(names(train), "DX")]
X_test <- test[, setdiff(names(test), "DX")]
y_train <- train$DX
y_test <- test$DX

# Logistic regression
baseline_logreg <- train(
    DX ~ .,
    data = cbind(X_train, DX = y_train),
    method = "svmLinear",
    trControl = trainControl(method = "none"),
    maxit = 1000
)

# predictions and confusion matrices
train_preds <- stats::predict(baseline_logreg, newdata = X_train, type = "raw")
test_preds <- stats::predict(baseline_logreg, newdata = X_test)


# assigning class based on closest integer or threshold
# as the predictions represent probabilities or scores, I considered
# using a threshold-based classification instead of rounding
train_preds_class <- factor(ifelse(train_preds < 1.5, 1,
                                   ifelse(train_preds < 2.5, 2, 3)))
cnf_matrix_tr <- confusionMatrix(train_preds_class, factor(y_train))

test_preds_class <- factor(ifelse(test_preds < 1.5, 1,
                                   ifelse(train_preds < 2.5, 2, 3)))
cnf_matrix_ts <- confusionMatrix(test_preds_class, factor(y_test))

plot_confusion_matrix <- function(conf_matrix, title, classes, normalize = FALSE) {
    df <- as.data.frame(conf_matrix$table)
    if (normalize) {
        df$Freq <- df$Freq / rowSums(conf_matrix$table)[df$Reference]
    }
    df$Prop <- df$Freq / rowSums(conf_matrix$table)[df$Prediction]
    ggplot(df, aes(x = Prediction, y = Reference, fill = Freq)) +
        geom_tile(color = "white") +
        geom_text(aes(label = sprintf(ifelse(normalize, "%.2f", "%d"), Freq)),
                  color = "black", size = 6) +
        scale_fill_gradient2(low = "white", high = "blue", midpoint = max(df$Freq) / 2) +
        labs(title = title, x = "Predicted Label", y = "True Label", fill = "Proportion") +
        scale_x_discrete(labels = classes) +
        scale_y_discrete(labels = classes) +
        theme_minimal(base_size = 15) +
        theme(axis.text.x = element_text(size = 13, angle = 0),
              axis.text.y = element_text(size = 13),
              plot.title = element_text(size = 18, face = "bold"))
}

class_labels <- c("CN", "MCI", "AD")

plot_train <- plot_confusion_matrix(cnf_matrix_tr,
                                    title = "Training Set",
                                    class_labels,
                                    normalize = F)
plot_test <- plot_confusion_matrix(cnf_matrix_ts,
                                   title = "Test Set",
                                   class_labels,
                                   normalize = F)

# Arrange plots side by side
grid.arrange(plot_train, plot_test, nrow = 1)


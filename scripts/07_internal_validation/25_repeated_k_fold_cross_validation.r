# Repeated K-fold cross-validation for the elastic-net age prediction model
# This repeats K-fold cross-validation multiple times to get more robust performance estimates
library(glmnet)
library(rsample)
source("scripts/common/elastic_net_alpha_tuning.r")

dir.create("results/internal_validation", recursive = TRUE, showWarnings = FALSE)

set.seed(123)

beta_matrix <- readRDS("data/GSE87571/beta_matrix_age_model.rds")
metadata <- read.csv("data/GSE87571/modelling_metadata_age_model.csv")

# glmnet expects samples as rows and CpG sites as columns
x <- t(beta_matrix[, match(metadata$sample_id, colnames(beta_matrix))])

metadata_folds <- vfold_cv(metadata, v = 10, repeats = 5, strata = age)

alpha_grid <- seq(0.05, 1, by = 0.05)

all_performance <- data.frame()
all_residuals <- data.frame()
all_selected_cpgs <- data.frame()
all_alpha_tuning <- data.frame()

for (i in seq_len(nrow(metadata_folds))) {
  train_metadata <- analysis(metadata_folds$splits[[i]])
  test_metadata <- assessment(metadata_folds$splits[[i]])

  # Match beta matrix columns to metadata rows for training and test sets
  x_train <- x[train_metadata$sample_id, ]
  y_train <- train_metadata$age

  x_test <- x[test_metadata$sample_id, ]
  y_test <- test_metadata$age 

  # Tune alpha and lambda using the training folds only
  alpha_tuned_model <- tune_alpha_model(x_train, y_train, alpha_grid)
  k_fold_model <- alpha_tuned_model$model

  # Save the CpGs selected by this fold model
  fold_selected_cpgs <- get_selected_cpgs(
    k_fold_model,
    "repeated_k_fold_cross_validation",
    metadata_folds$id[i]
  )
  fold_selected_cpgs$selected_alpha <- alpha_tuned_model$selected_alpha
  fold_selected_cpgs$lambda_min <- k_fold_model$lambda.min
  fold_selected_cpgs$lambda_1se <- k_fold_model$lambda.1se

  fold_alpha_tuning <- alpha_tuned_model$alpha_performance
  fold_alpha_tuning$validation_method <- "repeated_k_fold_cross_validation"
  fold_alpha_tuning$resample_id <- metadata_folds$id[i]

  # Predict age in the held-out fold
  predicted_age <- predict(
    k_fold_model,
    newx = x_test,
    s = "lambda.min"
  )

  predicted_age <- as.numeric(predicted_age)

  fold_residuals <- data.frame(
    validation_method = "repeated_k_fold_cross_validation",
    resample_id = metadata_folds$id[i],
    sample_id = test_metadata$sample_id,
    geo_accession = test_metadata$geo_accession,
    age = y_test,
    predicted_age = predicted_age,
    residual = predicted_age - y_test,
    absolute_error = abs(predicted_age - y_test)
  )

  fold_performance <- data.frame(
    fold = metadata_folds$id[i],
    training_samples = length(y_train),
    test_samples = length(y_test),
    input_cpgs = ncol(x),
    selected_alpha = alpha_tuned_model$selected_alpha,
    lambda_min = k_fold_model$lambda.min,
    lambda_1se = k_fold_model$lambda.1se,
    selected_cpgs = sum(coef(k_fold_model, s = "lambda.min")[-1, ] != 0),
    mae = mean(abs(predicted_age - y_test)),
    median_absolute_error = median(abs(predicted_age - y_test)),
    rmse = sqrt(mean((predicted_age - y_test)^2)),
    mean_error = mean(predicted_age - y_test),
    correlation = cor(predicted_age, y_test),
    r_squared = cor(predicted_age, y_test)^2
  )

  all_performance <- rbind(all_performance, fold_performance)
  all_residuals <- rbind(all_residuals, fold_residuals)
  all_selected_cpgs <- rbind(all_selected_cpgs, fold_selected_cpgs)
  all_alpha_tuning <- rbind(all_alpha_tuning, fold_alpha_tuning)
}

repeated_k_fold_summary <- data.frame(
  folds = nrow(all_performance),
  input_cpgs = ncol(x),
  mean_selected_cpgs = mean(all_performance$selected_cpgs),
  sd_selected_cpgs = sd(all_performance$selected_cpgs),
  min_selected_cpgs = min(all_performance$selected_cpgs),
  max_selected_cpgs = max(all_performance$selected_cpgs),
  mean_mae = mean(all_performance$mae),
  sd_mae = sd(all_performance$mae),
  mean_median_absolute_error = mean(all_performance$median_absolute_error),
  mean_rmse = mean(all_performance$rmse),
  sd_rmse = sd(all_performance$rmse),
  mean_error = mean(all_performance$mean_error),
  mean_correlation = mean(all_performance$correlation),
  mean_r_squared = mean(all_performance$r_squared)
)

write.csv(
  all_performance,
  "results/internal_validation/repeated_k_fold_per_fold_summary.csv",
  row.names = FALSE
)

write.csv(
  all_residuals,
  "results/internal_validation/repeated_k_fold_residuals.csv",
  row.names = FALSE
)

write.csv(
  all_selected_cpgs,
  "results/internal_validation/repeated_k_fold_selected_cpgs.csv",
  row.names = FALSE
)

write.csv(
  all_alpha_tuning,
  "results/internal_validation/repeated_k_fold_alpha_tuning.csv",
  row.names = FALSE
)

write.csv(
  repeated_k_fold_summary,
  "results/internal_validation/repeated_k_fold_summary.csv",
  row.names = FALSE
)

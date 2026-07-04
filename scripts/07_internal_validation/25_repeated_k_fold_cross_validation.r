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

alpha_grid <- c(0.25, 0.50, 0.75)

array_task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (array_task_id != "") {
  array_task_id <- as.integer(array_task_id)
  repeat_name <- paste0("Repeat", array_task_id)
  fold_indices <- grep(paste0("^", repeat_name), metadata_folds$id)

  if (length(fold_indices) == 0) {
    folds_per_chunk <- 10
    first_fold <- ((array_task_id - 1) * folds_per_chunk) + 1
    last_fold <- min(array_task_id * folds_per_chunk, nrow(metadata_folds))
    fold_indices <- first_fold:last_fold
  }

  output_suffix <- paste0("_chunk_", sprintf("%02d", array_task_id))
} else {
  fold_indices <- seq_len(nrow(metadata_folds))
  output_suffix <- ""
}

all_performance <- data.frame()
all_residuals <- data.frame()
all_selected_cpgs <- data.frame()
all_alpha_tuning <- data.frame()

for (i in fold_indices) {
  fold_id <- metadata_folds$id[i]
  if ("id2" %in% names(metadata_folds)) {
    fold_id <- paste(metadata_folds$id[i], metadata_folds$id2[i], sep = "_")
  }

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
    fold_id
  )
  fold_selected_cpgs$selected_alpha <- alpha_tuned_model$selected_alpha
  fold_selected_cpgs$lambda_min <- k_fold_model$lambda.min
  fold_selected_cpgs$lambda_1se <- k_fold_model$lambda.1se

  fold_alpha_tuning <- alpha_tuned_model$alpha_performance
  fold_alpha_tuning$validation_method <- "repeated_k_fold_cross_validation"
  fold_alpha_tuning$resample_id <- fold_id

  # Predict age in the held-out fold
  predicted_age <- predict(
    k_fold_model,
    newx = x_test,
    s = "lambda.min"
  )

  predicted_age <- as.numeric(predicted_age)

  fold_residuals <- data.frame(
    validation_method = "repeated_k_fold_cross_validation",
    resample_id = fold_id,
    sample_id = test_metadata$sample_id,
    geo_accession = test_metadata$geo_accession,
    age = y_test,
    predicted_age = predicted_age,
    residual = predicted_age - y_test,
    absolute_error = abs(predicted_age - y_test)
  )

  fold_performance <- data.frame(
    fold = fold_id,
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
  mean_selected_alpha = mean(all_performance$selected_alpha),
  sd_selected_alpha = sd(all_performance$selected_alpha),
  min_selected_alpha = min(all_performance$selected_alpha),
  max_selected_alpha = max(all_performance$selected_alpha),
  mean_lambda_min = mean(all_performance$lambda_min),
  mean_lambda_1se = mean(all_performance$lambda_1se),
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
  paste0(
    "results/internal_validation/repeated_k_fold_per_fold_summary",
    output_suffix,
    ".csv"
  ),
  row.names = FALSE
)

write.csv(
  all_residuals,
  paste0(
    "results/internal_validation/repeated_k_fold_residuals",
    output_suffix,
    ".csv"
  ),
  row.names = FALSE
)

write.csv(
  all_selected_cpgs,
  paste0(
    "results/internal_validation/repeated_k_fold_selected_cpgs",
    output_suffix,
    ".csv"
  ),
  row.names = FALSE
)

write.csv(
  all_alpha_tuning,
  paste0(
    "results/internal_validation/repeated_k_fold_alpha_tuning",
    output_suffix,
    ".csv"
  ),
  row.names = FALSE
)

write.csv(
  repeated_k_fold_summary,
  paste0(
    "results/internal_validation/repeated_k_fold_summary",
    output_suffix,
    ".csv"
  ),
  row.names = FALSE
)

# K-fold cross-validation for the elastic-net age prediction model
# This trains and tests the model across 10 age-stratified folds

library(glmnet)
library(rsample)

dir.create("results/internal_validation", recursive = TRUE, showWarnings = FALSE)

set.seed(123)

beta_matrix <- readRDS("data/GSE87571/beta_matrix_age_model.rds")
metadata <- read.csv("data/GSE87571/modelling_metadata_age_model.csv")

# glmnet expects samples as rows and CpG sites as columns
x <- t(beta_matrix[, match(metadata$sample_id, colnames(beta_matrix))])

# Create 10 folds, stratified by age
# v = 10 is the number of folds, and strata = age ensures similar age distributions in each fold
metadata_folds <- vfold_cv(metadata, v = 10, strata = age)

all_performance <- data.frame()
all_residuals <- data.frame()

# loop through each fold
for (i in seq_len(nrow(metadata_folds))) {
  train_metadata <- analysis(metadata_folds$splits[[i]])
  test_metadata <- assessment(metadata_folds$splits[[i]])

    # Match beta matrix columns to metadata rows for training and test sets
  x_train <- x[train_metadata$sample_id, ]
  y_train <- train_metadata$age

  x_test <- x[test_metadata$sample_id, ]
  y_test <- test_metadata$age 

  # Train the elastic-net model using the training folds only
  k_fold_model <- cv.glmnet(
    x = x_train,
    y = y_train,
    alpha = 0.5,
    family = "gaussian"
  )

  # Predict age in the held-out fold
  predicted_age <- predict(
    k_fold_model,
    newx = x_test,
    s = "lambda.min"
  )

  predicted_age <- as.numeric(predicted_age)

  fold_residuals <- data.frame(
    validation_method = "k_fold_cross_validation",
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
}

k_fold_summary <- data.frame(
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
  "results/internal_validation/k_fold_per_fold_summary.csv",
  row.names = FALSE
)

write.csv(
  all_residuals,
  "results/internal_validation/k_fold_residuals.csv",
  row.names = FALSE
)

write.csv(
  k_fold_summary,
  "results/internal_validation/k_fold_summary.csv",
  row.names = FALSE
)

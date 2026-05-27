# K-fold cross-validation for the elastic-net age prediction model
# This trains and tests the model across 10 age-stratified folds

library(glmnet)
library(rsample)

set.seed(123)

beta_matrix <- readRDS("data/GSE87571/beta_matrix_age_model.rds")
metadata <- read.csv("data/GSE87571/modelling_metadata_age_model.csv")

# glmnet expects samples as rows and CpG sites as columns
x <- t(beta_matrix[, match(metadata$sample_id, colnames(beta_matrix))])

# Create 10 folds, stratified by age
# v = 10 is the number of folds, and strata = age ensures similar age distributions in each fold
metadata_folds <- vfold_cv(metadata, v = 10, strata = age)

all_predictions <- data.frame()
all_performance <- data.frame()

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

  fold_predictions <- data.frame(
    fold = metadata_folds$id[i], 
    sample_id = test_metadata$sample_id,
    geo_accession = test_metadata$geo_accession,
    age = y_test,
    predicted_age = predicted_age,
    age_error = predicted_age - y_test,
    absolute_error = abs(predicted_age - y_test)
  )

  fold_performance <- data.frame(
    fold = metadata_folds$id[i],
    training_samples = length(y_train),
    test_samples = length(y_test),
    cpgs = ncol(x),
    mae = mean(abs(predicted_age - y_test)),
    median_absolute_error = median(abs(predicted_age - y_test)),
    rmse = sqrt(mean((predicted_age - y_test)^2)),
    mean_error = mean(predicted_age - y_test),
    correlation = cor(predicted_age, y_test),
    r_squared = cor(predicted_age, y_test)^2
  )

  all_predictions <- rbind(all_predictions, fold_predictions)
  all_performance <- rbind(all_performance, fold_performance)
}

k_fold_summary <- data.frame(
  folds = nrow(all_performance),
  cpgs = ncol(x),
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
  all_predictions,
  "results/internal_validation/k_fold_predictions.csv",
  row.names = FALSE
)

write.csv(
  all_performance,
  "results/internal_validation/k_fold_per_fold_summary.csv",
  row.names = FALSE
)

write.csv(
  k_fold_summary,
  "results/internal_validation/k_fold_summary.csv",
  row.names = FALSE
)

pdf("results/internal_validation/k_fold_predicted_vs_actual_age.pdf")
plot(
  all_predictions$age,
  all_predictions$predicted_age,
  xlab = "Chronological age",
  ylab = "Predicted age",
  main = "K-fold cross-validation",
  pch = 16
)
abline(0, 1, col = "red")
dev.off()

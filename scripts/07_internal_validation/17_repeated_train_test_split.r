# Repeated train-test split validation for the elastic-net age prediction model
# This repeats random 80:20 train-test splits using Monte Carlo cross-validation

library(glmnet)
library(rsample)

set.seed(123)

beta_matrix <- readRDS("data/GSE87571/beta_matrix_age_model.rds")
metadata <- read.csv("data/GSE87571/modelling_metadata_age_model.csv")

# glmnet expects samples as rows and CpG sites as columns
x <- t(beta_matrix[, match(metadata$sample_id, colnames(beta_matrix))])

# Create repeated 80:20 train-test splits
# Stratify by age so training and test sets have similar age distributions
metadata_splits <- mc_cv(metadata, prop = 0.8, times = 10, strata = age)

# data frame to store performance metrics for each split
all_performance <- data.frame()

#loop through each split
for (i in seq_len(nrow(metadata_splits))) {
  train_metadata <- analysis(metadata_splits$splits[[i]])
  test_metadata <- assessment(metadata_splits$splits[[i]])

# Match beta matrix columns to metadata rows for training and test sets
  x_train <- x[train_metadata$sample_id, ]
  y_train <- train_metadata$age

  x_test <- x[test_metadata$sample_id, ]
  y_test <- test_metadata$age

  # Train the elastic-net model using the training samples only
  repeated_train_test_model <- cv.glmnet(
    x = x_train,
    y = y_train,
    alpha = 0.5,
    family = "gaussian"
  )

  # Predict age in the held-out test samples
  predicted_age <- predict(
    repeated_train_test_model,
    newx = x_test,
    s = "lambda.min"
  )

  predicted_age <- as.numeric(predicted_age)

# Calculate performance metrics for this split and store in the data frame
  split_performance <- data.frame(
    split = metadata_splits$id[i],
    training_samples = length(y_train),
    test_samples = length(y_test),
    input_cpgs = ncol(x),
    selected_cpgs = sum(coef(repeated_train_test_model, s = "lambda.min")[-1, ] != 0),
    mae = mean(abs(predicted_age - y_test)),
    median_absolute_error = median(abs(predicted_age - y_test)),
    rmse = sqrt(mean((predicted_age - y_test)^2)),
    mean_error = mean(predicted_age - y_test),
    correlation = cor(predicted_age, y_test),
    r_squared = cor(predicted_age, y_test)^2
  )

# Append this split's performance to the overall performance data frame
  all_performance <- rbind(all_performance, split_performance)
}

# Summarise performance across all splits
performance_summary <- data.frame(
  repeats = nrow(all_performance),
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
  "results/internal_validation/repeated_train_test_split_per_split_summary.csv",
  row.names = FALSE
)

write.csv(
  performance_summary,
  "results/internal_validation/repeated_train_test_split_summary.csv",
  row.names = FALSE
)

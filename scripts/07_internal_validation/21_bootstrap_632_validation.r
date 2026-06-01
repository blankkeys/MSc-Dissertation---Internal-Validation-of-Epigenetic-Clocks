# Bootstrap and .632 bootstrap validation for the elastic-net age prediction model
# This trains the model on bootstrap samples and tests it on out-of-bag samples

library(glmnet)

set.seed(123)

beta_matrix <- readRDS("data/GSE87571/beta_matrix_age_model.rds")
metadata <- read.csv("data/GSE87571/modelling_metadata_age_model.csv")

# glmnet expects samples as rows and CpG sites as columns
x <- t(beta_matrix[, match(metadata$sample_id, colnames(beta_matrix))])
y <- metadata$age

# Use 100 bootstrap resamples
# Efron and Tibshirani recommend 50 to 200 bootstrap samples for this approach
n_bootstrap <- 100
n_samples <- nrow(x)

bootstrap_performance <- data.frame()

for (i in seq_len(n_bootstrap)) {
  # Sample individuals with replacement to create the bootstrap training set
  bootstrap_index <- sample(seq_len(n_samples), size = n_samples, replace = TRUE)

  # Samples not selected are the out-of-bag test set
  oob_index <- setdiff(seq_len(n_samples), unique(bootstrap_index))

  x_bootstrap <- x[bootstrap_index, , drop = FALSE]
  y_bootstrap <- y[bootstrap_index]

  x_oob <- x[oob_index, , drop = FALSE]
  y_oob <- y[oob_index]

  # Train the elastic-net model on the bootstrap sample
  bootstrap_model <- cv.glmnet(
    x = x_bootstrap,
    y = y_bootstrap,
    alpha = 0.5,
    family = "gaussian"
  )

  # Apparent predictions are made on the bootstrap training sample
  apparent_predicted_age <- predict(
    bootstrap_model,
    newx = x_bootstrap,
    s = "lambda.min"
  )

  # Out-of-bag predictions are made on samples not used for training
  oob_predicted_age <- predict(
    bootstrap_model,
    newx = x_oob,
    s = "lambda.min"
  )

  #calucalating bootsrap .632 estimates of performance

  apparent_predicted_age <- as.numeric(apparent_predicted_age)
  oob_predicted_age <- as.numeric(oob_predicted_age)

  apparent_mae <- mean(abs(apparent_predicted_age - y_bootstrap))
  oob_mae <- mean(abs(oob_predicted_age - y_oob))

  apparent_median_absolute_error <- median(abs(apparent_predicted_age - y_bootstrap))
  oob_median_absolute_error <- median(abs(oob_predicted_age - y_oob))

  apparent_mse <- mean((apparent_predicted_age - y_bootstrap)^2)
  oob_mse <- mean((oob_predicted_age - y_oob)^2)

  apparent_rmse <- sqrt(apparent_mse)
  oob_rmse <- sqrt(oob_mse)

  # The .632 estimate combines apparent and out-of-bag error
  bootstrap_632_mae <- (0.368 * apparent_mae) + (0.632 * oob_mae)
  bootstrap_632_median_absolute_error <- (0.368 * apparent_median_absolute_error) +
    (0.632 * oob_median_absolute_error)
  bootstrap_632_rmse <- sqrt((0.368 * apparent_mse) + (0.632 * oob_mse))

  split_performance <- data.frame(
    bootstrap = i,
    training_samples = length(bootstrap_index),
    unique_training_samples = length(unique(bootstrap_index)),
    out_of_bag_samples = length(oob_index),
    input_cpgs = ncol(x),
    selected_cpgs = sum(coef(bootstrap_model, s = "lambda.min")[-1, ] != 0),
    apparent_mae = apparent_mae,
    oob_mae = oob_mae,
    bootstrap_632_mae = bootstrap_632_mae,
    apparent_median_absolute_error = apparent_median_absolute_error,
    oob_median_absolute_error = oob_median_absolute_error,
    bootstrap_632_median_absolute_error = bootstrap_632_median_absolute_error,
    apparent_rmse = apparent_rmse,
    oob_rmse = oob_rmse,
    bootstrap_632_rmse = bootstrap_632_rmse,
    oob_mean_error = mean(oob_predicted_age - y_oob),
    oob_correlation = cor(oob_predicted_age, y_oob),
    oob_r_squared = cor(oob_predicted_age, y_oob)^2
  )

  bootstrap_performance <- rbind(bootstrap_performance, split_performance)
}

bootstrap_summary <- data.frame(
  bootstrap_resamples = nrow(bootstrap_performance),
  input_cpgs = ncol(x),
  mean_selected_cpgs = mean(bootstrap_performance$selected_cpgs),
  sd_selected_cpgs = sd(bootstrap_performance$selected_cpgs),
  min_selected_cpgs = min(bootstrap_performance$selected_cpgs),
  max_selected_cpgs = max(bootstrap_performance$selected_cpgs),
  mean_unique_training_samples = mean(bootstrap_performance$unique_training_samples),
  mean_out_of_bag_samples = mean(bootstrap_performance$out_of_bag_samples),
  mean_apparent_mae = mean(bootstrap_performance$apparent_mae),
  mean_oob_mae = mean(bootstrap_performance$oob_mae),
  mean_bootstrap_632_mae = mean(bootstrap_performance$bootstrap_632_mae),
  sd_bootstrap_632_mae = sd(bootstrap_performance$bootstrap_632_mae),
  mean_apparent_median_absolute_error = mean(
    bootstrap_performance$apparent_median_absolute_error
  ),
  mean_oob_median_absolute_error = mean(
    bootstrap_performance$oob_median_absolute_error
  ),
  mean_bootstrap_632_median_absolute_error = mean(
    bootstrap_performance$bootstrap_632_median_absolute_error
  ),
  mean_apparent_rmse = mean(bootstrap_performance$apparent_rmse),
  mean_oob_rmse = mean(bootstrap_performance$oob_rmse),
  mean_bootstrap_632_rmse = mean(bootstrap_performance$bootstrap_632_rmse),
  sd_bootstrap_632_rmse = sd(bootstrap_performance$bootstrap_632_rmse),
  mean_oob_error = mean(bootstrap_performance$oob_mean_error),
  mean_oob_correlation = mean(bootstrap_performance$oob_correlation),
  mean_oob_r_squared = mean(bootstrap_performance$oob_r_squared)
)

write.csv(
  bootstrap_performance,
  "results/internal_validation/bootstrap_632_per_resample_summary.csv",
  row.names = FALSE
)

write.csv(
  bootstrap_summary,
  "results/internal_validation/bootstrap_632_summary.csv",
  row.names = FALSE
)

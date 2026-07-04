# Bootstrap and .632 bootstrap validation for the elastic-net age prediction model
# This trains the model on bootstrap samples and tests it on out-of-bag samples

library(glmnet)
source("scripts/common/elastic_net_alpha_tuning.r")

dir.create("results/internal_validation", recursive = TRUE, showWarnings = FALSE)

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
alpha_grid <- c(0.25, 0.50, 0.75)

array_task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (array_task_id != "") {
  array_task_id <- as.integer(array_task_id)
  bootstraps_per_chunk <- 10
  first_bootstrap <- ((array_task_id - 1) * bootstraps_per_chunk) + 1
  last_bootstrap <- min(array_task_id * bootstraps_per_chunk, n_bootstrap)
  bootstrap_indices <- first_bootstrap:last_bootstrap
  output_suffix <- paste0("_chunk_", sprintf("%02d", array_task_id))
} else {
  bootstrap_indices <- seq_len(n_bootstrap)
  output_suffix <- ""
}

bootstrap_performance <- data.frame()

bootstrap_oob_residuals <- data.frame()
bootstrap_selected_cpgs <- data.frame()
bootstrap_alpha_tuning <- data.frame()

# bootstrap training set and out of bag test set
# out of bag is samples not used for training (on average 36.8% of samples)
for (i in bootstrap_indices) {
  set.seed(123 + i)

  # Sample individuals with replacement to create the bootstrap training set
  # randomply samples sample numbers from the dataset
  # with replacement means the same sample can be selected multiple times, 
  # and some samples may not be selected at all
  # eg , have 1,2,3,4,5,6,7,8,9,10
  # bootstrap sample could be 2,5,5,7,1,3,2,8,9,4
  # out of bag samples would be 6 and 10
  bootstrap_index <- sample(seq_len(n_samples), size = n_samples, replace = TRUE)

  # Samples not selected are the out-of-bag test set
  oob_index <- setdiff(seq_len(n_samples), unique(bootstrap_index))

  # if no samples are left for out-of-bag testing, skip this iteration
  if (length(oob_index) == 0) {
    next
  }

  # Create the bootstrap training methylation matrix for  bootstrap sample 
  x_bootstrap <- x[bootstrap_index, , drop = FALSE]
  # Create the bootstrap training age vector for bootstrap sample 
  # training age vector is the age of the samples in the bootstrap training set
  y_bootstrap <- y[bootstrap_index]

  # Create the out-of-bag test methylation matrix for out-of-bag samples 
  #and age vector for out-of-bag samples
  x_oob <- x[oob_index, , drop = FALSE]
  y_oob <- y[oob_index]

  # Tune alpha and lambda using the bootstrap training sample
  alpha_tuned_model <- tune_alpha_model(x_bootstrap, y_bootstrap, alpha_grid)
  bootstrap_model <- alpha_tuned_model$model

  # Save the CpGs selected by this bootstrap model
  resample_selected_cpgs <- get_selected_cpgs(
    bootstrap_model,
    "bootstrap",
    paste0("bootstrap_", i)
  )
  resample_selected_cpgs$selected_alpha <- alpha_tuned_model$selected_alpha
  resample_selected_cpgs$lambda_min <- bootstrap_model$lambda.min
  resample_selected_cpgs$lambda_1se <- bootstrap_model$lambda.1se

  resample_alpha_tuning <- alpha_tuned_model$alpha_performance
  resample_alpha_tuning$validation_method <- "bootstrap"
  resample_alpha_tuning$resample_id <- paste0("bootstrap_", i)

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

  oob_residuals <- data.frame(
    validation_method = "bootstrap_oob",
    # paste0 is used to create a unique identifier for each bootstrap resample
    resample_id = paste0("bootstrap_", i),
    sample_id = metadata$sample_id[oob_index],
    geo_accession = metadata$geo_accession[oob_index],
    age = y_oob,
    predicted_age = oob_predicted_age,
    residual = oob_predicted_age - y_oob,
    absolute_error = abs(oob_predicted_age - y_oob)
  )

  # Calculate performance metrics for apparent and out-of-bag predictions
  apparent_mae <- mean(abs(apparent_predicted_age - y_bootstrap))
  oob_mae <- mean(abs(oob_predicted_age - y_oob))

  # Median absolute error is less sensitive to outliers than mean absolute error, 
  # so we calculate both
  apparent_median_absolute_error <- median(abs(apparent_predicted_age - y_bootstrap))
  oob_median_absolute_error <- median(abs(oob_predicted_age - y_oob))

  # Mean squared error and root mean squared error are also calculated to provide a 
  # more comprehensive assessment of model performance
  apparent_mse <- mean((apparent_predicted_age - y_bootstrap)^2)
  oob_mse <- mean((oob_predicted_age - y_oob)^2)

  # Root mean squared error is on the same scale as the original age variable, 
  # making it easier to interpret than mean squared error
  apparent_rmse <- sqrt(apparent_mse)
  oob_rmse <- sqrt(oob_mse)

  apparent_mean_error <- mean(apparent_predicted_age - y_bootstrap)
  oob_mean_error <- mean(oob_predicted_age - y_oob)

  # The .632 estimate combines apparent and out-of-bag error
  bootstrap_632_mae <- (0.368 * apparent_mae) + (0.632 * oob_mae)
  bootstrap_632_median_absolute_error <- (0.368 * apparent_median_absolute_error) +
    (0.632 * oob_median_absolute_error)
  bootstrap_632_rmse <- sqrt((0.368 * apparent_mse) + (0.632 * oob_mse))
  bootstrap_632_mean_error <- (0.368 * apparent_mean_error) + (0.632 * oob_mean_error)

  # Summarise performance for this bootstrap resample
  split_performance <- data.frame(
    bootstrap = i,
    training_samples = length(bootstrap_index),
    unique_training_samples = length(unique(bootstrap_index)),
    out_of_bag_samples = length(oob_index),
    input_cpgs = ncol(x),
    selected_alpha = alpha_tuned_model$selected_alpha,
    lambda_min = bootstrap_model$lambda.min,
    lambda_1se = bootstrap_model$lambda.1se,
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
    apparent_mean_error = apparent_mean_error,
    oob_mean_error = oob_mean_error,
    bootstrap_632_mean_error = bootstrap_632_mean_error,
    apparent_correlation = cor(apparent_predicted_age, y_bootstrap),
    oob_correlation = cor(oob_predicted_age, y_oob),
    apparent_r_squared = cor(apparent_predicted_age, y_bootstrap)^2,
    oob_r_squared = cor(oob_predicted_age, y_oob)^2
  )

  bootstrap_performance <- rbind(bootstrap_performance, split_performance)
  bootstrap_oob_residuals <- rbind(bootstrap_oob_residuals, oob_residuals)
  bootstrap_selected_cpgs <- rbind(
    bootstrap_selected_cpgs,
    resample_selected_cpgs
  )
  bootstrap_alpha_tuning <- rbind(
    bootstrap_alpha_tuning,
    resample_alpha_tuning
  )
}

# Summarise bootstrap performance across all resamples
bootstrap_summary <- data.frame(
  bootstrap_resamples = nrow(bootstrap_performance),
  input_cpgs = ncol(x),
  mean_selected_cpgs = mean(bootstrap_performance$selected_cpgs),
  sd_selected_cpgs = sd(bootstrap_performance$selected_cpgs),
  min_selected_cpgs = min(bootstrap_performance$selected_cpgs),
  max_selected_cpgs = max(bootstrap_performance$selected_cpgs),
  mean_selected_alpha = mean(bootstrap_performance$selected_alpha),
  sd_selected_alpha = sd(bootstrap_performance$selected_alpha),
  min_selected_alpha = min(bootstrap_performance$selected_alpha),
  max_selected_alpha = max(bootstrap_performance$selected_alpha),
  mean_lambda_min = mean(bootstrap_performance$lambda_min),
  mean_lambda_1se = mean(bootstrap_performance$lambda_1se),
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
  mean_apparent_error = mean(bootstrap_performance$apparent_mean_error),
  mean_oob_error = mean(bootstrap_performance$oob_mean_error),
  mean_bootstrap_632_error = mean(bootstrap_performance$bootstrap_632_mean_error),
  mean_apparent_correlation = mean(bootstrap_performance$apparent_correlation),
  mean_oob_correlation = mean(bootstrap_performance$oob_correlation),
  mean_apparent_r_squared = mean(bootstrap_performance$apparent_r_squared),
  mean_oob_r_squared = mean(bootstrap_performance$oob_r_squared)
)

write.csv(
  bootstrap_performance,
  paste0(
    "results/internal_validation/bootstrap_632_per_resample_summary",
    output_suffix,
    ".csv"
  ),
  row.names = FALSE
)

write.csv(
  bootstrap_oob_residuals,
  paste0(
    "results/internal_validation/bootstrap_oob_residuals",
    output_suffix,
    ".csv"
  ),
  row.names = FALSE
)

write.csv(
  bootstrap_selected_cpgs,
  paste0(
    "results/internal_validation/bootstrap_selected_cpgs",
    output_suffix,
    ".csv"
  ),
  row.names = FALSE
)

write.csv(
  bootstrap_alpha_tuning,
  paste0(
    "results/internal_validation/bootstrap_alpha_tuning",
    output_suffix,
    ".csv"
  ),
  row.names = FALSE
)

write.csv(
  bootstrap_summary,
  paste0(
    "results/internal_validation/bootstrap_632_summary",
    output_suffix,
    ".csv"
  ),
  row.names = FALSE
)

# Combine bootstrap .632 validation chunk outputs

source("scripts/common/combine_validation_chunks.r")

dir.create("results/internal_validation", recursive = TRUE, showWarnings = FALSE)

bootstrap_performance <- read_validation_chunks(
  "^bootstrap_632_per_resample_summary_chunk_[0-9]+\\.csv$",
  10
)
bootstrap_oob_residuals <- read_validation_chunks(
  "^bootstrap_oob_residuals_chunk_[0-9]+\\.csv$",
  10
)
bootstrap_selected_cpgs <- read_validation_chunks(
  "^bootstrap_selected_cpgs_chunk_[0-9]+\\.csv$",
  10
)
bootstrap_alpha_tuning <- read_validation_chunks(
  "^bootstrap_alpha_tuning_chunk_[0-9]+\\.csv$",
  10
)

bootstrap_summary <- data.frame(
  bootstrap_resamples = nrow(bootstrap_performance),
  input_cpgs = bootstrap_performance$input_cpgs[1],
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
  "results/internal_validation/bootstrap_632_per_resample_summary.csv",
  row.names = FALSE
)

write.csv(
  bootstrap_oob_residuals,
  "results/internal_validation/bootstrap_oob_residuals.csv",
  row.names = FALSE
)

write.csv(
  bootstrap_selected_cpgs,
  "results/internal_validation/bootstrap_selected_cpgs.csv",
  row.names = FALSE
)

write.csv(
  bootstrap_alpha_tuning,
  "results/internal_validation/bootstrap_alpha_tuning.csv",
  row.names = FALSE
)

write.csv(
  bootstrap_summary,
  "results/internal_validation/bootstrap_632_summary.csv",
  row.names = FALSE
)

# Compare internal validation estimates with external GSE42861 performance
# The main question is whether internal validation estimates match independent external error
# Optimism gap means external error minus internal error

dir.create("results/analysis", recursive = TRUE, showWarnings = FALSE)

# Read a CSV if it exists and return NULL if the output has not been created yet
read_if_exists <- function(file) {
  if (!file.exists(file)) {
    warning("Missing file: ", file)
    return(NULL)
  }

  read.csv(file, stringsAsFactors = FALSE)
}

# Extract one metric from a summary file while allowing different column names across scripts
# For example some scripts use mae and repeated methods use mean_mae
summary_value <- function(summary_data, possible_names) {
  if (is.null(summary_data)) {
    return(NA_real_)
  }

  available_name <- possible_names[possible_names %in% names(summary_data)][1]

  if (is.na(available_name)) {
    return(NA_real_)
  }

  summary_data[[available_name]][1]
}

# Return a whole column if it exists, otherwise return NA values of the right length
# This keeps the final comparison working when older outputs lack newer metrics
optional_column <- function(summary_data, column_name) {
  if (column_name %in% names(summary_data)) {
    return(summary_data[[column_name]])
  }

  rep(NA_real_, nrow(summary_data))
}

# Return a character column if it exists, otherwise return NA character values
optional_character_column <- function(summary_data, column_name) {
  if (column_name %in% names(summary_data)) {
    return(summary_data[[column_name]])
  }

  rep(NA_character_, nrow(summary_data))
}

# Convert one internal validation summary into a standard format
# used_for_independent_validation is FALSE for apparent performance because it reuses training data
add_internal_method <- function(
  validation_method,
  estimate_type,
  summary_file,
  used_for_independent_validation
) {
  summary_data <- read_if_exists(summary_file)

  if (is.null(summary_data)) {
    return(NULL)
  }

  data.frame(
    validation_method = validation_method,
    estimate_type = estimate_type,
    used_for_independent_validation = used_for_independent_validation,
    internal_mae = summary_value(summary_data, c("mae", "mean_mae", "mean_oob_mae")),
    internal_rmse = summary_value(summary_data, c("rmse", "mean_rmse", "mean_oob_rmse")),
    internal_median_absolute_error = summary_value(
      summary_data,
      c("median_absolute_error", "mean_median_absolute_error", "mean_oob_median_absolute_error")
    ),
    internal_mean_error = summary_value(summary_data, c("mean_error", "mean_oob_error")),
    internal_correlation = summary_value(
      summary_data,
      c("correlation", "mean_correlation", "mean_oob_correlation")
    ),
    internal_r_squared = summary_value(
      summary_data,
      c("r_squared", "mean_r_squared", "mean_oob_r_squared")
    ),
    mean_selected_cpgs = summary_value(
      summary_data,
      c("selected_cpgs", "mean_selected_cpgs")
    ),
    mean_selected_alpha = summary_value(
      summary_data,
      c("selected_alpha", "mean_selected_alpha")
    ),
    mean_lambda_min = summary_value(summary_data, c("lambda_min", "mean_lambda_min")),
    mean_lambda_1se = summary_value(summary_data, c("lambda_1se", "mean_lambda_1se")),
    stringsAsFactors = FALSE
  )
}

# Build the internal-performance table across all validation methods
internal_methods <- do.call(
  rbind,
  list(
    add_internal_method(
      "full_data_tuned_alpha_baseline",
      "full_data_apparent_training_error",
      "results/internal_validation/apparent_performance_summary.csv",
      FALSE
    ),
    add_internal_method(
      "apparent_performance",
      "apparent_training_error",
      "results/internal_validation/apparent_performance_summary.csv",
      FALSE
    ),
    add_internal_method(
      "single_train_test_split",
      "single_holdout",
      "results/internal_validation/single_train_test_split_summary.csv",
      TRUE
    ),
    add_internal_method(
      "repeated_train_test_split",
      "repeated_holdout",
      "results/internal_validation/repeated_train_test_split_summary.csv",
      TRUE
    ),
    add_internal_method(
      "k_fold_cross_validation",
      "k_fold",
      "results/internal_validation/k_fold_summary.csv",
      TRUE
    ),
    add_internal_method(
      "repeated_k_fold_cross_validation",
      "repeated_k_fold",
      "results/internal_validation/repeated_k_fold_summary.csv",
      TRUE
    ),
    add_internal_method(
      "nested_cross_validation",
      "nested_outer_fold",
      "results/internal_validation/nested_cross_validation_summary.csv",
      TRUE
    ),
    add_internal_method(
      "bootstrap_oob",
      "bootstrap_out_of_bag",
      "results/internal_validation/bootstrap_632_summary.csv",
      TRUE
    )
  )
)

# Add richer internal calibration and age-acceleration summaries if they have been generated
# These come from 31_internal_validation_age_acceleration_analysis.r
internal_age_analysis <- read_if_exists(
  "results/analysis/internal_validation_aaa_raa_summary.csv"
)

if (!is.null(internal_age_analysis)) {
  # Map optional age-analysis columns to names used in the final comparison table
  age_analysis_name_map <- c(
    validation_method = "validation_method",
    calibration_intercept = "internal_calibration_intercept",
    calibration_slope = "internal_calibration_slope",
    percentile_95_absolute_error = "internal_percentile_95_absolute_error",
    max_absolute_error = "internal_max_absolute_error",
    aaa_empirical_95_lower = "internal_aaa_empirical_95_lower",
    aaa_empirical_95_upper = "internal_aaa_empirical_95_upper",
    raa_empirical_95_lower = "internal_raa_empirical_95_lower",
    raa_empirical_95_upper = "internal_raa_empirical_95_upper",
    aaa_age_slope = "internal_aaa_age_slope",
    absolute_error_age_slope = "internal_absolute_error_age_slope"
  )
  age_analysis_columns <- names(age_analysis_name_map)[
    names(age_analysis_name_map) %in% names(internal_age_analysis)
  ]
  internal_age_analysis <- internal_age_analysis[, age_analysis_columns, drop = FALSE]
  names(internal_age_analysis) <- age_analysis_name_map[age_analysis_columns]

  internal_methods <- merge(
    internal_methods,
    internal_age_analysis,
    by = "validation_method",
    all.x = TRUE
  )
}

# Read the standard external validation result for the final full-data clock
external_baseline <- read_if_exists(
  "results/external_validation/gse42861_external_validation_summary.csv"
)

if (!is.null(external_baseline)) {
  # Put the baseline external result into the same format as validation-informed clocks
  baseline_external <- data.frame(
    validation_method = "full_data_tuned_alpha_baseline",
    external_model_type = "final_full_data_clock",
    external_clock_type = "full_data_tuned_alpha",
    external_lambda_source = "cv.glmnet_lambda_min",
    external_samples = external_baseline$samples[1],
    external_mae = external_baseline$mae[1],
    external_rmse = external_baseline$rmse[1],
    external_median_absolute_error = external_baseline$median_absolute_error[1],
    external_mean_error = external_baseline$mean_error[1],
    external_correlation = external_baseline$correlation[1],
    external_r_squared = external_baseline$r_squared[1],
    external_calibration_intercept = summary_value(
      external_baseline,
      "calibration_intercept"
    ),
    external_calibration_slope = summary_value(
      external_baseline,
      "calibration_slope"
    ),
    external_percentile_95_absolute_error = summary_value(
      external_baseline,
      "percentile_95_absolute_error"
    ),
    external_max_absolute_error = summary_value(
      external_baseline,
      "max_absolute_error"
    ),
    stringsAsFactors = FALSE
  )
} else {
  baseline_external <- data.frame()
}

# Read external results for validation-informed clocks if those clocks have been fitted
validation_informed_external <- read_if_exists(
  "results/external_validation/validation_informed_clocks/gse42861_validation_informed_clock_performance.csv"
)

if (!is.null(validation_informed_external) && nrow(validation_informed_external) > 0) {
  # Standardise validation-informed external results for merging
  validation_informed_external <- data.frame(
    validation_method = validation_informed_external$validation_method,
    external_model_type = optional_character_column(
      validation_informed_external,
      "clock_type"
    ),
    external_clock_type = optional_character_column(
      validation_informed_external,
      "clock_type"
    ),
    external_lambda_source = optional_character_column(
      validation_informed_external,
      "lambda_source"
    ),
    external_samples = validation_informed_external$samples,
    external_mae = validation_informed_external$external_mae,
    external_rmse = validation_informed_external$external_rmse,
    external_median_absolute_error = validation_informed_external$median_absolute_error,
    external_mean_error = validation_informed_external$mean_AAA,
    external_correlation = validation_informed_external$correlation,
    external_r_squared = validation_informed_external$r_squared,
    external_calibration_intercept = optional_column(
      validation_informed_external,
      "calibration_intercept"
    ),
    external_calibration_slope = optional_column(
      validation_informed_external,
      "calibration_slope"
    ),
    external_percentile_95_absolute_error = optional_column(
      validation_informed_external,
      "percentile_95_absolute_error"
    ),
    external_max_absolute_error = optional_column(
      validation_informed_external,
      "max_absolute_error"
    ),
    stringsAsFactors = FALSE
  )
} else {
  validation_informed_external <- data.frame()
}

# Avoid duplicating the full-data baseline if it also appears in the validation-informed output
if (
  nrow(validation_informed_external) > 0 &&
    "full_data_tuned_alpha_baseline" %in% validation_informed_external$validation_method
) {
  external_results <- validation_informed_external
} else if (nrow(baseline_external) > 0 && nrow(validation_informed_external) > 0) {
  external_results <- rbind(baseline_external, validation_informed_external)
} else if (nrow(baseline_external) > 0) {
  external_results <- baseline_external
} else if (nrow(validation_informed_external) > 0) {
  external_results <- validation_informed_external
} else {
  external_results <- data.frame()
}

# Create an empty external table if no external results are available yet
if (nrow(external_results) == 0) {
  external_results <- data.frame(
    validation_method = character(),
    external_model_type = character(),
    external_clock_type = character(),
    external_lambda_source = character(),
    external_samples = numeric(),
    external_mae = numeric(),
    external_rmse = numeric(),
    external_median_absolute_error = numeric(),
    external_mean_error = numeric(),
    external_correlation = numeric(),
    external_r_squared = numeric(),
    external_calibration_intercept = numeric(),
    external_calibration_slope = numeric(),
    external_percentile_95_absolute_error = numeric(),
    external_max_absolute_error = numeric(),
    stringsAsFactors = FALSE
  )
}

# Merge internal and external metrics by validation method
comparison <- merge(
  internal_methods,
  external_results,
  by = "validation_method",
  all = TRUE
)

# Calculate internal-to-external disagreement
# A smaller absolute MAE gap means the internal estimate was closer to external performance
comparison$external_minus_internal_mae <- comparison$external_mae -
  comparison$internal_mae
comparison$absolute_mae_gap <- abs(comparison$external_minus_internal_mae)
comparison$external_minus_internal_rmse <- comparison$external_rmse -
  comparison$internal_rmse
comparison$absolute_rmse_gap <- abs(comparison$external_minus_internal_rmse)
comparison$calibration_slope_distance_from_1 <- abs(
  comparison$external_calibration_slope - 1
)
comparison$absolute_external_mean_error <- abs(comparison$external_mean_error)

# Rank methods only when both internal and external performance are available
ranking <- comparison[
  !is.na(comparison$external_mae) & !is.na(comparison$internal_mae),
]

if (nrow(ranking) > 0) {
  ranking <- ranking[order(
    ranking$absolute_mae_gap,
    ranking$external_mae,
    ranking$calibration_slope_distance_from_1,
    ranking$absolute_external_mean_error
  ), ]
  ranking$rank_by_external_agreement <- seq_len(nrow(ranking))
}

# Save the complete comparison and the ranked method table
write.csv(
  comparison,
  "results/analysis/final_internal_external_comparison.csv",
  row.names = FALSE
)

write.csv(
  ranking,
  "results/analysis/final_validation_method_ranking.csv",
  row.names = FALSE
)

print(comparison)

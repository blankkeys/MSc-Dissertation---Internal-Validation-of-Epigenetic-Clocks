# Compare internal validation estimates with external GSE42861 performance
# The main question is whether internal validation estimates match independent external error
# Optimism gap means external error minus internal error

summary_value <- function(summary_data, possible_names) {
  # Extract one metric while allowing different column names across scripts
  available_name <- possible_names[possible_names %in% names(summary_data)][1]

  if (is.na(available_name)) {
    return(NA_real_)
  }

  summary_data[[available_name]][1]
}

add_internal_method <- function(
  validation_method,
  estimate_type,
  summary_file,
  used_for_independent_validation
) {
  # Convert one internal validation summary into a standard format
  # Apparent performance is not independent validation because it reuses training data
  summary_data <- read.csv(summary_file, stringsAsFactors = FALSE)

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

add_internal_age_analysis <- function(internal_methods, internal_age_analysis) {
  # Add calibration and age-acceleration summaries from script 47 if available
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

  merge(
    internal_methods,
    internal_age_analysis,
    by = "validation_method",
    all.x = TRUE
  )
}

build_baseline_external <- function(external_baseline) {
  # Put the standard full-data external result into the same format as other clocks
  data.frame(
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
}

build_validation_informed_external <- function(external_data) {
  # Standardise validation-informed external results for merging
  data.frame(
    validation_method = external_data$validation_method,
    external_model_type = external_data$clock_type,
    external_clock_type = external_data$clock_type,
    external_lambda_source = external_data$lambda_source,
    external_samples = external_data$samples,
    external_mae = external_data$external_mae,
    external_rmse = external_data$external_rmse,
    external_median_absolute_error = external_data$median_absolute_error,
    external_mean_error = external_data$mean_AAA,
    external_correlation = external_data$correlation,
    external_r_squared = external_data$r_squared,
    external_calibration_intercept = external_data$calibration_intercept,
    external_calibration_slope = external_data$calibration_slope,
    external_percentile_95_absolute_error = external_data$percentile_95_absolute_error,
    external_max_absolute_error = external_data$max_absolute_error,
    stringsAsFactors = FALSE
  )
}

combine_external_results <- function(baseline_external, validation_informed_external) {
  # Avoid duplicating the full-data baseline if it also appears in validation-informed output
  if (
    nrow(validation_informed_external) > 0 &&
      "full_data_tuned_alpha_baseline" %in% validation_informed_external$validation_method
  ) {
    return(validation_informed_external)
  }

  rbind(baseline_external, validation_informed_external)
}

add_external_gap_metrics <- function(comparison) {
  # Smaller absolute gaps mean internal estimates were closer to external performance
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

  comparison
}

rank_validation_methods <- function(comparison) {
  # Rank only methods with both internal and external performance available
  ranking <- comparison[
    !is.na(comparison$external_mae) & !is.na(comparison$internal_mae),
  ]

  if (nrow(ranking) == 0) {
    return(ranking)
  }

  ranking <- ranking[order(
    ranking$absolute_mae_gap,
    ranking$external_mae,
    ranking$calibration_slope_distance_from_1,
    ranking$absolute_external_mean_error
  ), ]
  ranking$rank_by_external_agreement <- seq_len(nrow(ranking))

  ranking
}

dir.create("results/analysis", recursive = TRUE, showWarnings = FALSE)

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

internal_age_analysis <- read.csv(
  "results/analysis/internal_validation_aaa_raa_summary.csv",
  stringsAsFactors = FALSE
)
internal_methods <- add_internal_age_analysis(internal_methods, internal_age_analysis)

baseline_external <- build_baseline_external(
  read.csv(
    "results/external_validation/gse42861_external_validation_summary.csv",
    stringsAsFactors = FALSE
  )
)
validation_informed_external <- build_validation_informed_external(
  read.csv(
    "results/external_validation/validation_informed_clocks/gse42861_validation_informed_clock_performance.csv",
    stringsAsFactors = FALSE
  )
)

external_results <- combine_external_results(
  baseline_external,
  validation_informed_external
)

comparison <- merge(
  internal_methods,
  external_results,
  by = "validation_method",
  all = TRUE
)
comparison <- add_external_gap_metrics(comparison)
ranking <- rank_validation_methods(comparison)

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

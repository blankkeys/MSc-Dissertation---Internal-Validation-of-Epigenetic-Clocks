# Summarise uncertainty and subgroup performance for external GSE42861 validation
# This uses existing prediction files and does not refit any clocks

read_prediction_file <- function(file, analysis_type) {
  # Read an external prediction file and label whether it is strict or imputed
  predictions <- read.csv(file, stringsAsFactors = FALSE)

  predictions$analysis_type <- analysis_type
  predictions
}

add_external_sex <- function(predictions, external_metadata) {
  # Add sex labels by matching sample IDs and then GEO accessions if needed
  predictions$sex <- NA_character_

  if ("Sample_Name" %in% names(external_metadata)) {
    sample_matches <- match(predictions$sample_id, external_metadata$Sample_Name)
    matched <- !is.na(sample_matches)
    predictions$sex[matched] <- external_metadata$sex[sample_matches[matched]]
  }

  if ("geo_accession" %in% names(external_metadata) && "geo_accession" %in% names(predictions)) {
    geo_matches <- match(predictions$geo_accession, external_metadata$geo_accession)
    matched <- is.na(predictions$sex) & !is.na(geo_matches)
    predictions$sex[matched] <- external_metadata$sex[geo_matches[matched]]
  }

  predictions
}

add_age_acceleration_metrics <- function(predictions) {
  # AAA is predicted age minus chronological age
  # RAA is residual age acceleration after regressing predicted age on chronological age
  predictions$AAA <- predictions$predicted_age - predictions$age
  predictions$absolute_error <- abs(predictions$AAA)

  prediction_groups <- split(
    predictions,
    list(predictions$analysis_type, predictions$validation_method),
    drop = TRUE
  )

  prediction_groups <- lapply(prediction_groups, function(group_data) {
    group_data$RAA <- residuals(lm(predicted_age ~ age, data = group_data))
    group_data
  })

  predictions <- do.call(rbind, prediction_groups)
  rownames(predictions) <- NULL

  predictions
}

safe_cor <- function(x, y) {
  # Correlation is undefined if there are fewer than two values or no variation
  if (length(x) < 2 || sd(x) == 0 || sd(y) == 0) {
    return(NA_real_)
  }

  cor(x, y)
}

summarise_predictions <- function(data) {
  # Summarise external prediction accuracy, bias and calibration
  correlation <- safe_cor(data$predicted_age, data$age)
  calibration <- if (nrow(data) > 1 && sd(data$age) > 0) {
    coef(lm(predicted_age ~ age, data = data))
  } else {
    c("(Intercept)" = NA_real_, age = NA_real_)
  }

  data.frame(
    samples = nrow(data),
    mae = mean(data$absolute_error),
    median_absolute_error = median(data$absolute_error),
    rmse = sqrt(mean(data$AAA^2)),
    mean_error = mean(data$AAA),
    correlation = correlation,
    r_squared = correlation^2,
    calibration_intercept = calibration[1],
    calibration_slope = calibration[2],
    mean_AAA = mean(data$AAA),
    sd_AAA = sd(data$AAA),
    mean_RAA = mean(data$RAA),
    sd_RAA = sd(data$RAA),
    stringsAsFactors = FALSE
  )
}

bootstrap_ci <- function(data, n_bootstrap) {
  # Bootstrap confidence intervals resample external samples with replacement
  bootstrap_metrics <- replicate(n_bootstrap, {
    sampled_data <- data[sample(seq_len(nrow(data)), replace = TRUE), ]
    metrics <- summarise_predictions(sampled_data)

    c(
      mae = metrics$mae,
      rmse = metrics$rmse,
      correlation = metrics$correlation,
      calibration_intercept = metrics$calibration_intercept,
      calibration_slope = metrics$calibration_slope
    )
  })

  bootstrap_metrics <- t(bootstrap_metrics)

  data.frame(
    metric = colnames(bootstrap_metrics),
    lower_95_ci = apply(bootstrap_metrics, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper_95_ci = apply(bootstrap_metrics, 2, quantile, probs = 0.975, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

summarise_metric_ci <- function(prediction_groups, n_bootstrap) {
  # Calculate bootstrap confidence intervals for each external clock
  external_metric_ci <- do.call(rbind, lapply(prediction_groups, function(group_data) {
    point_estimate <- summarise_predictions(group_data)
    ci <- bootstrap_ci(group_data, n_bootstrap)

    data.frame(
      analysis_type = group_data$analysis_type[1],
      validation_method = group_data$validation_method[1],
      metric = ci$metric,
      point_estimate = as.numeric(unlist(point_estimate[1, ci$metric])),
      lower_95_ci = ci$lower_95_ci,
      upper_95_ci = ci$upper_95_ci,
      stringsAsFactors = FALSE
    )
  }))
  rownames(external_metric_ci) <- NULL

  external_metric_ci
}

compare_with_benchmark <- function(predictions, benchmark_method, n_bootstrap) {
  # Paired comparisons ask whether another clock has lower absolute error than the benchmark
  paired_results <- list()

  for (analysis_type in unique(predictions$analysis_type)) {
    analysis_data <- predictions[predictions$analysis_type == analysis_type, ]
    benchmark_data <- analysis_data[analysis_data$validation_method == benchmark_method, ]

    if (nrow(benchmark_data) == 0) {
      next
    }

    comparison_methods <- setdiff(unique(analysis_data$validation_method), benchmark_method)

    for (method in comparison_methods) {
      method_data <- analysis_data[analysis_data$validation_method == method, ]
      shared_samples <- intersect(benchmark_data$sample_id, method_data$sample_id)

      benchmark_error <- benchmark_data$absolute_error[
        match(shared_samples, benchmark_data$sample_id)
      ]
      method_error <- method_data$absolute_error[
        match(shared_samples, method_data$sample_id)
      ]
      error_difference <- method_error - benchmark_error

      bootstrap_difference <- replicate(n_bootstrap, {
        mean(error_difference[sample(seq_along(error_difference), replace = TRUE)])
      })

      paired_results[[length(paired_results) + 1]] <- data.frame(
        analysis_type = analysis_type,
        benchmark_method = benchmark_method,
        comparison_method = method,
        paired_samples = length(shared_samples),
        mean_absolute_error_difference = mean(error_difference),
        median_absolute_error_difference = median(error_difference),
        lower_95_ci = as.numeric(quantile(bootstrap_difference, probs = 0.025)),
        upper_95_ci = as.numeric(quantile(bootstrap_difference, probs = 0.975)),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(paired_results) == 0) {
    return(data.frame())
  }

  do.call(rbind, paired_results)
}

summarise_by_group <- function(predictions, group_column) {
  # Summarise performance within age bands or sex groups
  complete_group_data <- predictions[!is.na(predictions[[group_column]]), ]

  if (nrow(complete_group_data) == 0) {
    return(data.frame())
  }

  grouped_data <- split(
    complete_group_data,
    list(
      complete_group_data$analysis_type,
      complete_group_data$validation_method,
      complete_group_data[[group_column]]
    ),
    drop = TRUE
  )

  do.call(rbind, lapply(grouped_data, function(group_data) {
    metrics <- summarise_predictions(group_data)

    data.frame(
      analysis_type = group_data$analysis_type[1],
      validation_method = group_data$validation_method[1],
      group = group_data[[group_column]][1],
      metrics,
      stringsAsFactors = FALSE
    )
  }))
}

dir.create("results/analysis", recursive = TRUE, showWarnings = FALSE)

set.seed(20260707)
n_bootstrap <- 1000
benchmark_method <- "benchmark_alpha_0.5_cv_lambda"

strict_prediction_file <- paste0(
  "results/external_validation/validation_informed_clocks/",
  "gse42861_validation_informed_clock_predictions.csv"
)
imputed_prediction_file <- paste0(
  "results/external_validation/validation_informed_clocks/",
  "gse42861_mean_imputed_clock_predictions.csv"
)
external_metadata_file <- "data/GSE42861/gse42861_qc_ready_sample_sheet_controls.csv"

prediction_tables <- list(
  read_prediction_file(strict_prediction_file, "strict_primary"),
  read_prediction_file(imputed_prediction_file, "training_mean_imputed_sensitivity")
)

predictions <- do.call(rbind, prediction_tables)
external_metadata <- read.csv(external_metadata_file, stringsAsFactors = FALSE)
predictions <- add_external_sex(predictions, external_metadata)
predictions <- add_age_acceleration_metrics(predictions)

prediction_groups <- split(
  predictions,
  list(predictions$analysis_type, predictions$validation_method),
  drop = TRUE
)

external_metric_ci <- summarise_metric_ci(prediction_groups, n_bootstrap)
external_paired_comparisons <- compare_with_benchmark(
  predictions,
  benchmark_method,
  n_bootstrap
)

predictions$age_band <- cut(
  predictions$age,
  breaks = c(-Inf, 29, 44, 59, 74, Inf),
  labels = c("under_30", "30_44", "45_59", "60_74", "75_plus"),
  right = TRUE
)

external_age_band_performance <- summarise_by_group(predictions, "age_band")
external_sex_stratified_performance <- summarise_by_group(predictions, "sex")

write.csv(
  external_metric_ci,
  "results/analysis/external_metric_bootstrap_confidence_intervals.csv",
  row.names = FALSE
)

write.csv(
  external_paired_comparisons,
  "results/analysis/external_paired_clock_comparisons.csv",
  row.names = FALSE
)

write.csv(
  external_age_band_performance,
  "results/analysis/external_age_band_performance.csv",
  row.names = FALSE
)

write.csv(
  external_sex_stratified_performance,
  "results/analysis/external_sex_stratified_performance.csv",
  row.names = FALSE
)

print(external_metric_ci)
print(external_paired_comparisons)

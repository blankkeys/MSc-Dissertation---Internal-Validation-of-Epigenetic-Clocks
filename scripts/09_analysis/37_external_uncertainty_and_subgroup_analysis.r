# Summarise uncertainty and subgroup performance for external GSE42861 validation
# This uses existing prediction files and does not refit any clocks

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

read_prediction_file <- function(file, analysis_type) {
  if (!file.exists(file)) {
    return(NULL)
  }

  predictions <- read.csv(file, stringsAsFactors = FALSE)

  if (nrow(predictions) == 0) {
    return(NULL)
  }

  predictions$analysis_type <- analysis_type
  predictions
}

prediction_tables <- list(
  read_prediction_file(strict_prediction_file, "strict_primary"),
  read_prediction_file(imputed_prediction_file, "training_mean_imputed_sensitivity")
)
prediction_tables <- prediction_tables[!vapply(prediction_tables, is.null, logical(1))]

if (length(prediction_tables) == 0) {
  stop("No external prediction rows were available for uncertainty analysis")
}

predictions <- do.call(rbind, prediction_tables)
predictions$AAA <- predictions$predicted_age - predictions$age
predictions$absolute_error <- abs(predictions$AAA)

# Add sex labels from the external metadata for subgroup summaries
external_metadata <- read.csv(external_metadata_file, stringsAsFactors = FALSE)
predictions$sex <- external_metadata$sex[
  match(predictions$sample_id, external_metadata$Sample_Name)
]

# RAA is residual age acceleration after adjusting predicted age for chronological age
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

safe_cor <- function(x, y) {
  if (length(x) < 2 || sd(x) == 0 || sd(y) == 0) {
    return(NA_real_)
  }

  cor(x, y)
}

summarise_predictions <- function(data) {
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

bootstrap_ci <- function(data) {
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

# Bootstrap confidence intervals for each external clock
external_metric_ci <- do.call(rbind, lapply(prediction_groups, function(group_data) {
  point_estimate <- summarise_predictions(group_data)
  ci <- bootstrap_ci(group_data)

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

if (length(paired_results) > 0) {
  external_paired_comparisons <- do.call(rbind, paired_results)
} else {
  external_paired_comparisons <- data.frame()
}

summarise_by_group <- function(group_column) {
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

# Age-band and sex summaries check whether external error differs across subgroups
predictions$age_band <- cut(
  predictions$age,
  breaks = c(-Inf, 29, 44, 59, 74, Inf),
  labels = c("under_30", "30_44", "45_59", "60_74", "75_plus"),
  right = TRUE
)

external_age_band_performance <- summarise_by_group("age_band")
external_sex_stratified_performance <- summarise_by_group("sex")

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

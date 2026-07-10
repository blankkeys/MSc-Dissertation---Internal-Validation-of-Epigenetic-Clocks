# Compare uncertainty between internal validation methods
# This tests whether small MAE/RMSE differences between methods are meaningful
# Repeated methods are first collapsed to one average error per sample

read_internal_predictions <- function(validation_method, used_for_validation, file) {
  predictions <- read.csv(file, stringsAsFactors = FALSE)
  predictions$validation_method <- validation_method
  predictions$used_for_validation <- used_for_validation

  if (!"residual" %in% names(predictions)) {
    predictions$residual <- predictions$predicted_age - predictions$age
  }

  if (!"absolute_error" %in% names(predictions)) {
    predictions$absolute_error <- abs(predictions$residual)
  }

  predictions$squared_error <- predictions$residual^2
  predictions
}

summarise_sample_errors <- function(method_predictions) {
  # Collapse repeated predictions so each sample contributes once per method
  sample_errors <- aggregate(
    cbind(absolute_error, squared_error, residual) ~
      validation_method + used_for_validation + sample_id,
    data = method_predictions,
    FUN = mean
  )

  prediction_counts <- aggregate(
    method_predictions$sample_id,
    by = list(sample_id = method_predictions$sample_id),
    FUN = length
  )
  names(prediction_counts)[2] <- "prediction_events"

  merge(sample_errors, prediction_counts, by = "sample_id", all.x = TRUE)
}

summarise_method_error <- function(sample_errors) {
  data.frame(
    validation_method = sample_errors$validation_method[1],
    used_for_validation = sample_errors$used_for_validation[1],
    samples = nrow(sample_errors),
    mae = mean(sample_errors$absolute_error),
    rmse = sqrt(mean(sample_errors$squared_error)),
    mean_error = mean(sample_errors$residual),
    stringsAsFactors = FALSE
  )
}

bootstrap_method_ci <- function(sample_errors, n_bootstrap) {
  bootstrap_metrics <- replicate(n_bootstrap, {
    sampled_rows <- sample(seq_len(nrow(sample_errors)), replace = TRUE)
    sampled_data <- sample_errors[sampled_rows, ]

    c(
      mae = mean(sampled_data$absolute_error),
      rmse = sqrt(mean(sampled_data$squared_error)),
      mean_error = mean(sampled_data$residual)
    )
  })

  bootstrap_metrics <- t(bootstrap_metrics)
  point_estimate <- summarise_method_error(sample_errors)

  data.frame(
    validation_method = sample_errors$validation_method[1],
    used_for_validation = sample_errors$used_for_validation[1],
    metric = colnames(bootstrap_metrics),
    point_estimate = c(
      mae = point_estimate$mae,
      rmse = point_estimate$rmse,
      mean_error = point_estimate$mean_error
    ),
    lower_95_ci = apply(bootstrap_metrics, 2, quantile, probs = 0.025),
    upper_95_ci = apply(bootstrap_metrics, 2, quantile, probs = 0.975),
    stringsAsFactors = FALSE
  )
}

paired_method_comparison <- function(sample_level_errors, method_a, method_b, n_bootstrap) {
  method_a_errors <- sample_level_errors[
    sample_level_errors$validation_method == method_a,
  ]
  method_b_errors <- sample_level_errors[
    sample_level_errors$validation_method == method_b,
  ]

  shared_samples <- intersect(method_a_errors$sample_id, method_b_errors$sample_id)

  method_a_errors <- method_a_errors[
    match(shared_samples, method_a_errors$sample_id),
  ]
  method_b_errors <- method_b_errors[
    match(shared_samples, method_b_errors$sample_id),
  ]

  mae_difference <- method_a_errors$absolute_error - method_b_errors$absolute_error
  squared_error_difference <- method_a_errors$squared_error -
    method_b_errors$squared_error

  bootstrap_differences <- replicate(n_bootstrap, {
    sampled_rows <- sample(seq_along(shared_samples), replace = TRUE)

    c(
      mae_difference = mean(mae_difference[sampled_rows]),
      rmse_difference = sqrt(mean(method_a_errors$squared_error[sampled_rows])) -
        sqrt(mean(method_b_errors$squared_error[sampled_rows])),
      mean_squared_error_difference = mean(squared_error_difference[sampled_rows])
    )
  })

  bootstrap_differences <- t(bootstrap_differences)

  data.frame(
    method_a = method_a,
    method_b = method_b,
    shared_samples = length(shared_samples),
    mae_difference = mean(mae_difference),
    mae_difference_lower_95_ci = as.numeric(
      quantile(bootstrap_differences[, "mae_difference"], probs = 0.025)
    ),
    mae_difference_upper_95_ci = as.numeric(
      quantile(bootstrap_differences[, "mae_difference"], probs = 0.975)
    ),
    rmse_difference = sqrt(mean(method_a_errors$squared_error)) -
      sqrt(mean(method_b_errors$squared_error)),
    rmse_difference_lower_95_ci = as.numeric(
      quantile(bootstrap_differences[, "rmse_difference"], probs = 0.025)
    ),
    rmse_difference_upper_95_ci = as.numeric(
      quantile(bootstrap_differences[, "rmse_difference"], probs = 0.975)
    ),
    mean_squared_error_difference = mean(squared_error_difference),
    mean_squared_error_difference_lower_95_ci = as.numeric(
      quantile(bootstrap_differences[, "mean_squared_error_difference"], probs = 0.025)
    ),
    mean_squared_error_difference_upper_95_ci = as.numeric(
      quantile(bootstrap_differences[, "mean_squared_error_difference"], probs = 0.975)
    ),
    stringsAsFactors = FALSE
  )
}

dir.create("results/analysis", recursive = TRUE, showWarnings = FALSE)

set.seed(20260710)
n_bootstrap <- 1000

prediction_files <- data.frame(
  validation_method = c(
    "apparent_performance",
    "single_train_test_split",
    "repeated_train_test_split",
    "k_fold_cross_validation",
    "repeated_k_fold_cross_validation",
    "nested_cross_validation",
    "bootstrap_oob"
  ),
  used_for_validation = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  file = c(
    "results/internal_validation/apparent_performance_residuals.csv",
    "results/internal_validation/single_train_test_split_residuals.csv",
    "results/internal_validation/repeated_train_test_split_residuals.csv",
    "results/internal_validation/k_fold_residuals.csv",
    "results/internal_validation/repeated_k_fold_residuals.csv",
    "results/internal_validation/nested_cross_validation_residuals.csv",
    "results/internal_validation/bootstrap_oob_residuals.csv"
  ),
  stringsAsFactors = FALSE
)

internal_predictions <- do.call(
  rbind,
  lapply(seq_len(nrow(prediction_files)), function(i) {
    read_internal_predictions(
      prediction_files$validation_method[i],
      prediction_files$used_for_validation[i],
      prediction_files$file[i]
    )
  })
)

sample_level_errors <- do.call(
  rbind,
  lapply(
    split(internal_predictions, internal_predictions$validation_method),
    summarise_sample_errors
  )
)

method_error_summary <- do.call(
  rbind,
  lapply(
    split(sample_level_errors, sample_level_errors$validation_method),
    summarise_method_error
  )
)

method_metric_ci <- do.call(
  rbind,
  lapply(
    split(sample_level_errors, sample_level_errors$validation_method),
    bootstrap_method_ci,
    n_bootstrap = n_bootstrap
  )
)

validation_methods <- unique(
  sample_level_errors$validation_method[sample_level_errors$used_for_validation]
)

paired_comparisons <- data.frame()
for (i in seq_along(validation_methods)) {
  for (j in seq_along(validation_methods)) {
    if (i < j) {
      paired_comparisons <- rbind(
        paired_comparisons,
        paired_method_comparison(
          sample_level_errors,
          validation_methods[i],
          validation_methods[j],
          n_bootstrap
        )
      )
    }
  }
}

write.csv(
  sample_level_errors,
  "results/analysis/internal_validation_sample_level_errors.csv",
  row.names = FALSE
)

write.csv(
  method_error_summary,
  "results/analysis/internal_validation_method_error_summary.csv",
  row.names = FALSE
)

write.csv(
  method_metric_ci,
  "results/analysis/internal_validation_metric_bootstrap_confidence_intervals.csv",
  row.names = FALSE
)

write.csv(
  paired_comparisons,
  "results/analysis/internal_validation_paired_method_comparisons.csv",
  row.names = FALSE
)

print(method_metric_ci)
print(paired_comparisons)

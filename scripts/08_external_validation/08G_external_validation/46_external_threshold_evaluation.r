# Test internal-validation residual thresholds in the external GSE51032 dataset
# External residuals are classified as expected error, acceleration, or deceleration

threshold_file <- "results/internal_validation/internal_validation_residual_thresholds.csv"
external_predictions_file <- "results/external_validation/gse51032_external_validation_predictions.csv"

if (!file.exists(threshold_file)) {
  stop("Internal validation residual threshold file was not found")
}

if (!file.exists(external_predictions_file)) {
  stop("GSE51032 external validation prediction file was not found")
}

dir.create("results/external_validation", recursive = TRUE, showWarnings = FALSE)

thresholds <- read.csv(threshold_file)
external_predictions <- read.csv(external_predictions_file)

if (!"residual" %in% names(external_predictions)) {
  if ("age_error" %in% names(external_predictions)) {
    external_predictions$residual <- external_predictions$age_error
  } else {
    stop("External prediction file does not contain residual values")
  }
}

if (any(is.na(external_predictions$residual))) {
  stop("Missing residual values found in external validation predictions")
}

if (any(abs(external_predictions$residual -
  (external_predictions$predicted_age - external_predictions$age)) > 1e-8)) {
  stop("External residual values do not match predicted_age - age")
}

thresholds <- thresholds[thresholds$used_for_threshold == TRUE, ]

if (nrow(thresholds) == 0) {
  stop("No internal validation thresholds were available for external evaluation")
}

all_classifications <- data.frame()
coverage_summary <- data.frame()

for (i in seq_len(nrow(thresholds))) {
  method_thresholds <- thresholds[i, ]

  threshold_types <- data.frame(
    threshold_type = c("empirical_95", "parametric_95"),
    lower = c(
      method_thresholds$empirical_95_lower,
      method_thresholds$parametric_95_lower
    ),
    upper = c(
      method_thresholds$empirical_95_upper,
      method_thresholds$parametric_95_upper
    )
  )

  for (j in seq_len(nrow(threshold_types))) {
    lower_threshold <- threshold_types$lower[j]
    upper_threshold <- threshold_types$upper[j]

    if (is.na(lower_threshold) || is.na(upper_threshold)) {
      next
    }

    classification <- rep(
      "within_expected_error",
      nrow(external_predictions)
    )
    classification[external_predictions$residual > upper_threshold] <- "age_accelerated"
    classification[external_predictions$residual < lower_threshold] <- "age_decelerated"

    method_classifications <- data.frame(
      validation_method = method_thresholds$validation_method,
      threshold_type = threshold_types$threshold_type[j],
      lower_threshold = lower_threshold,
      upper_threshold = upper_threshold,
      sample_id = external_predictions$sample_id,
      geo_accession = external_predictions$geo_accession,
      age = external_predictions$age,
      predicted_age = external_predictions$predicted_age,
      residual = external_predictions$residual,
      absolute_error = external_predictions$absolute_error,
      threshold_classification = classification
    )

    external_samples <- nrow(method_classifications)
    within_count <- sum(classification == "within_expected_error")
    accelerated_count <- sum(classification == "age_accelerated")
    decelerated_count <- sum(classification == "age_decelerated")

    method_coverage <- data.frame(
      validation_method = method_thresholds$validation_method,
      threshold_type = threshold_types$threshold_type[j],
      external_samples = external_samples,
      within_expected_error_count = within_count,
      within_expected_error_proportion = within_count / external_samples,
      age_accelerated_count = accelerated_count,
      age_accelerated_proportion = accelerated_count / external_samples,
      age_decelerated_count = decelerated_count,
      age_decelerated_proportion = decelerated_count / external_samples,
      expected_coverage = 0.95,
      coverage_difference_from_95 = (within_count / external_samples) - 0.95
    )

    all_classifications <- rbind(
      all_classifications,
      method_classifications
    )
    coverage_summary <- rbind(coverage_summary, method_coverage)
  }
}

write.csv(
  all_classifications,
  "results/external_validation/gse51032_threshold_classifications.csv",
  row.names = FALSE
)

write.csv(
  coverage_summary,
  "results/external_validation/gse51032_internal_threshold_coverage_summary.csv",
  row.names = FALSE
)

print(coverage_summary)

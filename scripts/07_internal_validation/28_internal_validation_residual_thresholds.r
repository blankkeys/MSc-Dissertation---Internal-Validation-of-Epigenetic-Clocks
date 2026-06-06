# Summarise internal-validation residuals as age-acceleration thresholds
# Held-out residuals estimate expected prediction error in unseen samples

dir.create("results/internal_validation", recursive = TRUE, showWarnings = FALSE)

residual_files <- c(
  "results/internal_validation/apparent_performance_residuals.csv",
  "results/internal_validation/single_train_test_split_residuals.csv",
  "results/internal_validation/repeated_train_test_split_residuals.csv",
  "results/internal_validation/k_fold_residuals.csv",
  "results/internal_validation/repeated_k_fold_residuals.csv",
  "results/internal_validation/nested_cross_validation_residuals.csv",
  "results/internal_validation/bootstrap_oob_residuals.csv"
)

missing_files <- residual_files[!file.exists(residual_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing residual files: ",
    paste(missing_files, collapse = ", ")
  )
}

required_columns <- c(
  "validation_method",
  "resample_id",
  "sample_id",
  "geo_accession",
  "age",
  "predicted_age",
  "residual",
  "absolute_error"
)

all_residuals <- data.frame()

for (residual_file in residual_files) {
  residual_data <- read.csv(residual_file)

  missing_columns <- setdiff(required_columns, names(residual_data))

  if (length(missing_columns) > 0) {
    stop(
      "Missing required columns in ",
      residual_file,
      ": ",
      paste(missing_columns, collapse = ", ")
    )
  }

  residual_data <- residual_data[, required_columns]
  all_residuals <- rbind(all_residuals, residual_data)
}

if (any(is.na(all_residuals$residual))) {
  stop("Missing residual values found in internal validation residual files")
}

if (any(abs(all_residuals$residual - (all_residuals$predicted_age - all_residuals$age)) > 1e-8)) {
  stop("Residual values do not match predicted_age - age")
}

validation_methods <- unique(all_residuals$validation_method)
threshold_summary <- data.frame()

for (method_name in validation_methods) {
  method_residuals <- all_residuals[all_residuals$validation_method == method_name, ]
  residual <- method_residuals$residual
  used_for_threshold <- method_name != "apparent_performance"

  if (used_for_threshold) {
    empirical_limits <- quantile(residual, probs = c(0.025, 0.975), na.rm = TRUE)
    parametric_lower <- mean(residual) - (1.96 * sd(residual))
    parametric_upper <- mean(residual) + (1.96 * sd(residual))
  } else {
    empirical_limits <- c(NA, NA)
    parametric_lower <- NA
    parametric_upper <- NA
  }

  method_summary <- data.frame(
    validation_method = method_name,
    used_for_threshold = used_for_threshold,
    residual_count = nrow(method_residuals),
    unique_samples = length(unique(method_residuals$sample_id)),
    residual_mean = mean(residual),
    residual_sd = sd(residual),
    residual_median = median(residual),
    residual_iqr = IQR(residual),
    empirical_95_lower = as.numeric(empirical_limits[1]),
    empirical_95_upper = as.numeric(empirical_limits[2]),
    parametric_95_lower = parametric_lower,
    parametric_95_upper = parametric_upper,
    mae = mean(abs(residual)),
    median_absolute_error = median(abs(residual)),
    rmse = sqrt(mean(residual^2))
  )

  threshold_summary <- rbind(threshold_summary, method_summary)
}

write.csv(
  all_residuals,
  "results/internal_validation/internal_validation_residuals_combined.csv",
  row.names = FALSE
)

write.csv(
  threshold_summary,
  "results/internal_validation/internal_validation_residual_thresholds.csv",
  row.names = FALSE
)

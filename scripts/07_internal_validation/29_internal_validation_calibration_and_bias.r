# Calibration and age-bias checks for internal validation residuals
# This checks whether residual thresholds behave consistently across age

dir.create("results/internal_validation", recursive = TRUE, showWarnings = FALSE)

residual_file <- "results/internal_validation/internal_validation_residuals_combined.csv"

if (!file.exists(residual_file)) {
  stop("Combined internal validation residual file was not found")
}

residuals <- read.csv(residual_file)

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

missing_columns <- setdiff(required_columns, names(residuals))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (any(is.na(residuals$age)) || any(is.na(residuals$predicted_age))) {
  stop("Missing age or predicted age values were found")
}

if (any(abs(residuals$residual - (residuals$predicted_age - residuals$age)) > 1e-8)) {
  stop("Residual values do not match predicted_age - age")
}

# Create broad age bins to check whether prediction error changes across age
residuals$age_bin <- cut(
  residuals$age,
  breaks = c(-Inf, 29, 44, 59, 74, Inf),
  labels = c("under_30", "30_to_44", "45_to_59", "60_to_74", "75_plus"),
  right = TRUE
)

validation_methods <- unique(residuals$validation_method)
calibration_summary <- data.frame()
age_bin_summary <- data.frame()

for (method_name in validation_methods) {
  method_residuals <- residuals[residuals$validation_method == method_name, ]

  calibration_model <- lm(predicted_age ~ age, data = method_residuals)
  residual_age_model <- lm(residual ~ age, data = method_residuals)
  absolute_error_age_model <- lm(absolute_error ~ age, data = method_residuals)

  calibration_coefficients <- coef(calibration_model)
  residual_age_coefficients <- summary(residual_age_model)$coefficients
  absolute_error_age_coefficients <- summary(absolute_error_age_model)$coefficients

  method_calibration <- data.frame(
    validation_method = method_name,
    used_for_threshold = method_name != "apparent_performance",
    prediction_count = nrow(method_residuals),
    unique_samples = length(unique(method_residuals$sample_id)),
    calibration_intercept = calibration_coefficients[1],
    calibration_slope = calibration_coefficients[2],
    calibration_r_squared = summary(calibration_model)$r.squared,
    residual_age_slope = residual_age_coefficients["age", "Estimate"],
    residual_age_p_value = residual_age_coefficients["age", "Pr(>|t|)"],
    absolute_error_age_slope = absolute_error_age_coefficients["age", "Estimate"],
    absolute_error_age_p_value = absolute_error_age_coefficients["age", "Pr(>|t|)"],
    mean_residual = mean(method_residuals$residual),
    mae = mean(method_residuals$absolute_error),
    rmse = sqrt(mean(method_residuals$residual^2))
  )

  calibration_summary <- rbind(calibration_summary, method_calibration)

  for (bin_name in levels(residuals$age_bin)) {
    bin_residuals <- method_residuals[method_residuals$age_bin == bin_name, ]

    if (nrow(bin_residuals) == 0) {
      next
    }

    method_bin_summary <- data.frame(
      validation_method = method_name,
      used_for_threshold = method_name != "apparent_performance",
      age_bin = bin_name,
      prediction_count = nrow(bin_residuals),
      unique_samples = length(unique(bin_residuals$sample_id)),
      minimum_age = min(bin_residuals$age),
      maximum_age = max(bin_residuals$age),
      mean_residual = mean(bin_residuals$residual),
      median_residual = median(bin_residuals$residual),
      mae = mean(bin_residuals$absolute_error),
      median_absolute_error = median(bin_residuals$absolute_error),
      rmse = sqrt(mean(bin_residuals$residual^2))
    )

    age_bin_summary <- rbind(age_bin_summary, method_bin_summary)
  }
}

write.csv(
  calibration_summary,
  "results/internal_validation/internal_validation_calibration_summary.csv",
  row.names = FALSE
)

write.csv(
  age_bin_summary,
  "results/internal_validation/internal_validation_age_bin_error_summary.csv",
  row.names = FALSE
)

# Save simple diagnostic plots for each validation method
pdf("results/internal_validation/internal_validation_residual_vs_age.pdf")
for (method_name in validation_methods) {
  method_residuals <- residuals[residuals$validation_method == method_name, ]

  plot(
    method_residuals$age,
    method_residuals$residual,
    xlab = "Chronological age",
    ylab = "Residual",
    main = paste(method_name, "residual vs age"),
    pch = 16
  )
  abline(h = 0, col = "red")
  abline(lm(residual ~ age, data = method_residuals), col = "blue")
}
dev.off()

pdf("results/internal_validation/internal_validation_predicted_vs_age.pdf")
for (method_name in validation_methods) {
  method_residuals <- residuals[residuals$validation_method == method_name, ]

  plot(
    method_residuals$age,
    method_residuals$predicted_age,
    xlab = "Chronological age",
    ylab = "Predicted age",
    main = paste(method_name, "predicted age vs chronological age"),
    pch = 16
  )
  abline(0, 1, col = "red")
  abline(lm(predicted_age ~ age, data = method_residuals), col = "blue")
}
dev.off()

print(calibration_summary)

# Analyse age-acceleration and prediction-error behaviour across internal validation outputs
# AAA means absolute age acceleration, predicted DNAm age minus chronological age
# RAA means relative age acceleration, predicted DNAm age adjusted for chronological age using linear regression
# Calibration checks whether predicted age follows chronological age with slope close to 1 and intercept close to 0

empirical_lower <- function(x) {
  as.numeric(quantile(x, probs = 0.025, na.rm = TRUE))
}

empirical_upper <- function(x) {
  as.numeric(quantile(x, probs = 0.975, na.rm = TRUE))
}

read_residual_file <- function(validation_method, used_for_validation, file) {
  residuals <- read.csv(file, stringsAsFactors = FALSE)
  residuals$validation_method <- validation_method
  residuals$used_for_validation <- used_for_validation

  if (!"residual" %in% names(residuals)) {
    residuals$residual <- residuals$predicted_age - residuals$age
  }

  if (!"absolute_error" %in% names(residuals)) {
    residuals$absolute_error <- abs(residuals$residual)
  }

  residuals
}

add_age_acceleration <- function(predictions) {
  # AAA is the direct signed difference between predicted DNAm age and chronological age
  predictions$AAA <- predictions$residual

  # RAA is the residual after regressing predicted DNAm age on chronological age
  predictions$RAA <- NA_real_

  for (method in unique(predictions$validation_method)) {
    method_rows <- predictions$validation_method == method
    calibration_model <- lm(predicted_age ~ age, data = predictions[method_rows, ])
    predictions$RAA[method_rows] <- residuals(calibration_model)
  }

  predictions
}

add_age_bins <- function(predictions) {
  predictions$age_bin <- cut(
    predictions$age,
    breaks = c(-Inf, 29, 44, 59, 74, Inf),
    labels = c("14-29", "30-44", "45-59", "60-74", "75+"),
    right = TRUE
  )

  predictions
}

summarise_method <- function(method_predictions) {
  # calibration_model estimates predicted age = intercept + slope * chronological age
  calibration_model <- lm(predicted_age ~ age, data = method_predictions)

  # aaa_age_model checks whether signed error changes with chronological age
  aaa_age_model <- lm(AAA ~ age, data = method_predictions)

  # absolute_error_age_model checks whether error size changes with chronological age
  absolute_error_age_model <- lm(absolute_error ~ age, data = method_predictions)

  data.frame(
    validation_method = method_predictions$validation_method[1],
    used_for_validation = method_predictions$used_for_validation[1],
    predictions = nrow(method_predictions),
    unique_samples = length(unique(method_predictions$sample_id)),
    mae = mean(method_predictions$absolute_error),
    median_absolute_error = median(method_predictions$absolute_error),
    percentile_95_absolute_error = as.numeric(
      quantile(method_predictions$absolute_error, probs = 0.95, na.rm = TRUE)
    ),
    max_absolute_error = max(method_predictions$absolute_error),
    rmse = sqrt(mean(method_predictions$AAA^2)),
    correlation = cor(method_predictions$predicted_age, method_predictions$age),
    r_squared = cor(method_predictions$predicted_age, method_predictions$age)^2,
    calibration_intercept = coef(calibration_model)[1],
    calibration_slope = coef(calibration_model)[2],
    aaa_mean = mean(method_predictions$AAA),
    aaa_sd = sd(method_predictions$AAA),
    aaa_median = median(method_predictions$AAA),
    aaa_iqr = IQR(method_predictions$AAA),
    aaa_empirical_95_lower = empirical_lower(method_predictions$AAA),
    aaa_empirical_95_upper = empirical_upper(method_predictions$AAA),
    aaa_parametric_95_lower = mean(method_predictions$AAA) - 1.96 * sd(method_predictions$AAA),
    aaa_parametric_95_upper = mean(method_predictions$AAA) + 1.96 * sd(method_predictions$AAA),
    raa_mean = mean(method_predictions$RAA),
    raa_sd = sd(method_predictions$RAA),
    raa_median = median(method_predictions$RAA),
    raa_iqr = IQR(method_predictions$RAA),
    raa_empirical_95_lower = empirical_lower(method_predictions$RAA),
    raa_empirical_95_upper = empirical_upper(method_predictions$RAA),
    raa_parametric_95_lower = mean(method_predictions$RAA) - 1.96 * sd(method_predictions$RAA),
    raa_parametric_95_upper = mean(method_predictions$RAA) + 1.96 * sd(method_predictions$RAA),
    aaa_age_slope = coef(aaa_age_model)[2],
    absolute_error_age_slope = coef(absolute_error_age_model)[2],
    stringsAsFactors = FALSE
  )
}

summarise_age_bin <- function(age_bin_predictions) {
  data.frame(
    validation_method = age_bin_predictions$validation_method[1],
    used_for_validation = age_bin_predictions$used_for_validation[1],
    age_bin = age_bin_predictions$age_bin[1],
    predictions = nrow(age_bin_predictions),
    unique_samples = length(unique(age_bin_predictions$sample_id)),
    mean_age = mean(age_bin_predictions$age),
    mae = mean(age_bin_predictions$absolute_error),
    median_absolute_error = median(age_bin_predictions$absolute_error),
    percentile_95_absolute_error = as.numeric(
      quantile(age_bin_predictions$absolute_error, probs = 0.95, na.rm = TRUE)
    ),
    max_absolute_error = max(age_bin_predictions$absolute_error),
    rmse = sqrt(mean(age_bin_predictions$AAA^2)),
    mean_AAA = mean(age_bin_predictions$AAA),
    sd_AAA = sd(age_bin_predictions$AAA),
    mean_RAA = mean(age_bin_predictions$RAA),
    sd_RAA = sd(age_bin_predictions$RAA),
    stringsAsFactors = FALSE
  )
}

dir.create("results/analysis", recursive = TRUE, showWarnings = FALSE)

# List the prediction/residual files created by each internal validation method
# Apparent performance is included for comparison, but it is not independent validation
residual_files <- data.frame(
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

all_predictions <- do.call(
  rbind,
  lapply(seq_len(nrow(residual_files)), function(i) {
    read_residual_file(
      residual_files$validation_method[i],
      residual_files$used_for_validation[i],
      residual_files$file[i]
    )
  })
)

all_predictions <- add_age_acceleration(all_predictions)
all_predictions <- add_age_bins(all_predictions)

age_acceleration_summary <- do.call(
  rbind,
  lapply(
    split(all_predictions, all_predictions$validation_method),
    summarise_method
  )
)

age_bin_summary <- do.call(
  rbind,
  lapply(
    split(
      all_predictions,
      list(all_predictions$validation_method, all_predictions$age_bin),
      drop = TRUE
    ),
    summarise_age_bin
  )
)

write.csv(
  age_acceleration_summary,
  "results/analysis/internal_validation_aaa_raa_summary.csv",
  row.names = FALSE
)

write.csv(
  age_bin_summary,
  "results/analysis/internal_validation_age_bin_error_summary.csv",
  row.names = FALSE
)

write.csv(
  all_predictions,
  "results/analysis/internal_validation_age_acceleration_predictions.csv",
  row.names = FALSE
)

print(age_acceleration_summary)

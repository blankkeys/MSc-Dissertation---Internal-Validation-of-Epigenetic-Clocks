# External validation of the GSE87571 elastic-net clock in GSE42861
# This applies the trained model coefficients to GSE42861 without retraining it
# External validation tests whether the clock generalises to an independent dataset
library(glmnet)

# Load the externally processed GSE42861 beta-value matrix
external_beta_matrix <- readRDS("data/GSE42861/beta_matrix.rds")
# Load the GSE42861 control sample metadata
external_metadata <- read.csv("data/GSE42861/gse42861_qc_ready_sample_sheet_controls.csv")
# Load the final elastic-net clock trained in GSE87571
elastic_net_model <- readRDS("results/modelling/elastic_net_final_model.rds")
# Store how many CpGs are available in the external beta matrix before model matching
external_input_cpgs <- nrow(external_beta_matrix)
# Extract the CpGs selected by the final GSE87571 elastic-net model
# lambda.min is the penalty value chosen by cv.glmnet with the lowest cross-validation error
clock_coefficients <- as.matrix(coef(elastic_net_model, s = "lambda.min"))

# Put the model coefficients into a table
# lists selected CpGs and their coefficients
clock_coefficients <- data.frame(
  cpg = rownames(clock_coefficients),
  coefficient = as.numeric(clock_coefficients[, 1])
)

# Keep only CpGs with non-zero coefficients plus the intercept
clock_coefficients <- clock_coefficients[clock_coefficients$coefficient != 0, ]

# Store the model intercept separately because it is not a CpG
intercept <- clock_coefficients$coefficient[clock_coefficients$cpg == "(Intercept)"]

# Keep only the selected CpG coefficients
selected_coefficients <- clock_coefficients[clock_coefficients$cpg != "(Intercept)", ]

# Store the selected CpG IDs
selected_cpgs <- selected_coefficients$cpg

# Check whether any selected clock CpGs are missing from the external beta matrix
missing_cpgs <- setdiff(selected_cpgs, rownames(external_beta_matrix))

# Stop if the external data cannot support the trained clock
if (length(missing_cpgs) > 0) {
  stop("Some selected clock CpGs are missing from the GSE42861 beta matrix")
}

# Match metadata samples to the external beta matrix
external_metadata <- external_metadata[
  external_metadata$Sample_Name %in% colnames(external_beta_matrix) &
    !is.na(external_metadata$age),
]

# Stop if no external samples have both beta values and age metadata
if (nrow(external_metadata) == 0) {
  stop("No GSE42861 samples with matched beta values and age metadata were found")
}

# Use the selected CpGs in the same order as the clock coefficients
external_beta_matrix <- external_beta_matrix[
  selected_cpgs,
  match(external_metadata$Sample_Name, colnames(external_beta_matrix)),
  drop = FALSE
]

# Change the external beta matrix to have samples as rows and CpG sites as columns
# expected layout for glmnet 
x_external <- t(external_beta_matrix)

# Store chronological age as the external outcome
y_external <- external_metadata$age

# Predict chronological age in the independent GSE42861 samples
predicted_age <- intercept + x_external %*% selected_coefficients$coefficient

# Convert predictions from a matrix to a numeric vector
predicted_age <- as.numeric(predicted_age)

# Calculate prediction error as predicted age minus chronological age
# Positive residual means predicted DNAm age is older than chronological age
# Negative residual means predicted DNAm age is younger than chronological age
residual <- predicted_age - y_external

# Calculate absolute prediction error for each sample
absolute_error <- abs(residual)

# Calculate relative age acceleration as the residual after adjusting for age
# AAA is the direct residual, while RAA removes the linear relationship with chronological age
# RAA is useful because DNAm age can be systematically related to chronological age
calibration_model <- lm(predicted_age ~ y_external)
RAA <- residuals(calibration_model)

# Summarise external validation performance
# Calibration slope close to 1 and intercept close to 0 suggest well-calibrated predictions
# 95th percentile absolute error shows the size of larger but non-maximum errors
external_validation_summary <- data.frame(
  samples = length(y_external),
  external_input_cpgs = external_input_cpgs,
  selected_cpgs = length(selected_cpgs),
  mae = mean(absolute_error),
  median_absolute_error = median(absolute_error),
  percentile_95_absolute_error = as.numeric(
    quantile(absolute_error, probs = 0.95, na.rm = TRUE)
  ),
  max_absolute_error = max(absolute_error),
  rmse = sqrt(mean(residual^2)),
  mean_error = mean(residual),
  correlation = cor(predicted_age, y_external),
  r_squared = cor(predicted_age, y_external)^2,
  calibration_intercept = coef(calibration_model)[1],
  calibration_slope = coef(calibration_model)[2],
  mean_AAA = mean(residual),
  sd_AAA = sd(residual),
  mean_RAA = mean(RAA),
  sd_RAA = sd(RAA)
)

# Save one row per external sample for checking individual prediction errors
external_validation_predictions <- data.frame(
  sample_id = external_metadata$Sample_Name,
  geo_accession = external_metadata$geo_accession,
  age = y_external,
  predicted_age = predicted_age,
  residual = residual,
  AAA = residual,
  RAA = RAA,
  absolute_error = absolute_error
)

# Split samples into age groups to check whether performance changes across age ranges
external_validation_predictions$age_bin <- cut(
  external_validation_predictions$age,
  breaks = c(-Inf, 29, 44, 59, 74, Inf),
  labels = c("14-29", "30-44", "45-59", "60-74", "75+"),
  right = TRUE
)

# Summarise external prediction error separately within each age group
summarise_age_bin <- function(age_bin_predictions) {
  data.frame(
    age_bin = age_bin_predictions$age_bin[1],
    samples = nrow(age_bin_predictions),
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
    sd_RAA = sd(age_bin_predictions$RAA)
  )
}

external_age_bin_summary <- do.call(
  rbind,
  lapply(
    split(external_validation_predictions, external_validation_predictions$age_bin),
    summarise_age_bin
  )
)

dir.create("results/external_validation", recursive = TRUE, showWarnings = FALSE)

# Save the external validation summary metrics
write.csv(
  external_validation_summary,
  "results/external_validation/gse42861_external_validation_summary.csv",
  row.names = FALSE
)

# Save the sample-level external predictions
write.csv(
  external_validation_predictions,
  "results/external_validation/gse42861_external_validation_predictions.csv",
  row.names = FALSE
)

# Save age-bin external performance
write.csv(
  external_age_bin_summary,
  "results/external_validation/gse42861_external_validation_age_bin_summary.csv",
  row.names = FALSE
)

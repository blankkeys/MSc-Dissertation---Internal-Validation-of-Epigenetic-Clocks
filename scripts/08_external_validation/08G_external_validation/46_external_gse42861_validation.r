# External validation of the GSE87571 elastic-net clock in GSE42861
# This applies the trained model coefficients to GSE42861 without retraining it
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
residual <- predicted_age - y_external

# Calculate absolute prediction error for each sample
absolute_error <- abs(residual)

# Summarise external validation performance
external_validation_summary <- data.frame(
  samples = length(y_external),
  external_input_cpgs = external_input_cpgs,
  selected_cpgs = length(selected_cpgs),
  mae = mean(absolute_error),
  median_absolute_error = median(absolute_error),
  rmse = sqrt(mean(residual^2)),
  mean_error = mean(residual),
  correlation = cor(predicted_age, y_external),
  r_squared = cor(predicted_age, y_external)^2
)

external_validation_predictions <- data.frame(
  sample_id = external_metadata$Sample_Name,
  geo_accession = external_metadata$geo_accession,
  age = y_external,
  predicted_age = predicted_age,
  residual = residual,
  absolute_error = absolute_error
)


dir.create("results/external_validation", recursive = TRUE, showWarnings = FALSE)

# Save the external validation summary metrics
write.csv(
  external_validation_summary,
  "results/external_validation/gse42861_external_validation_summary.csv",
  row.names = FALSE
)


write.csv(
  external_validation_predictions,
  "results/external_validation/gse42861_external_validation_predictions.csv",
  row.names = FALSE
)

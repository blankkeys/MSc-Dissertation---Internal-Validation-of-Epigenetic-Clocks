# External validation of the GSE87571 elastic-net clock in GSE51032
# This applies the trained model coefficients to GSE51032 without retraining it

external_beta_matrix <- readRDS("data/GSE51032/beta_matrix.rds")
external_metadata <- read.csv("data/GSE51032/gse51032_qc_ready_sample_sheet_cancer_free.csv")
clock_coefficients <- read.csv("results/modelling/elastic_net_selected_cpg_coefficients.csv")

# Use the CpGs selected by the final GSE87571 elastic-net model
intercept <- clock_coefficients$coefficient[clock_coefficients$cpg == "(Intercept)"]
selected_coefficients <- clock_coefficients[clock_coefficients$cpg != "(Intercept)", ]
selected_cpgs <- selected_coefficients$cpg
missing_cpgs <- setdiff(selected_cpgs, rownames(external_beta_matrix))

if (length(missing_cpgs) > 0) {
  stop("Some selected clock CpGs are missing from the GSE51032 beta matrix")
}

# Match metadata samples to the external beta matrix
external_metadata$sample_id <- external_metadata$Sample_Name
external_metadata <- external_metadata[
  external_metadata$sample_id %in% colnames(external_beta_matrix) &
    !is.na(external_metadata$age),
]

if (nrow(external_metadata) == 0) {
  stop("No GSE51032 samples with matched beta values and age metadata were found")
}

# Use the selected CpGs in the same order as the clock coefficients
external_beta_matrix <- external_beta_matrix[
  selected_cpgs,
  match(external_metadata$sample_id, colnames(external_beta_matrix))
]

x_external <- t(external_beta_matrix)
y_external <- external_metadata$age

# Predict chronological age in the independent GSE51032 samples
predicted_age <- intercept + x_external %*% selected_coefficients$coefficient
predicted_age <- as.numeric(predicted_age)

external_validation_summary <- data.frame(
  samples = length(y_external),
  external_input_cpgs = nrow(external_beta_matrix),
  selected_cpgs = length(selected_cpgs),
  mae = mean(abs(predicted_age - y_external)),
  median_absolute_error = median(abs(predicted_age - y_external)),
  rmse = sqrt(mean((predicted_age - y_external)^2)),
  mean_error = mean(predicted_age - y_external),
  correlation = cor(predicted_age, y_external),
  r_squared = cor(predicted_age, y_external)^2
)

external_validation_predictions <- data.frame(
  sample_id = external_metadata$sample_id,
  geo_accession = external_metadata$geo_accession,
  age = y_external,
  predicted_age = predicted_age,
  age_error = predicted_age - y_external,
  absolute_error = abs(predicted_age - y_external)
)

dir.create("results/external_validation", recursive = TRUE, showWarnings = FALSE)

write.csv(
  external_validation_summary,
  "results/external_validation/gse51032_external_validation_summary.csv",
  row.names = FALSE
)

write.csv(
  external_validation_predictions,
  "results/external_validation/gse51032_external_validation_predictions.csv",
  row.names = FALSE
)

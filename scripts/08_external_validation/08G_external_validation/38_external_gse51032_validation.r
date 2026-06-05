# External validation of the GSE87571 elastic-net clock in GSE51032
# This applies the trained model to GSE51032 without retraining it

library(glmnet)

training_beta_matrix <- readRDS("data/GSE87571/beta_matrix_age_model.rds")
external_beta_matrix <- readRDS("data/GSE51032/beta_matrix.rds")
external_metadata <- read.csv("data/GSE51032/gse51032_qc_ready_sample_sheet_cancer_free.csv")
elastic_net_model <- readRDS("results/modelling/elastic_net_final_model.rds")

# Use the same CpG order used to train the GSE87571 model
training_cpgs <- rownames(training_beta_matrix)
external_cpgs <- rownames(external_beta_matrix)

missing_cpgs <- setdiff(training_cpgs, external_cpgs)

if (length(missing_cpgs) > 0) {
  stop("Some training CpGs are missing from the GSE51032 beta matrix")
}

# Match metadata samples to the external beta matrix
external_metadata$basename_id <- basename(external_metadata$Basename)

if (sum(external_metadata$Sample_Name %in% colnames(external_beta_matrix)) > 0) {
  external_metadata$sample_id <- external_metadata$Sample_Name
} else if (sum(external_metadata$basename_id %in% colnames(external_beta_matrix)) > 0) {
  external_metadata$sample_id <- external_metadata$basename_id
} else if (sum(external_metadata$geo_accession %in% colnames(external_beta_matrix)) > 0) {
  external_metadata$sample_id <- external_metadata$geo_accession
} else {
  stop("Could not match GSE51032 metadata samples to beta matrix columns")
}

external_metadata <- external_metadata[
  external_metadata$sample_id %in% colnames(external_beta_matrix) &
    !is.na(external_metadata$age),
]

if (nrow(external_metadata) == 0) {
  stop("No GSE51032 samples with matched beta values and age metadata were found")
}

# glmnet expects samples as rows and CpG sites as columns
x_external <- t(external_beta_matrix[
  training_cpgs,
  match(external_metadata$sample_id, colnames(external_beta_matrix))
])
y_external <- external_metadata$age

# Predict chronological age in the independent GSE51032 samples
predicted_age <- predict(
  elastic_net_model,
  newx = x_external,
  s = "lambda.min"
)

predicted_age <- as.numeric(predicted_age)

selected_cpgs <- sum(coef(elastic_net_model, s = "lambda.min")[-1, ] != 0)

external_validation_summary <- data.frame(
  samples = length(y_external),
  input_cpgs = ncol(x_external),
  selected_cpgs = selected_cpgs,
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

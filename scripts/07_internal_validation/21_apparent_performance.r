# Apparent performance for the elastic-net age prediction model
# This evaluates the model on the same samples used to train it.

library(glmnet)

dir.create("results/internal_validation", recursive = TRUE, showWarnings = FALSE)

beta_matrix <- readRDS("data/GSE87571/beta_matrix_age_model.rds")
metadata <- read.csv("data/GSE87571/modelling_metadata_age_model.csv")
elastic_net_model <- readRDS("results/modelling/elastic_net_final_model.rds")

# glmnet expects samples as rows and CpG sites as columns.
x <- t(beta_matrix[, match(metadata$sample_id, colnames(beta_matrix))])
y <- metadata$age

# Predict age using the trained elastic-net model (from 14) on the 
# training samples (apparent performance).
predicted_age <- predict(
  elastic_net_model,
  newx = x, 
  s = "lambda.min" # use the penalty value with the lowest cv error
)

predicted_age <- as.numeric(predicted_age)

apparent_residuals <- data.frame(
  validation_method = "apparent_performance",
  resample_id = "apparent",
  sample_id = metadata$sample_id,
  geo_accession = metadata$geo_accession,
  age = y,
  predicted_age = predicted_age,
  residual = predicted_age - y,
  absolute_error = abs(predicted_age - y)
)

apparent_performance <- data.frame(
  samples = length(y),
  input_cpgs = ncol(x),
  selected_cpgs = sum(coef(elastic_net_model, s = "lambda.min")[-1, ] != 0),
  mae = mean(abs(predicted_age - y)),
  median_absolute_error = median(abs(predicted_age - y)),
  rmse = sqrt(mean((predicted_age - y)^2)),
  mean_error = mean(predicted_age - y),
  correlation = cor(predicted_age, y),
  r_squared = cor(predicted_age, y)^2
)

# Save the selected CpGs from the final full-data model
selected_cpgs <- as.matrix(coef(elastic_net_model, s = "lambda.min"))
selected_cpgs <- data.frame(
  validation_method = "apparent_performance",
  resample_id = "apparent",
  cpg = rownames(selected_cpgs),
  coefficient = as.numeric(selected_cpgs[, 1])
)

selected_cpgs <- selected_cpgs[
  selected_cpgs$cpg != "(Intercept)" & selected_cpgs$coefficient != 0,
]

model_hyperparameters <- read.csv(
  "results/modelling/elastic_net_final_model_hyperparameters.csv"
)
selected_cpgs$selected_alpha <- model_hyperparameters$selected_alpha
selected_cpgs$lambda_min <- model_hyperparameters$lambda_min
selected_cpgs$lambda_1se <- model_hyperparameters$lambda_1se

write.csv(
  apparent_performance,
  "results/internal_validation/apparent_performance_summary.csv",
  row.names = FALSE
)

write.csv(
  apparent_residuals,
  "results/internal_validation/apparent_performance_residuals.csv",
  row.names = FALSE
)

write.csv(
  selected_cpgs,
  "results/internal_validation/apparent_performance_selected_cpgs.csv",
  row.names = FALSE
)

# Apparent performance for the elastic-net age prediction model
# This evaluates the model on the same samples used to train it.

# need to do feature selection, need to select the best ones, 
# need to do the same for the repeated k-fold cross-validation, need to do the same for the single train-test split validation, need to do the same for the leave-one-out cross-validation, need to do the same for the bootstrap validation, need to do the same for the nested cross-validation, need to do the same for the external validation, need to do the same for the independent validation, need to do the same for the prospective validation, need to do the same for the retrospective validation, need to do the same for the real-world validation, need to do the same for the clinical validation, need to do the same for the translational validation, need to do the same for the implementation validation, need to do the same for the impact validation, need to do the same for the cost-effectiveness validation, need to do the same for the ethical validation, need to do the same for the social validation, need to do the same for the regulatory validation, need to do the same for the policy validation, need to do the same for the public health validation.
# horvath has ~300 cpg after feature selection

# check list
# cross reactive
# sex (maybe remove or make sex based clocks)
# feature selection 
# consider python

library(glmnet)

beta_matrix <- readRDS("data/GSE87571/beta_matrix_age_model.rds")
metadata <- read.csv("data/GSE87571/modelling_metadata_age_model.csv")
elastic_net_model <- readRDS("results/modelling/elastic_net_final_model.rds")

# glmnet expects samples as rows and CpG sites as columns.
x <- t(beta_matrix[, match(metadata$sample_id, colnames(beta_matrix))])
y <- metadata$age

# Predict age using the trained elastic-net model on the 
# training samples (apparent performance).
predicted_age <- predict(
  elastic_net_model,
  newx = x,
  s = "lambda.min"
)

predicted_age <- as.numeric(predicted_age)

apparent_predictions <- data.frame(
  sample_id = metadata$sample_id,
  geo_accession = metadata$geo_accession,
  age = y,
  predicted_age = predicted_age,
  age_error = predicted_age - y,
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

write.csv(
  apparent_predictions,
  "results/internal_validation/apparent_performance_predictions.csv",
  row.names = FALSE
)

write.csv(
  apparent_performance,
  "results/internal_validation/apparent_performance_summary.csv",
  row.names = FALSE
)

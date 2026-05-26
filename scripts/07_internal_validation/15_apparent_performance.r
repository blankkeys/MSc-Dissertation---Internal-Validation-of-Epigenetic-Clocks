# Apparent performance for the elastic-net age prediction model
# This evaluates the model on the same samples used to train it.

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
  cpgs = ncol(x),
  mae = mean(abs(predicted_age - y)),
  rmse = sqrt(mean((predicted_age - y)^2)),
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

pdf("results/internal_validation/apparent_performance_predicted_vs_actual_age.pdf")
plot(
  y,
  predicted_age,
  xlab = "Chronological age",
  ylab = "Predicted age",
  main = "Apparent performance",
  pch = 16
)
abline(0, 1, col = "red")
dev.off()

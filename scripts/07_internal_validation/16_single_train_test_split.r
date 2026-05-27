# Single train-test split validation for the elastic-net age prediction model
# This trains the model on one subset and tests it on held-out samples.

library(glmnet)
library(rsample)

set.seed(123)

beta_matrix <- readRDS("data/GSE87571/beta_matrix_age_model.rds")
metadata <- read.csv("data/GSE87571/modelling_metadata_age_model.csv")

# glmnet expects samples as rows and CpG sites as columns.
x <- t(beta_matrix[, match(metadata$sample_id, colnames(beta_matrix))])

# Use 80% of samples for training and 20% for testing.
data_split <- initial_split(metadata, prop = 0.8)
train_metadata <- training(data_split)
test_metadata <- testing(data_split)

x_train <- x[train_metadata$sample_id, ]
y_train <- train_metadata$age

x_test <- x[test_metadata$sample_id, ]
y_test <- test_metadata$age

# Train the elastic-net model using the training samples only.
train_test_model <- cv.glmnet(
  x = x_train,
  y = y_train,
  alpha = 0.5,
  family = "gaussian"
)

# Predict age in the held-out test samples.
predicted_age <- predict(
  train_test_model,
  newx = x_test,
  s = "lambda.min"
)

predicted_age <- as.numeric(predicted_age)

train_test_predictions <- data.frame(
  sample_id = test_metadata$sample_id,
  geo_accession = test_metadata$geo_accession,
  age = y_test,
  predicted_age = predicted_age,
  age_error = predicted_age - y_test,
  absolute_error = abs(predicted_age - y_test)
)

train_test_performance <- data.frame(
  training_samples = length(y_train),
  test_samples = length(y_test),
  cpgs = ncol(x),
  mae = mean(abs(predicted_age - y_test)),
  rmse = sqrt(mean((predicted_age - y_test)^2)),
  correlation = cor(predicted_age, y_test),
  r_squared = cor(predicted_age, y_test)^2
)

write.csv(
  train_test_predictions,
  "results/internal_validation/single_train_test_split_predictions.csv",
  row.names = FALSE
)

write.csv(
  train_test_performance,
  "results/internal_validation/single_train_test_split_summary.csv",
  row.names = FALSE
)

pdf("results/internal_validation/single_train_test_split_predicted_vs_actual_age.pdf")
plot(
  y_test,
  predicted_age,
  xlab = "Chronological age",
  ylab = "Predicted age",
  main = "Single train-test split validation",
  pch = 16
)
abline(0, 1, col = "red")
dev.off()

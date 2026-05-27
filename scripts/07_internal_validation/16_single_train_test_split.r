# Single train-test split validation for the elastic-net age prediction model
# This trains the model on one subset and tests it on held-out samples

library(glmnet)
library(rsample)

set.seed(123)

beta_matrix <- readRDS("data/GSE87571/beta_matrix_age_model.rds")
metadata <- read.csv("data/GSE87571/modelling_metadata_age_model.csv")

# glmnet expects samples as rows and CpG sites as columns
x <- t(beta_matrix[, match(metadata$sample_id, colnames(beta_matrix))])

# Use 80% of samples for training and 20% for testing
# Stratify by age so training and test sets have similar age distributions
metadata_split <- initial_split(metadata, prop = 0.8, strata = age)
train_metadata <- training(metadata_split)
test_metadata <- testing(metadata_split)

# Match beta matrix columns to metadata rows for training and test sets
x_train <- x[train_metadata$sample_id, ] #x_train is the training data matrix with samples as rows and CpG sites as columns
y_train <- train_metadata$age # y_train is the vector of ages corresponding to the training samples

x_test <- x[test_metadata$sample_id, ] # x_test is the test data matrix with samples as rows and CpG sites as columns
y_test <- test_metadata$age # y_test is the vector of ages corresponding to the test samples

# Train the elastic-net model using the training samples only.
train_test_model <- cv.glmnet(
  x = x_train,
  y = y_train,
  alpha = 0.5,
  family = "gaussian"
)

# Predict age in the held-out test samples
predicted_age <- predict(
  train_test_model,
  newx = x_test,
  s = "lambda.min"
)

# Convert predicted age to numeric vector for performance calculations
predicted_age <- as.numeric(predicted_age)

# Create a data frame with predictions and errors for the test set
train_test_predictions <- data.frame(
  sample_id = test_metadata$sample_id,
  geo_accession = test_metadata$geo_accession,
  age = y_test,
  predicted_age = predicted_age,
  age_error = predicted_age - y_test,
  absolute_error = abs(predicted_age - y_test)
)

# Create a data frame summarising performance metrics for this single train-test split
train_test_performance <- data.frame(
  training_samples = length(y_train),
  test_samples = length(y_test),
  input_cpgs = ncol(x),
  selected_cpgs = sum(coef(train_test_model, s = "lambda.min")[-1, ] != 0),
  mae = mean(abs(predicted_age - y_test)),
  median_absolute_error = median(abs(predicted_age - y_test)),
  rmse = sqrt(mean((predicted_age - y_test)^2)),
  mean_error = mean(predicted_age - y_test),
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

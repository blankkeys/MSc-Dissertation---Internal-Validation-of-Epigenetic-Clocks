# Single train-test split validation for the elastic-net age prediction model
# This trains the model on one subset and tests it on held-out samples

library(glmnet)
library(rsample)
source("scripts/common/elastic_net_alpha_tuning.r")

dir.create("results/internal_validation", recursive = TRUE, showWarnings = FALSE)

# make the randomness reproducible
# (uses the same samples everytime)
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

# Tune alpha and lambda using the training samples only
alpha_grid <- c(0.25, 0.50, 0.75)
alpha_tuned_model <- tune_alpha_model(x_train, y_train, alpha_grid)
train_test_model <- alpha_tuned_model$model
alpha_tuning <- alpha_tuned_model$alpha_performance
alpha_tuning$validation_method <- "single_train_test_split"
alpha_tuning$resample_id <- "split_1"

# Save the CpGs selected by the training model
selected_cpgs <- get_selected_cpgs(
  train_test_model,
  "single_train_test_split",
  "split_1"
)
selected_cpgs$selected_alpha <- alpha_tuned_model$selected_alpha
selected_cpgs$lambda_min <- train_test_model$lambda.min
selected_cpgs$lambda_1se <- train_test_model$lambda.1se

# Predict age in the held-out test samples
predicted_age <- predict(
  train_test_model,
  newx = x_test,
  s = "lambda.min"
)

# Convert predicted age to numeric vector for performance calculations
predicted_age <- as.numeric(predicted_age)

# Save held-out residuals for age-acceleration threshold estimation
train_test_residuals <- data.frame(
  validation_method = "single_train_test_split",
  resample_id = "split_1",
  sample_id = test_metadata$sample_id,
  geo_accession = test_metadata$geo_accession,
  age = y_test,
  predicted_age = predicted_age,
  residual = predicted_age - y_test,
  absolute_error = abs(predicted_age - y_test)
)

# Create a data frame summarising performance metrics for this single train-test split
train_test_performance <- data.frame(
  training_samples = length(y_train),
  test_samples = length(y_test),
  input_cpgs = ncol(x),
  selected_alpha = alpha_tuned_model$selected_alpha,
  lambda_min = train_test_model$lambda.min,
  lambda_1se = train_test_model$lambda.1se,
  #how many cpg selected by elastic net model
  # coef(..)gets model coefficients at best lamda values chosen by cv
  # lambda.min menas use lamda value giving lowest cv error
  # -1 mens remove first row (the intercept, not a cpg just the model baseline)
  # !=0 checks which coe. are not zero
  # sum counts how many TRUE values there are
  selected_cpgs = sum(coef(train_test_model, s = "lambda.min")[-1, ] != 0),
  mae = mean(abs(predicted_age - y_test)),
  median_absolute_error = median(abs(predicted_age - y_test)),
  rmse = sqrt(mean((predicted_age - y_test)^2)),
  mean_error = mean(predicted_age - y_test),
  correlation = cor(predicted_age, y_test),
  r_squared = cor(predicted_age, y_test)^2
)

write.csv(
  train_test_performance,
  "results/internal_validation/single_train_test_split_summary.csv",
  row.names = FALSE
)

write.csv(
  train_test_residuals,
  "results/internal_validation/single_train_test_split_residuals.csv",
  row.names = FALSE
)

write.csv(
  selected_cpgs,
  "results/internal_validation/single_train_test_split_selected_cpgs.csv",
  row.names = FALSE
)

write.csv(
  alpha_tuning,
  "results/internal_validation/single_train_test_split_alpha_tuning.csv",
  row.names = FALSE
)

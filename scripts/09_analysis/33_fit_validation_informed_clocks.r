# Fit final clocks from validation-derived parameters
# Benchmark clocks use fixed alpha values with lambda selected by cv.glmnet
# Validation-informed clocks use alpha and lambda carried forward from internal validation
# glmnet is used for fixed-lambda clocks so lambda is not re-selected

library(glmnet)

set.seed(123)

dir.create(
  "results/modelling/validation_informed_clocks",
  recursive = TRUE,
  showWarnings = FALSE
)

selected_hyperparameters <- read.csv(
  "results/analysis/validation_informed_hyperparameter_selection.csv",
  stringsAsFactors = FALSE
)

beta_matrix <- readRDS("data/GSE87571/beta_matrix_age_model.rds")
metadata <- read.csv("data/GSE87571/modelling_metadata_age_model.csv")

x <- t(beta_matrix[, match(metadata$sample_id, colnames(beta_matrix))])
y <- metadata$age

# Use the same folds for the benchmark cv.glmnet fit each time this script runs
foldid <- sample(rep(seq_len(10), length.out = length(y)))

clock_models <- list()
clock_coefficients <- data.frame()
clock_summary <- data.frame()

# Save non-zero coefficients and metadata for one fitted clock
add_clock_outputs <- function(
  validation_method,
  clock_type,
  lambda_source,
  selected_alpha,
  lambda_min,
  lambda_1se,
  coefficients,
  internal_estimated_mae,
  internal_estimated_rmse
) {
  coefficients <- data.frame(
    validation_method = validation_method,
    clock_type = clock_type,
    lambda_source = lambda_source,
    cpg = rownames(coefficients),
    coefficient = as.numeric(coefficients[, 1]),
    selected_alpha = selected_alpha,
    lambda_min = lambda_min,
    lambda_1se = lambda_1se,
    stringsAsFactors = FALSE
  )
  coefficients <- coefficients[coefficients$coefficient != 0, ]

  list(
    coefficients = coefficients,
    summary = data.frame(
      validation_method = validation_method,
      clock_type = clock_type,
      lambda_source = lambda_source,
      selected_alpha = selected_alpha,
      lambda_min = lambda_min,
      lambda_1se = lambda_1se,
      selected_cpgs = sum(coefficients$cpg != "(Intercept)"),
      internal_estimated_mae = internal_estimated_mae,
      internal_estimated_rmse = internal_estimated_rmse,
      stringsAsFactors = FALSE
    )
  )
}

# Conventional benchmark clock using alpha 0.5
# This tests the standard balanced elastic-net approach against validation-derived clocks
benchmark_alpha_0.5_model <- cv.glmnet(
  x = x,
  y = y,
  alpha = 0.5,
  family = "gaussian",
  foldid = foldid
)

benchmark_alpha_0.5_outputs <- add_clock_outputs(
  validation_method = "benchmark_alpha_0.5_cv_lambda",
  clock_type = "benchmark_alpha_0.5",
  lambda_source = "cv.glmnet_lambda_min",
  selected_alpha = 0.5,
  lambda_min = benchmark_alpha_0.5_model$lambda.min,
  lambda_1se = benchmark_alpha_0.5_model$lambda.1se,
  coefficients = as.matrix(coef(benchmark_alpha_0.5_model, s = "lambda.min")),
  internal_estimated_mae = NA_real_,
  internal_estimated_rmse = NA_real_
)

clock_models[["benchmark_alpha_0.5_cv_lambda"]] <- benchmark_alpha_0.5_model
clock_coefficients <- rbind(clock_coefficients, benchmark_alpha_0.5_outputs$coefficients)
clock_summary <- rbind(clock_summary, benchmark_alpha_0.5_outputs$summary)

# Conventional benchmark clock using alpha 0.25
# This tests whether the tuned alpha value helps even when lambda is selected by cv.glmnet
benchmark_alpha_0.25_model <- cv.glmnet(
  x = x,
  y = y,
  alpha = 0.25,
  family = "gaussian",
  foldid = foldid
)

benchmark_alpha_0.25_outputs <- add_clock_outputs(
  validation_method = "benchmark_alpha_0.25_cv_lambda",
  clock_type = "benchmark_alpha_0.25",
  lambda_source = "cv.glmnet_lambda_min",
  selected_alpha = 0.25,
  lambda_min = benchmark_alpha_0.25_model$lambda.min,
  lambda_1se = benchmark_alpha_0.25_model$lambda.1se,
  coefficients = as.matrix(coef(benchmark_alpha_0.25_model, s = "lambda.min")),
  internal_estimated_mae = NA_real_,
  internal_estimated_rmse = NA_real_
)

clock_models[["benchmark_alpha_0.25_cv_lambda"]] <- benchmark_alpha_0.25_model
clock_coefficients <- rbind(clock_coefficients, benchmark_alpha_0.25_outputs$coefficients)
clock_summary <- rbind(clock_summary, benchmark_alpha_0.25_outputs$summary)

# Validation-informed clocks
# Each clock uses the alpha/lambda pair selected from one internal validation method
for (i in seq_len(nrow(selected_hyperparameters))) {
  method <- selected_hyperparameters$validation_method[i]
  selected_alpha <- selected_hyperparameters$selected_alpha[i]
  selected_lambda <- selected_hyperparameters$selected_lambda_min[i]

  clock_model <- glmnet(
    x = x,
    y = y,
    alpha = selected_alpha,
    lambda = selected_lambda,
    family = "gaussian"
  )

  clock_outputs <- add_clock_outputs(
    validation_method = method,
    clock_type = "validation_parameter_clock",
    lambda_source = "internal_validation_median_lambda_min",
    selected_alpha = selected_alpha,
    lambda_min = selected_lambda,
    lambda_1se = selected_hyperparameters$selected_lambda_1se[i],
    coefficients = as.matrix(coef(clock_model, s = selected_lambda)),
    internal_estimated_mae = selected_hyperparameters$internal_estimated_mae[i],
    internal_estimated_rmse = selected_hyperparameters$internal_estimated_rmse[i]
  )

  clock_models[[method]] <- clock_model
  clock_coefficients <- rbind(clock_coefficients, clock_outputs$coefficients)
  clock_summary <- rbind(clock_summary, clock_outputs$summary)
}

saveRDS(
  clock_models,
  "results/modelling/validation_informed_clocks/validation_informed_clock_models.rds"
)

write.csv(
  clock_coefficients,
  "results/modelling/validation_informed_clocks/validation_informed_clock_coefficients.csv",
  row.names = FALSE
)

write.csv(
  clock_summary,
  "results/modelling/validation_informed_clocks/validation_informed_clock_summary.csv",
  row.names = FALSE
)

print(clock_summary)

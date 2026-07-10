# Fit final clocks from validation-derived parameters
# Benchmark clocks use fixed alpha values with lambda selected by cv.glmnet
# Validation-informed clocks use alpha and lambda carried forward from internal validation
# glmnet is used for fixed-lambda clocks so lambda is not re-selected

library(glmnet)

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

fit_cv_benchmark_clock <- function(
  x,
  y,
  foldid,
  selected_alpha,
  validation_method,
  clock_type
) {
  clock_model <- cv.glmnet(
    x = x,
    y = y,
    alpha = selected_alpha,
    family = "gaussian",
    foldid = foldid
  )

  clock_outputs <- add_clock_outputs(
    validation_method = validation_method,
    clock_type = clock_type,
    lambda_source = "cv.glmnet_lambda_min",
    selected_alpha = selected_alpha,
    lambda_min = clock_model$lambda.min,
    lambda_1se = clock_model$lambda.1se,
    coefficients = as.matrix(coef(clock_model, s = "lambda.min")),
    internal_estimated_mae = NA_real_,
    internal_estimated_rmse = NA_real_
  )

  list(model = clock_model, outputs = clock_outputs)
}

fit_validation_parameter_clock <- function(x, y, selected_hyperparameters, row_index) {
  method <- selected_hyperparameters$validation_method[row_index]
  selected_alpha <- selected_hyperparameters$selected_alpha[row_index]
  selected_lambda <- selected_hyperparameters$selected_lambda_min[row_index]

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
    lambda_1se = selected_hyperparameters$selected_lambda_1se[row_index],
    coefficients = as.matrix(coef(clock_model, s = selected_lambda)),
    internal_estimated_mae = selected_hyperparameters$internal_estimated_mae[row_index],
    internal_estimated_rmse = selected_hyperparameters$internal_estimated_rmse[row_index]
  )

  list(model = clock_model, outputs = clock_outputs, method = method)
}

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
foldid <- sample(rep(seq_len(10), length.out = length(y)))

clock_models <- list()
clock_coefficients <- data.frame()
clock_summary <- data.frame()

# Conventional benchmark clock using alpha 0.5
benchmark_alpha_0.5 <- fit_cv_benchmark_clock(
  x = x,
  y = y,
  foldid = foldid,
  selected_alpha = 0.5,
  validation_method = "benchmark_alpha_0.5_cv_lambda",
  clock_type = "benchmark_alpha_0.5"
)
clock_models[["benchmark_alpha_0.5_cv_lambda"]] <- benchmark_alpha_0.5$model
clock_coefficients <- rbind(clock_coefficients, benchmark_alpha_0.5$outputs$coefficients)
clock_summary <- rbind(clock_summary, benchmark_alpha_0.5$outputs$summary)

# Conventional benchmark clock using alpha 0.25
benchmark_alpha_0.25 <- fit_cv_benchmark_clock(
  x = x,
  y = y,
  foldid = foldid,
  selected_alpha = 0.25,
  validation_method = "benchmark_alpha_0.25_cv_lambda",
  clock_type = "benchmark_alpha_0.25"
)
clock_models[["benchmark_alpha_0.25_cv_lambda"]] <- benchmark_alpha_0.25$model
clock_coefficients <- rbind(clock_coefficients, benchmark_alpha_0.25$outputs$coefficients)
clock_summary <- rbind(clock_summary, benchmark_alpha_0.25$outputs$summary)

# Validation-informed clocks
for (i in seq_len(nrow(selected_hyperparameters))) {
  validation_clock <- fit_validation_parameter_clock(x, y, selected_hyperparameters, i)
  clock_models[[validation_clock$method]] <- validation_clock$model
  clock_coefficients <- rbind(clock_coefficients, validation_clock$outputs$coefficients)
  clock_summary <- rbind(clock_summary, validation_clock$outputs$summary)
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

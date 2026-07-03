# Refit one final full-data clock for each validation-informed alpha

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

# Use the same folds for all full-data refits so models differ only by alpha
foldid <- sample(rep(seq_len(10), length.out = length(y)))

clock_models <- list()
clock_coefficients <- data.frame()
clock_summary <- data.frame()

for (i in seq_len(nrow(selected_hyperparameters))) {
  method <- selected_hyperparameters$validation_method[i]
  selected_alpha <- selected_hyperparameters$selected_alpha[i]

  clock_model <- cv.glmnet(
    x = x,
    y = y,
    alpha = selected_alpha,
    family = "gaussian",
    foldid = foldid
  )

  coefficients <- as.matrix(coef(clock_model, s = "lambda.min"))
  coefficients <- data.frame(
    validation_method = method,
    cpg = rownames(coefficients),
    coefficient = as.numeric(coefficients[, 1]),
    selected_alpha = selected_alpha,
    lambda_min = clock_model$lambda.min,
    lambda_1se = clock_model$lambda.1se,
    stringsAsFactors = FALSE
  )
  coefficients <- coefficients[coefficients$coefficient != 0, ]

  clock_models[[method]] <- clock_model
  clock_coefficients <- rbind(clock_coefficients, coefficients)
  clock_summary <- rbind(
    clock_summary,
    data.frame(
      validation_method = method,
      selected_alpha = selected_alpha,
      lambda_min = clock_model$lambda.min,
      lambda_1se = clock_model$lambda.1se,
      selected_cpgs = sum(coefficients$cpg != "(Intercept)"),
      internal_estimated_mae = selected_hyperparameters$internal_estimated_mae[i],
      internal_estimated_rmse = selected_hyperparameters$internal_estimated_rmse[i],
      stringsAsFactors = FALSE
    )
  )
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

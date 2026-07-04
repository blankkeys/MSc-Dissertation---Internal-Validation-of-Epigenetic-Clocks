# Elastic-net age prediction model for GSE87571
# Train an epigenetic clock using age-available samples.

library(glmnet)
source("scripts/common/elastic_net_alpha_tuning.r")

# set seed makes the random processes in glmnet reproducible,
# so we get the same model each time we run this script
set.seed(123)

beta_matrix <- readRDS("data/GSE87571/beta_matrix_age_model.rds")
metadata <- read.csv("data/GSE87571/modelling_metadata_age_model.csv")

# Match beta matrix columns to metadata rows.
# glmnet expects samples as rows and CpG sites as columns.
x <- t(beta_matrix[, match(metadata$sample_id, colnames(beta_matrix))])
y <- metadata$age

# Test ridge-like, balanced and lasso-like elastic-net alpha values.
# Lambda is selected by cv.glmnet within each alpha value.
alpha_grid <- c(0.25, 0.50, 0.75)
alpha_tuned_model <- tune_alpha_model(x, y, alpha_grid)
elastic_net_model <- alpha_tuned_model$model

dir.create("results/modelling", recursive = TRUE, showWarnings = FALSE)

saveRDS(
  elastic_net_model,
  "results/modelling/elastic_net_final_model.rds"
)

write.csv(
  alpha_tuned_model$alpha_performance,
  "results/modelling/elastic_net_alpha_tuning_summary.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    selected_alpha = alpha_tuned_model$selected_alpha,
    lambda_min = elastic_net_model$lambda.min,
    lambda_1se = elastic_net_model$lambda.1se
  ),
  "results/modelling/elastic_net_final_model_hyperparameters.csv",
  row.names = FALSE
)

# tells me which CpG sites were selected by the model and their coefficients
coefficients <- as.matrix(coef(elastic_net_model, s = "lambda.min"))
coefficients <- data.frame(
  cpg = rownames(coefficients),
  coefficient = as.numeric(coefficients[, 1]) #take coef values from first column, turn into numeric values
)

# only rows wher coef is not 0
coefficients <- coefficients[coefficients$coefficient != 0, ]

write.csv(
  coefficients,
  "results/modelling/elastic_net_selected_cpg_coefficients.csv",
  row.names = FALSE
)

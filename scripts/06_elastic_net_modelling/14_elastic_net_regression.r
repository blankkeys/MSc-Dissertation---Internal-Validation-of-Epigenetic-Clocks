# Elastic-net age prediction model for GSE87571
# Train an epigenetic clock using age-available samples.

library(glmnet)

# set seed makes the random processes in glmnet reproducible,
# so we get the same model each time we run this script
set.seed(123)

beta_matrix <- readRDS("data/GSE87571/beta_matrix_age_model.rds")
metadata <- read.csv("data/GSE87571/modelling_metadata_age_model.csv")

# Match beta matrix columns to metadata rows.
# glmnet expects samples as rows and CpG sites as columns.
x <- t(beta_matrix[, match(metadata$sample_id, colnames(beta_matrix))])
y <- metadata$age

# Train the elastic-net model.
# cv.glmnet uses cross-validation to choose the lambda penalty value.
elastic_net_model <- cv.glmnet(
  x = x,
  y = y,
  alpha = 0.5, # hovarth clock used alpha 0.5
  family = "gaussian" # gaussian family for regression (age prediction)
)

dir.create("results/modelling", recursive = TRUE, showWarnings = FALSE)

saveRDS(
  elastic_net_model,
  "results/modelling/elastic_net_final_model.rds"
)

# tells me which CpG sites were selected by the model and their coefficients
coefficients <- as.matrix(coef(elastic_net_model, s = "lambda.min"))
coefficients <- data.frame(
  cpg = rownames(coefficients),
  coefficient = as.numeric(coefficients[, 1])
)

coefficients <- coefficients[coefficients$coefficient != 0, ]

write.csv(
  coefficients,
  "results/modelling/elastic_net_selected_cpg_coefficients.csv",
  row.names = FALSE
)

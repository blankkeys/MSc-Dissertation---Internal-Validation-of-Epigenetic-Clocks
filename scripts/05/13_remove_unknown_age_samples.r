# Remove unknown-age samples for GSE87571 age modelling
# Keep only samples with chronological age available.

beta_matrix <- readRDS("data/GSE87571/beta_matrix.rds")
modelling_metadata <- read.csv("data/GSE87571/modelling_metadata.csv")

# Keep samples that can be used for supervised age prediction.
modelling_metadata_age <- modelling_metadata[
  modelling_metadata$included_in_age_modelling == TRUE,
]

# Reorder beta matrix columns to match the metadata rows.
sample_match <- match(modelling_metadata_age$sample_id, colnames(beta_matrix))
beta_matrix_age <- beta_matrix[, sample_match]

saveRDS(
  beta_matrix_age,
  "data/GSE87571/beta_matrix_age_model.rds"
)

write.csv(
  modelling_metadata_age,
  "data/GSE87571/modelling_metadata_age_model.csv",
  row.names = FALSE
)

dim(beta_matrix_age)

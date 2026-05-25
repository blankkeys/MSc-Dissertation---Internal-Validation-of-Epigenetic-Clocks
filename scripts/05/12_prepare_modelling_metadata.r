# Prepare modelling metadata for GSE87571
# Match age and sex metadata to the final beta matrix samples.

library(GEOquery)

# Load the final beta matrix to get the sample IDs used for modelling.
beta_matrix <- readRDS("data/GSE87571/beta_matrix.rds")

# Extract the GEO sample accession from each beta matrix column name.
modelling_metadata <- data.frame(
  sample_id = colnames(beta_matrix),
  geo_accession = sub("_.*$", "", colnames(beta_matrix))
)

# Read age and sex metadata from the local GEO Series Matrix.
gse <- getGEO(
  filename = "data/GSE87571/GSE87571_series_matrix.txt.gz",
  getGPL = FALSE
)
geo_metadata <- pData(gse)

series_metadata <- data.frame(
  geo_accession = rownames(geo_metadata),
  age = as.numeric(sub("age: ", "", geo_metadata$characteristics_ch1.1)),
  sex = tolower(sub("gender: ", "", geo_metadata$characteristics_ch1))
)

# Match age and sex to the beta matrix samples using the GEO accession.
metadata_match <- match(modelling_metadata$geo_accession, series_metadata$geo_accession)
modelling_metadata$age <- series_metadata$age[metadata_match]
modelling_metadata$sex <- series_metadata$sex[metadata_match]

# Samples with age available can be used for age prediction modelling.
modelling_metadata$included_in_age_modelling <- !is.na(modelling_metadata$age)

write.csv(
  modelling_metadata,
  "data/GSE87571/modelling_metadata.csv",
  row.names = FALSE
)

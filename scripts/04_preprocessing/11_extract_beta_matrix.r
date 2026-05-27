# This script extracts the beta matrix from the filtered MethylationSet 
# object and saves it for downstream analysis
# It also performs some basic checks on the beta matrix and saves a summary of its properties

library (minfi)

# Load the filtered MethylationSet object after detection, annotation,
# cross-reactive, and autosomal probe filtering.
filtered_mset_file <- readRDS("data/GSE87571/mset_normalised_filtered_annotation_crossreactive_autosomal.rds")

# Check that the loaded object is a MethylationSet and has the expected dimensions
print(dim(filtered_mset_file))

# Extract the beta matrix from the filtered MethylationSet object
beta_matrix <- getBeta(filtered_mset_file)

# Perform checks on the beta matrix
missing_value_count <- sum(is.na(beta_matrix))
duplicate_probe_count <- sum(duplicated(rownames(beta_matrix)))
duplicate_sample_count <- sum(duplicated(colnames(beta_matrix)))
infinite_value_count <- sum(is.infinite(beta_matrix))
outside_range_count <- sum(beta_matrix < 0 | beta_matrix > 1, na.rm = TRUE)

min_beta_value <- min(beta_matrix, na.rm = TRUE)
max_beta_value <- max(beta_matrix, na.rm = TRUE)

# Save the beta matrix as an RDS file for downstream analysis
saveRDS(beta_matrix, "data/GSE87571/beta_matrix.rds")

# Save summary of the final matrix.
final_matrix_summary <- data.frame(
  number_of_probes = nrow(beta_matrix),
  number_of_samples = ncol(beta_matrix),
  missing_value_count = missing_value_count,
  duplicate_probe_id_count = duplicate_probe_count,
  duplicate_sample_id_count = duplicate_sample_count,
  infinite_value_count = infinite_value_count, # saftey, if 0 then beta matrix is valid
  values_outside_0_1_count = outside_range_count, # saftey, if 0 then beta matrix is valid
  minimum_beta_value = min_beta_value,
  maximum_beta_value = max_beta_value
)

write.csv(
  final_matrix_summary,
  "results/qc/final_matrix_summary.csv",
  row.names = FALSE
)

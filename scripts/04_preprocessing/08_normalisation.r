# Normalisation for GSE87571


# preprocessNoob is used here as a single fixed normalisation method.
# It performs background correction and dye-bias correction for Illumina
# methylation array data. This is simpler and less distribution-forcing than
# quantile normalisation, while still being a standard minfi preprocessing method.

library(minfi)

# Load the raw imported RGChannelSet.
# Since no samples failed QC, this is also the sample-QC-pass object.
rgSet <- readRDS("data/GSE87571/rgset_raw.rds")

# Apply Noob normalisation/background correction.
# The output is a MethylSet containing processed methylated and unmethylated
# signals.
mSet_noob <- preprocessNoob(rgSet)

# Save the normalised object for post-normalisation exploration and filtering.
saveRDS(mSet_noob, "data/GSE87571/mset_normalised.rds")

# Save a simple summary of the normalisation step.
normalisation_summary <- data.frame(
  method = "preprocessNoob",
  input_object = "data/GSE87571/rgset_raw.rds",
  output_object = "data/GSE87571/mset_normalised.rds",
  samples = ncol(mSet_noob),
  features = nrow(mSet_noob)
)

write.csv(
  normalisation_summary,
  "results/qc/normalisation_summary.csv",
  row.names = FALSE
)

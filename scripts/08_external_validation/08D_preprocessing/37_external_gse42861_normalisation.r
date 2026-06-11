# Normalisation for GSE42861

# preprocessNoob is used here as a single fixed normalisation method.
# It performs background correction and dye-bias correction for Illumina
# methylation array data. This is simpler and less distribution-forcing than
# quantile normalisation, while still being a standard minfi preprocessing method.

library(minfi)

# Load the sample-QC-pass RGChannelSet.
rgSet <- readRDS("data/GSE42861/rgset_qc_pass.rds")

# Apply Noob normalisation/background correction.
# The output is a MethylSet containing processed methylated and unmethylated
# signals.
mSet_noob <- preprocessNoob(rgSet)

# Save the normalised object for post-normalisation exploration and filtering.
saveRDS(mSet_noob, "data/GSE42861/mset_normalised.rds")

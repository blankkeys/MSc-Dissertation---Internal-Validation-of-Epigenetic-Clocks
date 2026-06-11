# Import IDAT files into R with minfi for GSE42861

# IDAT files are imported into an RGchannel set that stores the raw methylation data for each sample,
# including the red and green channel intensities

library(minfi) # package for analyzing methylation array data, including functions to read IDAT files and create RGchannel sets

# Read the matched normal-control sample sheet created in 08B
targets <- read.csv(
    "data/GSE42861/gse42861_qc_ready_sample_sheet_controls.csv",
    check.names = FALSE # prevent R from changing column names
)

# Read the IDAT files and create an RGchannel set
rgSet <- read.metharray.exp(
    targets = targets
)

# Save the RGchannel set as an RDS file for later
saveRDS(rgSet, "data/GSE42861/rgset_raw.rds")

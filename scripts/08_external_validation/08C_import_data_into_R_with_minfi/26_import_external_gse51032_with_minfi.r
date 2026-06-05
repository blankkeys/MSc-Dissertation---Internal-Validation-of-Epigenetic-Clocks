# Import IDAT files into R with minfi for GSE51032

# IDAT files are imported into an RGchannel set that stores the raw methylation data for each sample,
# including the red and green channel intensities

library(minfi) # package for analyzing methylation array data, including functions to read IDAT files and create RGchannel sets

# Read the matched cancer-free sample sheet created in 08B
targets <- read.csv(
    "data/GSE51032/gse51032_qc_ready_sample_sheet_cancer_free.csv",
    check.names = FALSE # prevent R from changing column names
)

# Read the IDAT files and create an RGchannel set
rgSet <- read.metharray.exp(
    targets = targets
)

# Save the RGchannel set as an RDS file for later
saveRDS(rgSet, "data/GSE51032/rgset_raw.rds")

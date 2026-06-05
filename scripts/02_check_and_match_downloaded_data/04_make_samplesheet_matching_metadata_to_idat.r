# Match metadata to IDAT files for GSE87571

# This script creates a simple sample sheet
# by matching metadata rows to IDAT files using the shared GSM accession.

library(readxl) # package to read Excel files

# file locations
metadata_file <- "data/GSE87571/GSE87571_additional_sample_characteristics.xlsx"
idat_file <- "data/GSE87571"

# Read metadata. The real column names start on row 3, so skip the first two rows.
metadata <- read_excel(metadata_file, sheet = 1, skip = 2)

# List one IDAT file per sample.
# Stage 3 already checked that every green file has a matching red file.
green_idat_files <- list.files(
  idat_file,
  pattern = "_Grn\\.idat\\.gz$",
  recursive = TRUE,
  full.names = TRUE
)

# Extract GSM IDs from IDAT filenames.
idat_sample_ids <- sub("_.*$", "", basename(green_idat_files))

# Create the Basename column needed by minfi.
# Basename is the IDAT path without the red/green file ending.
idat_basenames <- gsub("_Grn\\.idat\\.gz$", "", green_idat_files)

# Link GSM IDs to IDAT basenames.
idat_table <- data.frame(
  geo_accession = idat_sample_ids,
  Basename = idat_basenames
)

# Add the matching IDAT Basename to each metadata row.
matched_metadata <- merge(
  metadata,
  idat_table,
  by.x = "geo accession", # the column in the metadata file with GSM IDs
  by.y = "geo_accession", # the column in the IDAT table with GSM IDs
  all.x = TRUE, # keep all metadata rows even if they don't have a matching IDAT file
  sort = FALSE  # keep the original order of metadata rows, otherwise they will be sorted by GSM ID which is not ideal for checking the match
)

# Check whether any metadata samples failed to match an IDAT file.
sum(is.na(matched_metadata$Basename)) # if this is greater than 0, it means some metadata rows did not have a matching IDAT file, which could indicate a problem with the matching process or missing data

# Save the matched sample sheet for minfi import.
write.csv(
  matched_metadata,
  "data/GSE87571/qc_ready_sample_sheet.csv",
  row.names = FALSE # don't include row names in the output CSV file, otherwise minfi will get confused when reading the sample sheet
)

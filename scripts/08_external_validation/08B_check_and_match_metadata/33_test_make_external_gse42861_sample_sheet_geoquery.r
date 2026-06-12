# Match GEOquery metadata to IDAT files for GSE42861

# This script tests whether GEOquery can make the same sample sheet more simply
# by reading the local GEO series matrix as a metadata table

library(GEOquery)

# file locations
series_matrix_file <- "data/GSE42861/GSE42861_series_matrix.txt.gz"
idat_file <- "data/GSE42861/GSE42861_RAW"

# read GEO metadata from the local series matrix file
gse <- getGEO(filename = series_matrix_file)
metadata <- pData(gse)

# create a simple metadata table
metadata <- data.frame(
  geo_accession = rownames(metadata),
  Sample_Name = metadata$title,
  source_name = metadata$source_name_ch1,
  age = as.numeric(metadata$age.ch1),
  sex = tolower(metadata$gender.ch1),
  disease_status = metadata$`disease status:ch1`
)

# identify normal controls
metadata$control_status <- ifelse(
  grepl("^Normal|Normal genomic DNA", metadata$Sample_Name, ignore.case = TRUE) |
    grepl("normal|control", metadata$disease_status, ignore.case = TRUE),
  "normal_control",
  "case_or_other"
)

# list one IDAT file per sample
green_idat_files <- list.files(
  idat_file,
  pattern = "_Grn\\.idat\\.gz$",
  recursive = TRUE,
  full.names = TRUE
)

# extract GSM IDs from IDAT filenames
idat_sample_ids <- sub("_.*$", "", basename(green_idat_files))

# create the Basename column needed by minfi
idat_basenames <- gsub("_Grn\\.idat\\.gz$", "", green_idat_files)

# link GSM IDs to IDAT basenames
idat_table <- data.frame(
  geo_accession = idat_sample_ids,
  Basename = idat_basenames
)

# add the matching IDAT Basename to each metadata row
matched_metadata <- merge(
  metadata,
  idat_table,
  by = "geo_accession",
  all.x = TRUE,
  sort = FALSE
)

# keep normal controls for external validation
control_metadata <- matched_metadata[
  matched_metadata$control_status == "normal_control" &
    !is.na(matched_metadata$age),
]

# save the GEOquery-derived sample sheet for comparison
write.csv(
  control_metadata,
  "data/GSE42861/gse42861_qc_ready_sample_sheet_controls_geoquery.csv",
  row.names = FALSE
)

dim(control_metadata)

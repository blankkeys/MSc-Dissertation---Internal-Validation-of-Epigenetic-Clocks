# Match metadata to IDAT files for GSE51032

# This script creates a sample sheet for external validation
# by matching GEO series matrix metadata to IDAT files using the shared GSM accession

series_matrix_file <- "data/GSE51032/GSE51032_series_matrix.txt.gz"
idat_dir <- "data/GSE51032/GSE51032_RAW"

if (!file.exists(series_matrix_file)) {
  stop("The GSE51032 series matrix file was not found")
}

if (!dir.exists(idat_dir)) {
  stop("The GSE51032 IDAT directory was not found")
}

# Read the GEO series matrix metadata
series_lines <- readLines(gzfile(series_matrix_file))

# Extract one sample metadata row from the series matrix
extract_sample_row <- function(row) {
  if (is.na(row)) {
    stop("A required GSE51032 metadata row was not found")
  }

  values <- strsplit(row, "\t")[[1]][-1]
  gsub('"', "", values)
}

geo_accession <- extract_sample_row(series_lines[grepl("^!Sample_geo_accession\t", series_lines)][1])
sample_title <- extract_sample_row(series_lines[grepl("^!Sample_title\t", series_lines)][1])
source_name <- extract_sample_row(series_lines[grepl("^!Sample_source_name_ch1\t", series_lines)][1])

characteristic_rows <- series_lines[grepl("^!Sample_characteristics_ch1\t", series_lines)]

gender <- extract_sample_row(characteristic_rows[grepl("gender:|sex:", characteristic_rows, ignore.case = TRUE)][1])
age <- extract_sample_row(characteristic_rows[grepl("age:", characteristic_rows, ignore.case = TRUE)][1])
cancer_characteristic <- extract_sample_row(characteristic_rows[grepl("cancer", characteristic_rows, ignore.case = TRUE)][1])
time_to_diagnosis <- extract_sample_row(characteristic_rows[grepl("time to diagnosis", characteristic_rows, ignore.case = TRUE)][1])

# Clean metadata values needed for external validation
gender <- tolower(trimws(gsub("gender:|sex:", "", gender, ignore.case = TRUE)))
age <- as.numeric(trimws(gsub("age:", "", age, ignore.case = TRUE)))
cancer_characteristic <- trimws(gsub("characteristics_ch1:", "", cancer_characteristic))
time_to_diagnosis <- trimws(gsub("time to diagnosis:", "", time_to_diagnosis, ignore.case = TRUE))

metadata <- data.frame(
  geo_accession = geo_accession,
  Sample_Name = sample_title,
  source_name = source_name,
  sex = gender,
  age = age,
  cancer_characteristic = cancer_characteristic,
  time_to_diagnosis = time_to_diagnosis
)

# Samples without a time to diagnosis are treated as the cancer-free external validation cohort
metadata$cancer_status <- ifelse(
  metadata$time_to_diagnosis == "",
  "cancer_free",
  "cancer_case"
)

# List one IDAT file per sample
# Stage 24 already checked that every green file has a matching red file
green_idat_files <- list.files(
  idat_dir,
  pattern = "_Grn\\.idat\\.gz$",
  recursive = TRUE,
  full.names = TRUE
)

# Extract GSM IDs from IDAT filenames
idat_sample_ids <- sub("_.*$", "", basename(green_idat_files))

# Create the Basename column needed by minfi
# Basename is the IDAT path without the red/green file ending
idat_basenames <- gsub("_Grn\\.idat\\.gz$", "", green_idat_files)

# Link GSM IDs to IDAT basenames
idat_table <- data.frame(
  geo_accession = idat_sample_ids,
  Basename = idat_basenames
)

# Add the matching IDAT Basename to each metadata row
matched_metadata <- merge(
  metadata,
  idat_table,
  by = "geo_accession",
  all.x = TRUE,
  sort = FALSE
)

# Check whether any metadata samples failed to match an IDAT file
if (sum(is.na(matched_metadata$Basename)) > 0) {
  stop("Some GSE51032 metadata samples did not match an IDAT file")
}

# Keep cancer-free samples for the primary external validation cohort
cancer_free_metadata <- matched_metadata[
  matched_metadata$cancer_status == "cancer_free",
]

if (nrow(cancer_free_metadata) == 0) {
  stop("No cancer-free GSE51032 samples were identified")
}

# Save the matched sample sheets for minfi import
write.csv(
  matched_metadata,
  "data/GSE51032/gse51032_qc_ready_sample_sheet_all_samples.csv",
  row.names = FALSE
)

write.csv(
  cancer_free_metadata,
  "data/GSE51032/gse51032_qc_ready_sample_sheet_cancer_free.csv",
  row.names = FALSE
)

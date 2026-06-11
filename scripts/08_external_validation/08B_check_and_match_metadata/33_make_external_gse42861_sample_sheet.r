# Match metadata to IDAT files for GSE42861

# This script creates a sample sheet for external validation
# by matching GEO series matrix metadata to IDAT files using the shared GSM accession

# file locations
series_matrix_file <- "data/GSE42861/GSE42861_series_matrix.txt.gz"
idat_file <- "data/GSE42861/GSE42861_RAW"

if (!file.exists(series_matrix_file)) {
  stop("The GSE42861 series matrix file was not found")
}

if (!dir.exists(idat_file)) {
  stop("The GSE42861 IDAT directory was not found")
}

# read the GEO series matrix metadata
series_lines <- readLines(gzfile(series_matrix_file))

# extract one sample metadata row from the series matrix
extract_sample_row <- function(row) {
  if (is.na(row)) {
    stop("A required GSE42861 metadata row was not found")
  }

  values <- strsplit(row, "\t")[[1]][-1]
  gsub('"', "", values)
}

extract_optional_sample_row <- function(row, sample_count) {
  if (is.na(row)) {
    return(rep(NA, sample_count))
  }

  extract_sample_row(row)
}

# extract sample IDs and metadata rows from the series matrix
geo_accession <- extract_sample_row(series_lines[grepl("^!Sample_geo_accession\t", series_lines)][1])
sample_title <- extract_sample_row(series_lines[grepl("^!Sample_title\t", series_lines)][1])
source_name <- extract_sample_row(series_lines[grepl("^!Sample_source_name_ch1\t", series_lines)][1])

characteristic_rows <- series_lines[grepl("^!Sample_characteristics_ch1\t", series_lines)]

gender <- extract_optional_sample_row(
  characteristic_rows[grepl("gender:|sex:", characteristic_rows, ignore.case = TRUE)][1],
  length(geo_accession)
)
age <- extract_optional_sample_row(
  characteristic_rows[grepl("age:", characteristic_rows, ignore.case = TRUE)][1],
  length(geo_accession)
)
disease_status <- extract_optional_sample_row(
  characteristic_rows[grepl("disease|diagnosis|status|phenotype|case|control", characteristic_rows, ignore.case = TRUE)][1],
  length(geo_accession)
)

# clean metadata values needed for external validation
gender <- tolower(trimws(gsub("gender:|sex:", "", gender, ignore.case = TRUE)))
age <- as.numeric(trimws(gsub("age:", "", age, ignore.case = TRUE)))
disease_status <- trimws(gsub("disease state:|diagnosis:|status:|phenotype:", "", disease_status, ignore.case = TRUE))

metadata <- data.frame(
  geo_accession = geo_accession,
  Sample_Name = sample_title,
  source_name = source_name,
  sex = gender,
  age = age,
  disease_status = disease_status
)

# GSE42861 titles identify rheumatoid arthritis patients and normal controls
metadata$control_status <- ifelse(
  grepl("^Normal|Normal genomic DNA", metadata$Sample_Name, ignore.case = TRUE) |
    grepl("normal|control", metadata$disease_status, ignore.case = TRUE),
  "normal_control",
  "case_or_other"
)

# list one IDAT file per sample
# Stage 32 already checked that every green file has a matching red file
green_idat_files <- list.files(
  idat_file,
  pattern = "_Grn\\.idat\\.gz$",
  recursive = TRUE,
  full.names = TRUE
)

# extract GSM IDs from IDAT filenames
idat_sample_ids <- sub("_.*$", "", basename(green_idat_files))

# create the Basename column needed by minfi
# Basename is the IDAT path without the red/green file ending
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

# check whether any metadata samples failed to match an IDAT file
if (any(is.na(matched_metadata$Basename))) {
  stop("Some GSE42861 metadata samples did not match an IDAT file")
}

# keep normal controls for the primary external validation cohort
control_metadata <- matched_metadata[
  matched_metadata$control_status == "normal_control" &
    !is.na(matched_metadata$age),
]

if (nrow(control_metadata) == 0) {
  stop("No GSE42861 normal controls with age metadata were identified")
}

dir.create("data/GSE42861", recursive = TRUE, showWarnings = FALSE)
dir.create("results/external_validation/qc", recursive = TRUE, showWarnings = FALSE)

# save the matched sample sheets for minfi import
write.csv(
  matched_metadata,
  "data/GSE42861/gse42861_qc_ready_sample_sheet_all_samples.csv",
  row.names = FALSE
)

write.csv(
  control_metadata,
  "data/GSE42861/gse42861_qc_ready_sample_sheet_controls.csv",
  row.names = FALSE
)

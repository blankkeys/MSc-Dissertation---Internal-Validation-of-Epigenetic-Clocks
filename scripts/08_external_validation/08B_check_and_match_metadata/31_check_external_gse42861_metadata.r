# check metadata file for GSE42861

# GSE42861 is being used for external validation
# the series matrix contains information about the samples like:
# ID, age, sex, and disease/control status
# need to check whether R can read the metadata file before matching metadata to IDAT files

# file location
series_matrix_file <- "data/GSE42861/GSE42861_series_matrix.txt.gz"

if (!file.exists(series_matrix_file)) {
  stop("The GSE42861 series matrix file was not found")
}

# read the GEO series matrix
series_lines <- readLines(gzfile(series_matrix_file))

# keep useful sample metadata lines
# these lines will be checked before making the sample sheet
metadata_lines <- series_lines[
  grepl("^!Sample_", series_lines) &
    grepl(
      "age|gender|sex|disease|diagnosis|case|control|status|title|geo_accession|characteristics",
      series_lines,
      ignore.case = TRUE
    )
]

dir.create("results/external_validation", recursive = TRUE, showWarnings = FALSE)

# save useful metadata lines so they can be inspected if needed
writeLines(
  metadata_lines,
  "results/external_validation/gse42861_series_metadata_lines.txt"
)

# show how many useful metadata lines were found
length(metadata_lines)

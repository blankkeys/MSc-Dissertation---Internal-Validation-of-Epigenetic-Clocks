# check downloaded metadata for GSE42861 external validation

# GSE42861 is being used for external validation
# This metadata will be used to identify age and normal-control samples

series_matrix_file <- "data/GSE42861/GSE42861_series_matrix.txt.gz"

dir.create("results/external_validation", recursive = TRUE, showWarnings = FALSE)

if (!file.exists(series_matrix_file)) {
  stop("The GSE42861 series matrix file was not found")
}

# read useful metadata lines from the GEO series matrix
series_lines <- readLines(gzfile(series_matrix_file))

metadata_lines <- series_lines[
  grepl("^!Sample_", series_lines) &
    grepl(
      "age|gender|sex|disease|diagnosis|case|control|status|title|geo_accession|characteristics",
      series_lines,
      ignore.case = TRUE
    )
]

writeLines(
  metadata_lines,
  "results/external_validation/gse42861_series_metadata_lines.txt"
)

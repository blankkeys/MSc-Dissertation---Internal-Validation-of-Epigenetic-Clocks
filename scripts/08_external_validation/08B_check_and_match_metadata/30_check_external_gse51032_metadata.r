# check downloaded metadata for GSE51032 external validation

# GSE51032 is being used for external validation
# This metadata will be used to identify age and cancer-free samples

series_matrix_file <- "data/GSE51032/GSE51032_series_matrix.txt.gz"

dir.create("results/external_validation", recursive = TRUE, showWarnings = FALSE)

if (!file.exists(series_matrix_file)) {
  stop("The GSE51032 series matrix file was not found")
}

# read useful metadata lines from the GEO series matrix
series_lines <- readLines(gzfile(series_matrix_file))

metadata_lines <- series_lines[
  grepl("^!Sample_", series_lines) &
    grepl(
      "age|gender|sex|cancer|case|control|status|title|geo_accession|characteristics",
      series_lines,
      ignore.case = TRUE
    )
]

writeLines(
  metadata_lines,
  "results/external_validation/gse51032_series_metadata_lines.txt"
)

# urls to download GEO beta matrix and metadata for GSE40279 external validation

options(timeout = 3600) # increase timeout for the large beta matrix file

dir.create("data/GSE40279", recursive = TRUE, showWarnings = FALSE)

beta_matrix_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/suppl/GSE40279_average_beta.txt.gz"
sample_key_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/suppl/GSE40279_sample_key.txt.gz"
series_matrix_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/matrix/GSE40279_series_matrix.txt.gz"

# download processed beta matrix
if (!file.exists("data/GSE40279/GSE40279_average_beta.txt.gz")) {
  download.file(beta_matrix_url,
    destfile = "data/GSE40279/GSE40279_average_beta.txt.gz",
    method = "auto",
    mode = "wb"
  )
}

# download sample key metadata file
if (!file.exists("data/GSE40279/GSE40279_sample_key.txt.gz")) {
  download.file(sample_key_url,
    destfile = "data/GSE40279/GSE40279_sample_key.txt.gz",
    method = "auto",
    mode = "wb"
  )
}

# download GEO series matrix metadata file
if (!file.exists("data/GSE40279/GSE40279_series_matrix.txt.gz")) {
  download.file(series_matrix_url,
    destfile = "data/GSE40279/GSE40279_series_matrix.txt.gz",
    method = "auto",
    mode = "wb"
  )
}

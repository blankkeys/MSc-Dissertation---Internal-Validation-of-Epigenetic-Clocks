# urls to download GEO IDAT and metadata for GSE51032 external validation

options(timeout = 7200) # large raw IDAT archive needs a longer timeout

dir.create("data/GSE51032", recursive = TRUE, showWarnings = FALSE)

idat_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE51032&format=file"
series_matrix_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE51nnn/GSE51032/matrix/GSE51032_series_matrix.txt.gz"

# download IDAT files
if (!file.exists("data/GSE51032/GSE51032_RAW.tar")) {
  download.file(idat_url,
  destfile = "data/GSE51032/GSE51032_RAW.tar",
  method = "auto",
  mode = "wb")
}

# download GEO series matrix metadata file
if (!file.exists("data/GSE51032/GSE51032_series_matrix.txt.gz")) {
  download.file(series_matrix_url,
  destfile = "data/GSE51032/GSE51032_series_matrix.txt.gz",
  method = "auto",
  mode = "wb")
}

# extract IDAT data
if (!dir.exists("data/GSE51032/GSE51032_RAW")) {
  untar("data/GSE51032/GSE51032_RAW.tar", exdir = "data/GSE51032/GSE51032_RAW")
}

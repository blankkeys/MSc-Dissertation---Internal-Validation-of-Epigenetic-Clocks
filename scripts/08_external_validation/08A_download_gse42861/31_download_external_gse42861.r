# urls to download GEO IDAT and metadata for GSE42861 external validation

options(timeout = 7200) # large raw IDAT archive needs a longer timeout

dir.create("data/GSE42861", recursive = TRUE, showWarnings = FALSE)

idat_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE42861&format=file"
series_matrix_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42861/matrix/GSE42861_series_matrix.txt.gz"

# download IDAT files
if (!file.exists("data/GSE42861/GSE42861_RAW.tar")) {
  download.file(idat_url,
  destfile = "data/GSE42861/GSE42861_RAW.tar",
  method = "auto",
  mode = "wb")
}

# download GEO series matrix metadata file
if (!file.exists("data/GSE42861/GSE42861_series_matrix.txt.gz")) {
  download.file(series_matrix_url,
  destfile = "data/GSE42861/GSE42861_series_matrix.txt.gz",
  method = "auto",
  mode = "wb")
}

# extract IDAT data
if (!dir.exists("data/GSE42861/GSE42861_RAW")) {
  untar("data/GSE42861/GSE42861_RAW.tar", exdir = "data/GSE42861/GSE42861_RAW")
}

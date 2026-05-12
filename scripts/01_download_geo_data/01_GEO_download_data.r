# urls to download GEO IDAT and metadata for GSE87571

options(timeout = 3600)

idat_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE87571&format=file"
metadata_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE87571&format=file&file=GSE87571%5Fadditional%5Fsample%5Fchararcteristics%2Exlsx"

# download IDAT files
if (!file.exists("GSE87571_RAW.tar")) {
  download.file(idat_url, destfile = "GSE87571_RAW.tar", method = "auto", mode = "wb")
}

# download metadata Excel file
if (!file.exists("GSE87571_additional_sample_characteristics.xlsx")) {
  download.file(metadata_url, destfile = "GSE87571_additional_sample_characteristics.xlsx", method = "auto", mode = "wb")
}

# extract IDAT data
if (!dir.exists("GSE87571_RAW")) {
  untar("GSE87571_RAW.tar", exdir = "GSE87571_RAW")
}

# location of IDAT files and metadata file
idat_directory <- "GSE87571_RAW"
metadata_file <- "GSE87571_additional_sample_characteristics.xlsx"

# list IDAT files
idat_files <- list.files(
  idat_directory,
  recursive = TRUE
)

idat_files

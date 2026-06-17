# urls to download GEO IDAT and metadata for GSE42861

options(timeout = 7200) # 7200s as external raw archive is large

idat_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE42861&format=file"
series_matrix_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42861/matrix/GSE42861_series_matrix.txt.gz"

dir.create("data/GSE42861", recursive = TRUE, showWarnings = FALSE)

# download IDAT files
if (!file.exists("data/GSE42861/GSE42861_RAW.tar")) {
  download.file(idat_url, 
  destfile = "data/GSE42861/GSE42861_RAW.tar", # name of the file to save the downloaded data as
  method = "auto", # method to download the file, "auto" will choose the best method available on the system
  mode = "wb") # mode = "wb" is used to write the file in binary mode, which is important for non-text files like tar archives
}

#download series matrix file
if (!file.exists("data/GSE42861/GSE42861_series_matrix.txt.gz")) {
  download.file(series_matrix_url, 
  destfile = "data/GSE42861/GSE42861_series_matrix.txt.gz", 
  method = "auto", 
  mode = "wb") 
}


# extract IDAT data
if (!dir.exists("data/GSE42861/GSE42861_RAW")) {
  untar("data/GSE42861/GSE42861_RAW.tar", exdir = "data/GSE42861/GSE42861_RAW") # exdir specifies the directory to extract the files to, in this case "GSE42861_RAW"
}


# urls to download GEO IDAT and metadata for GSE87571

options(timeout = 3600) # 3600s as 60s caused timeout

idat_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE87571&format=file"
metadata_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE87571&format=file&file=GSE87571%5Fadditional%5Fsample%5Fchararcteristics%2Exlsx"
series_matrix_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/GSE87571/matrix/GSE87571_series_matrix.txt.gz"
# download IDAT files
if (!file.exists("data/GSE87571/GSE87571_RAW.tar")) {
  download.file(idat_url, 
  destfile = "data/GSE87571/GSE87571_RAW.tar", # name of the file to save the downloaded data as
  method = "auto", # method to download the file, "auto" will choose the best method available on the system
  mode = "wb") # mode = "wb" is used to write the file in binary mode, which is important for non-text files like tar archives
}

# download metadata Excel file
if (!file.exists("data/GSE87571/GSE87571_additional_sample_characteristics.xlsx")) {
  download.file(metadata_url, 
  destfile = "data/GSE87571/GSE87571_additional_sample_characteristics.xlsx", 
  method = "auto", 
  mode = "wb") # mode = "wb" is used to write the file in binary mode, which is important for non-text files like Excel files
}

#download series matrix file
if (!file.exists("data/GSE87571/GSE87571_series_matrix.txt.gz")) {
  download.file(series_matrix_url, 
  destfile = "data/GSE87571/GSE87571_series_matrix.txt.gz", 
  method = "auto", 
  mode = "wb") 
}


# extract IDAT data
if (!dir.exists("data/GSE87571/GSE87571_RAW")) {
  untar("data/GSE87571/GSE87571_RAW.tar", exdir = "data/GSE87571/GSE87571_RAW") # exdir specifies the directory to extract the files to, in this case "GSE87571_RAW"
}



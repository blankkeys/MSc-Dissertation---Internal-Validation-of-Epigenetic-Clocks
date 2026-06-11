# check IDAT files for GSE42861

# GSE42861 uses the 450K platform and provides raw IDAT files
# each methylation array sample should have two IDAT files:
# red and green channel files
# this script checks that both channels are present as minfi requires both to work

# file locations
raw_tar_file <- "data/GSE42861/GSE42861_RAW.tar"
idat_file <- "data/GSE42861/GSE42861_RAW"

if (!file.exists(raw_tar_file)) {
  stop("The GSE42861 raw IDAT archive was not found")
}

if (!dir.exists(idat_file)) {
  stop("The GSE42861 IDAT directory was not found")
}

# list green and red IDAT files
green_idat_files <- list.files(
  idat_file,
  pattern = "_Grn\\.idat\\.gz$",
  recursive = TRUE, # read files in subdirectories too
  full.names = TRUE # include the folder path in the file names
)
red_idat_files <- list.files(
  idat_file,
  pattern = "_Red\\.idat\\.gz$",
  recursive = TRUE,
  full.names = TRUE
)

# number of green and red IDAT files should be the same
length(green_idat_files)
length(red_idat_files)

# remove colour channel endings to get sample base names to compare red and green files
green_sample_names <- gsub("_Grn\\.idat\\.gz", "", basename(green_idat_files))
red_sample_names <- gsub("_Red\\.idat\\.gz", "", basename(red_idat_files))

# check which red/green files are missing by comparing sample names
missing_green <- setdiff(red_sample_names, green_sample_names)
missing_red <- setdiff(green_sample_names, red_sample_names)

if (length(missing_red) > 0 || length(missing_green) > 0) {
  stop("Some GSE42861 IDAT files do not have matching red/green pairs")
}

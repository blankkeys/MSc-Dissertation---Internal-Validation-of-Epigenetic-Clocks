#check IDAT file for GSE87571

# each methylation array sample should have two IDAT files:
# red and green channels files
# this script checks that both channels are present as minfi requires both to work 

# file location
idat_file <- "data/GSE87571"

# list green and red IDAT files
green_idat_files <- list.files(
  idat_file,
  pattern = "_Grn\\.idat\\.gz$",
  recursive = TRUE, # read files in subdirectories too, if FALSE will just read given directory
  full.names = TRUE # include the folder path in the file names, if FALSE will just give file names without path
)
red_idat_files <- list.files(
  idat_file,
  pattern = "_Red\\.idat\\.gz$", # regex pattern to match red channel IDAT files, the $ at the end ensures it matches the end of the file name, so it won't match files that have additional characters after "Red.idat.gz"
  recursive = TRUE,
  full.names = TRUE
)

# number of green and red IDAT files should be the same (checking no. of each)
length(green_idat_files)
length(red_idat_files)


# remove colour channels endings to get sample base names to compare 'red' and 'green' files 
green_sample_names <- gsub("_Grn\\.idat\\.gz", "", basename(green_idat_files))
red_sample_names <- gsub("_Red\\.idat\\.gz", "", basename(red_idat_files))

# check which red/green files are missing by comparing sample names of red and green files
missing_green <- setdiff(red_sample_names, green_sample_names)
missing_red <- setdiff(green_sample_names, red_sample_names)

# Stop if any IDAT pairs are incomplete.
if (length(missing_red) > 0 || length(missing_green) > 0) {
  stop("Some IDAT files do not have matching red/green pairs.")
}

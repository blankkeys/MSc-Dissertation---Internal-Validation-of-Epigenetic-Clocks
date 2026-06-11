# check downloaded IDAT files for GSE42861 external validation

# GSE42861 uses the 450K platform and provides raw IDAT files
# This check confirms red and green IDAT files are available as matched pairs

raw_tar_file <- "data/GSE42861/GSE42861_RAW.tar"
idat_dir <- "data/GSE42861/GSE42861_RAW"

dir.create("results/external_validation", recursive = TRUE, showWarnings = FALSE)

if (!file.exists(raw_tar_file)) {
  stop("The GSE42861 raw IDAT archive was not found")
}

if (!dir.exists(idat_dir)) {
  stop("The GSE42861 IDAT directory was not found")
}

# list red and green IDAT files
idat_files <- list.files(idat_dir, pattern = "idat.gz$", recursive = TRUE)
green_idats <- idat_files[grepl("_Grn.idat.gz$", idat_files)]
red_idats <- idat_files[grepl("_Red.idat.gz$", idat_files)]

if (length(green_idats) == 0 || length(red_idats) == 0) {
  stop("No GSE42861 red/green IDAT files were found")
}

# check that red and green IDAT files form matched pairs
green_basenames <- sub("_Grn.idat.gz$", "", basename(green_idats))
red_basenames <- sub("_Red.idat.gz$", "", basename(red_idats))

missing_red <- setdiff(green_basenames, red_basenames)
missing_green <- setdiff(red_basenames, green_basenames)

if (length(missing_red) > 0 || length(missing_green) > 0) {
  stop("Some GSE42861 IDAT files do not have matching red/green pairs")
}

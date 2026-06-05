# check meta data file for GSE87571

# meta data file contains information about the samples like:
# ID, age, sex
# need to check wether R can read the metadata file before mathcing metadata to IDAT files

library(readxl) #package to read Excel files

# file location
metadata_file <- "data/GSE87571/GSE87571_additional_sample_characteristics.xlsx"

# show sheet names in the metadata file 
excel_sheets(metadata_file)

# read the first sheet of the metadata file
# skip = 2 to start on row 3 where GEO accession numbers start,
# otherwise the first two rows of the file will be read as metadata which is not correct
metadata <- read_excel(metadata_file, sheet = 1, skip = 2)

# look at the metadata
head(metadata)

# show column names 
#important because we will need to match the sample IDs in the metadata to the sample IDs in the IDAT files
names(metadata)

#show size of metadata (no.rows, no.columns) 
#to make sure everything expected to be read was read
dim(metadata)

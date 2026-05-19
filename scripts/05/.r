# Prepare modelling metadata for GSE87571
# This script extracts age from the GEO series matrix and 
# matches it to the final beta matrix samples.
library(GEOquery) # for reading GEO series matrix files

# Load the GEO series matrix file for GSE87571
gse <- getGEO(
    filename = "data/GSE87571/GSE87571_series_matrix.txt.gz",
)

#extract sample meta data
metadata <- pData(gse)

#only keep GEO sample ID and age
age_metadata <- data.frame(
    geo_accession = rownames(metadata),
    age = metadata$characteristics_ch1
)

#clean values so instead of age: 72 we just get 72
age_metadata$age <- gsub("age: ", "", age_metadata$age) #rmv age: 
age_metadata$age <- as.numeric(age_metadata$age) #convert to numeric

#save age as metadata
write.csv(
    age_metadata,
    "data/GSE87571/age_metadata.csv",
    row.names = FALSE  #FALSE to avoid writing into the csv file
)
# Download updated HM450 probe annotation from Zhou/InfiniumAnnotation.
# This is used later for annotation-based probe filtering.
# It is not methylation data; it is a platform annotation resource.

annotation_dir <- "data/annotation"

dir.create(annotation_dir, recursive = TRUE, showWarnings = FALSE)

zhou_annotation_url <- "https://zhouserver.research.chop.edu/InfiniumAnnotation/current/HM450/HM450.hg19.manifest.tsv.gz"

zhou_annotation_file <- file.path(
  annotation_dir,
  "HM450.hg19.manifest.tsv.gz"
)

if (!file.exists(zhou_annotation_file)) {
  download.file(
    zhou_annotation_url,
    destfile = zhou_annotation_file,
    method = "auto",
    mode = "wb"
  )
}
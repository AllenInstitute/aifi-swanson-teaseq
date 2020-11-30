installed <- .packages(all.available = TRUE)

cran_dep <- c(
  "BiocManager",
  "cowplot",
  "data.table",
  "dendextend",
  "dplyr",
  "ggplot2",
  "ggbeeswarm",
  "ggrastr",
  "jsonlite",
  "Matrix",
  "optparse",
  "purrr",
  "readr"
  "Signac",
  "uwot"
)

cran_missing <- setdiff(cran_dep, installed)

install.packages(cran_missing)

bioc_dep <- c(
  "GenomicRanges",
  "rhdf5",
  "S4Vectors"
)

bioc_missing <- setdiff(bioc_dep, installed)

BiocManager::install(bioc_missing)

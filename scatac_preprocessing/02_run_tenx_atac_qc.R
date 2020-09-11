library(optparse)

option_list <- list(
  make_option(opt_str = c("-p","--in_pre"),
              type = "character",
              default = NULL,
              help = "Input ATAC preprocessed directory",
              metavar = "character"),
  make_option(opt_str = c("-s","--in_sin"),
              type = "character",
              default = NULL,
              help = "Input singlecell.csv file",
              metavar = "character"),
  make_option(opt_str = c("-u", "--in_sum"),
              type = "character",
              default = NULL,
              help = "Input summary.csv file",
              metavar = "character"),
  make_option(opt_str = c("-k", "--in_key"),
              type = "character",
              default = NULL,
              help = "Input SampleSheet.csv file",
              metavar = "character"),
  make_option(opt_str = c("-w","--in_well"),
              type = "character",
              default = NULL,
              help = "Well",
              metavar = "character"),
  make_option(opt_str = c("-g","--genome"),
              type = "character",
              default = "hg38",
              help = "Genome (hg38 or hg19)",
              metavar = "character"),
  make_option(opt_str = c("-d","--out_dir"),
              type = "character",
              default = NULL,
              help = "Output directory",
              metavar = "character"),
  make_option(opt_str = c("-o","--out_html"),
              type = "character",
              default = NULL,
              help = "Output HTML run summary file",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)

args <- parse_args(opt_parser)

if(is.null(args$out_html)) {
  print_help(opt_parser)
  stop("No parameters supplied.")
}

rmd_loc <- "rmarkdown/tenx_atac_qc.Rmd"

rmarkdown::render(
  input = rmd_loc,
  params = list(in_pre = args$in_pre,
                in_sin = args$in_sin,
                in_sum = args$in_sum,
                in_key = args$in_key,
                in_well = args$in_well,
                genome = args$genome,
                out_dir = args$out_dir),
  output_file = args$out_html,
  quiet = TRUE
)

file.remove(rmd_loc)

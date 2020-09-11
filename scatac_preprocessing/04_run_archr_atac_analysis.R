library(optparse)

option_list <- list(
  make_option(opt_str = c("-f","--in_frag"),
              type = "character",
              default = NULL,
              help = "Input filtered fragments.tsv.gz",
              metavar = "character"),
  make_option(opt_str = c("-s", "--in_sample"),
              type = "character",
              default = NULL,
              help = "Input SampleID",
              metavar = "character"),
  make_option(opt_str = c("-g","--genome"),
              type = "character",
              default = "hg38",
              help = "Genome (hg38 or hg19)",
              metavar = "character"),
  make_option(opt_str = c("-t","--n_threads"),
              type = "character",
              default = NULL,
              help = "N Parallel Threads",
              metavar = "numeric"),
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

rmd_loc <- "rmarkdown/tenx_atac_archr_analysis.Rmd"

rmarkdown::render(
  input = rmd_loc,
  params = list(in_frag = args$in_frag,
                in_sample = args$in_sample,
                genome = args$genome,
                n_threads = args$n_threads,
                out_dir = args$out_dir),
  output_file = args$out_html,
  quiet = TRUE
)

file.remove(rmd_loc)

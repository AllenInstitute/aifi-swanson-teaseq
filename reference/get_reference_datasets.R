library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(tidyr)
library(Seurat)
source("../common/hdf5_helpers.R")

## UCSC Genome Browser chrom.sizes files:
genomes <- c("hg38")

for(genome in genomes) {
  out_file <- paste0(genome, ".chrom.sizes")
  download.file(paste0("http://hgdownload.soe.ucsc.edu/goldenPath/",genome,"/bigZips/",genome,".chrom.sizes"),
                out_file)
}

## UCSC Genome Browser chain files:
download.file("http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
              "hg38ToHg19.over.chain.gz")
download.file("http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
              "hg19ToHg38.over.chain.gz")

## Altius DHS Index from ENCODE:
temp_file <- tempfile(fileext = ".tsv")
download.file("https://www.encodeproject.org/files/ENCFF503GCK/@@download/ENCFF503GCK.tsv",
              temp_file)

alt_idx_raw <- fread(temp_file)

alt_idx_bed <- alt_idx_raw[,c("seqname","start","end","identifier")]
fwrite(alt_idx_bed,
       "hg38_altius.bed.gz",
       sep = "\t",
       quote = FALSE,
       col.names = FALSE,
       row.names = FALSE)

alt_idx_gr <- GRanges(seqnames = alt_idx_raw$seqname,
                      IRanges(start = alt_idx_raw$start,
                              end = alt_idx_raw$end),
                      identifier = alt_idx_raw$identifier)

saveRDS(alt_idx_gr,
        "inst/reference/hg38_altius_gr.rds")

## Buenrostro/GSE123577 peak set
temp_chain <- tempfile(fileext = ".over.chain")
R.utils::gunzip("hg19ToHg38.over.chain.gz",
                temp_chain,
                remove = FALSE)

ch_to_19 <- import.chain(temp_chain)

temp_file <- tempfile(fileext = ".bed.gz")
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE123577&format=file&file=GSE123577%5Fpbmc%5Fpeaks%2Ebed%2Egz",
              temp_file)

hg19_peaks <- fread(temp_file)
names(hg19_peaks) <- c("chr","start","end")
hg19_gr <- convert_fragments_gr(list(hg19_peaks))[[1]]

hg38_lo <- liftOver(hg19_gr, ch_to_19)
n_lo <- elementNROWS(hg38_lo)

keep_lo <- n_lo == 1

hg38_gr <- unlist(hg38_lo[keep_lo])
hg38_gr <- sort(hg38_gr, ignore.strand = TRUE)

saveRDS(hg38_gr, "inst/reference/hg38_peaks_gr.rds")

hg38_bed <- as.data.frame(hg38_gr)[,1:3]
fwrite(hg38_bed,
       "inst/reference/hg38_peaks.bed.gz",
       sep = "\t",
       quote = FALSE,
       col.names = FALSE,
       row.names = FALSE)

## ENSEMBL Hg38 v93
## Couldn't download directly, so downloaded through browser from:
## ftp://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/

hg38_chrom_sizes <- read_chrom_sizes("hg38")

gtf <- fread("C:/Users/lucasg/Downloads/Homo_sapiens.GRCh38.93.gtf.gz")
names(gtf) <- c("seqname","source","feature","start","end","score","strand","frame","attribute")
gtf <- gtf[feature == "gene"]

gtf$attribute <- gsub(" ?[a-z|_]+ \"([A-Za-z0-9|_|-]+)\"",
                      "\\1",
                      gtf$attribute)

gtf <- separate(gtf,
                attribute,
                sep = ";",
                into = c("gene_id","gene_version","gene_name","gene_source","gene_biotype"))

gtf$gene_name <- sub('gene_name "([^"]+)"', "\\1", gtf$gene_name)
gtf <- gtf[seqname != "MT"]
gtf$seqname <- paste0("chr", gtf$seqname)
gtf <- gtf[seqname %in% hg38_chrom_sizes$chr]

## Filter for genes used in the 10x Genomics ENSEMBLv93 reference for compatibility
## with scRNA-seq analysis

## TODO: Replace with download from 10x website
tenx_feat <- fread(system.file("reference/GRCh38_10x_gene_metadata.csv.gz",
                               package = "H5weaver"))
tenx_ensembl <- tenx_feat$id

keep_gtf <- gtf[gene_id %in% tenx_ensembl]

fwrite(keep_gtf,
       "hg38_ensemble93_tenx_genes.tsv.gz")

# Gene Bodies
gene_gr <- GRanges(seqnames = keep_gtf$seqname,
                   ranges = IRanges(start = keep_gtf$start,
                           end = keep_gtf$end),
                   strand = keep_gtf$strand,
                   gene_id = keep_gtf$gene_id,
                   gene_name = keep_gtf$gene_name)
gene_gr <- sort(gene_gr, ignore.strand = TRUE)

saveRDS(gene_gr,
        "hg38_gene_bodies_gr.rds")

gene_bed <- as.data.frame(gene_gr)[,c(1:3, 6)]
fwrite(gene_bed,
       "hg38_gene_bodies.bed.gz",
       sep = "\t",
       quote = FALSE,
       col.names = FALSE,
       row.names = FALSE)

# TSS Regions
tss_2kb_gr <- resize(gene_gr,
                     width = 2e3,
                     fix = "start")
tss_2kb_gr <- resize(tss_2kb_gr,
                     width = 4e3,
                     fix = "end")

tss_2kb_gr <- GenomicRanges::sort(tss_2kb_gr, ignore.strand = TRUE)

# Rename S4Vectors DFrame to DataFrame for backwards compatibility
class(tss_2kb_gr@elementMetadata) <- "DataFrame"

saveRDS(tss_2kb_gr,
        "hg38_tss_gr.rds")

tss_2kb_bed <- as.data.frame(tss_2kb_gr)[,c(1:3, 6)]
fwrite(tss_2kb_bed,
       "hg38_tss.bed.gz",
       sep = "\t",
       quote = FALSE,
       col.names = FALSE,
       row.names = FALSE)

### Seurat RNA-seq reference
tenx_genes <- fread("hg38_ensemble93_tenx_genes.tsv.gz")
tenx_genes$symbols <- make.unique(tenx_genes$gene_name)

if(!file.exists("pbmc_10k_v3.rds")) {
  download.file("https://www.dropbox.com/s/3f3p5nxrn5b3y4y/pbmc_10k_v3.rds?dl=1",
                "pbmc_10k_v3.rds")
}

pbmc_so <- readRDS("pbmc_10k_v3.rds")
pbmc_mat <- pbmc_so@assays$RNA@counts

pbmc_genes <- rownames(pbmc_mat)
keep_genes <- pbmc_genes %in% tenx_genes$symbols
pbmc_mat <- pbmc_mat[keep_genes,]

pbmc_genes <- rownames(pbmc_mat)
pbmc_ensembl <- tenx_genes$gene_id[match(pbmc_genes, tenx_genes$symbols)]

pbmc_meta <- pbmc_so@meta.data
pbmc_meta$original_barcodes <- rownames(pbmc_meta)

pbmc_meta <- pbmc_meta[pbmc_meta$celltype != "Platelets",]
pbmc_mat <- pbmc_mat[,pbmc_meta$original_barcodes]

pbmc_meta <- lapply(pbmc_meta,
                    function(x) {
                      if(class(x) == "factor") {
                        as.character(x)
                      } else {
                        x
                      }
                    })

name_conversion <- c(B.Naive = "pre-B cell",
                     B.Activated = "B cell progenitor",
                     T.CD4.Naive = "CD4 Naive",
                     T.CD8.Naive = "CD8 Naive",
                     T.CD4.Memory = "CD4 Memory",
                     T.CD8.Effector = "CD8 effector",
                     T.DoubleNegative = "Double negative T cell",
                     NK = "NK cell",
                     DC.Plasmacytoid = "pDC",
                     DC.Myeloid = "Dendritic cell",
                     Mono.CD14 = "CD14+ Monocytes",
                     Mono.CD16 = "CD16+ Monocytes")

pbmc_meta$celltype <- names(name_conversion[match(pbmc_meta$celltype, name_conversion)])

pbmc_list <- list(matrix_dgCMatrix = pbmc_mat,
                  matrix = list(observations = pbmc_meta,
                                features = list(name = pbmc_genes))
                  )

pbmc_list <- h5_list_convert_from_dgCMatrix(pbmc_list, "matrix")

write_h5_list(pbmc_list, "seurat_rna_pbmc_ref.h5")

file.remove("pbmc_10k_v3.rds")

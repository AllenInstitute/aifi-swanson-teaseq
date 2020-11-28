library(data.table)
library(dplyr)
library(GenomicRanges)

system("curl https://resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0/hg38.archetype_motifs.v1.0.bed.gz > hg38.archetype_motifs.v1.0.bed.gz")
system("zcat hg38.archetype_motifs.v1.0.bed.gz | grep CTCF > hg38.CTCF_motifs.v1.0.bed.gz")

ctcf_motifs <- fread("hg38.CTCF_motifs.v1.0.bed.gz")
names(ctcf_motifs) <- c("seqnames","start","end",
                        "name","score","strand",
                        "best_model","num_models")


altius_gr <- readRDS(system.file("../reference/hg38_altius_gr.rds", package = "ATAComb"))

ctcf_gr <- GRanges(seqnames = ctcf_motifs$seqnames,
                   IRanges(start = ctcf_motifs$start,
                           end = ctcf_motifs$end))

ctcf_ol <- findOverlaps(ctcf_gr,
                        altius_gr)

keep_ctcf <- unique(S4Vectors::queryHits(ctcf_ol))

ctcf_motifs <- ctcf_motifs[keep_ctcf,]

ctcf_motifs <- ctcf_motifs %>%
  arrange(desc(score)) %>%
  head(100000) %>%
  arrange(seqnames, start, end)

fwrite(ctcf_motifs,
       "hg38.CTCF_motifs.v1.0_altius.overlap_top100k.bed.gz")

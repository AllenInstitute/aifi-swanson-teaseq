library(ArchR)
library(purrr)
library(dplyr)
library(ggrastr)
options(stringsAsFactors = FALSE)
set.seed(3030)

addArchRThreads(8)
addArchRGenome("hg38")
macs2 <- "/home/lucas.graybuck/miniconda3/bin/macs2"

adt_types <- read.csv("../../X044_adt_type_labels.csv")

drive_mount <- "/mnt/x044-atac-cite/ICICLE_V6/pipeline_analysis_2020-08-17/"

well_ids <- c("X044-AP0C1W1","X044-AP0C1W2","X044-AP0C1W3")
well_labels <- c("X044_icicle-3","X044_icicle-4","X044_icicle-5")

batch_ids <- sub("_.+","",well_labels)

qc_dirs <- paste(drive_mount,
                 well_ids, 
                 "atac_qc", sep = "/")

frag_files <- map_chr(qc_dirs,
                      function(qc_dir) {
                        list.files(qc_dir, 
                        pattern = "filtered_fragments.tsv.gz$",
                        full.names = TRUE)
                      })

newArrowFile <- createArrowFiles(inputFiles = frag_files,
                                 sampleName = well_labels,
                                 filterTSS = 0,
                                 filterFrags = 250,
                                 addTileMat = TRUE,
                                 TileMatParams = list(tileSize = 20000,
                                                      binarize = TRUE),
                                 addGeneScoreMat = TRUE)

doubScores <- addDoubletScores(input = newArrowFile,
                               k = 20,
                               knnMethod = "UMAP")

proj <- ArchRProject(ArrowFiles = newArrowFile,
                     outputDirectory = paste0(batch_ids[1],"_ArchRProject/"),
                     copyArrows = TRUE)

# proj <- filterDoublets(ArchRProj = proj,
#                        cutEnrich = 1)

proj_meta <- as.data.frame(getCellColData(proj))
proj_meta$archr_names <- rownames(proj_meta)
proj_meta$barcodes <- sub(".+#","",rownames(proj_meta))
proj_meta <- left_join(proj_meta, adt_types)

proj <- subsetArchRProject(proj,
                           cells = proj_meta$archr_names[!is.na(proj_meta$cluster_label)],
                           outputDirectory = paste0(batch_ids[1],"_filtered_ArchRProject/"))


proj <- addIterativeLSI(ArchRProj = proj,
                        useMatrix = "TileMatrix",
                        name = "IterativeLSI",
                        force = TRUE)

proj <- addClusters(input = proj,
                    reducedDims = "IterativeLSI",
                    knnAssign = 7,
                    force = TRUE)

proj <- addUMAP(ArchRProj = proj, 
                reducedDims = "IterativeLSI",
                nNeighbors = 21,
                force = TRUE)

umap_clusters_plot <- plotEmbedding(ArchRProj = proj, 
                                    colorBy = "cellColData", 
                                    name = "Clusters", 
                                    embedding = "UMAP")

ggsave("umap_clusters_plot.png",
       umap_clusters_plot,
       width = 8, height = 8)

archr_meta <- getCellColData(proj)
archr_meta$archr_name <- rownames(archr_meta)
archr_meta$barcodes <- sub(".+#","",archr_meta$archr_name)
archr_meta <- archr_meta %>%
  as.data.frame() %>%
  left_join(adt_types)

proj <- addCellColData(proj,
                       data = archr_meta$cluster_label,
                       name = "cluster_label",
                       cells = archr_meta$archr_name)

plotEmbedding(proj,
              name = "cluster_label")

proj <- addGroupCoverages(proj,
                          groupBy = "cluster_label")

proj <- addReproduciblePeakSet(proj,
                               groupBy = "cluster_label",
                               pathToMacs2 = macs2)

proj <- addPeakMatrix(proj)

proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")

markersPeaks <- getMarkerFeatures(proj,
                                  useMatrix = "PeakMatrix",
                                  groupBy = "cluster_label",
                                  bias = c("TSSEnrichment","log10(nFrags)"),
                                  testMethod = "wilcoxon")

heatmapPeaks <- plotMarkerHeatmap(seMarker = markersPeaks,
                              cutOff = "FDR <= 0.1",
                              transpose = TRUE)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1"
)

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, cutOff = 1, pMax = 10, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")


markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
markerList <- map(as.list(markerList),
                  as.data.frame)

saveArchRProject(proj,
                 paste0(batch_ids[1],"_filtered_ArchRProject/"))

## Add cell type labels from label transfer
reference_source <- system.file("reference/seurat_rna_pbmc_ref.h5",
                                package = "ATAComb")
file.copy(reference_source, "seurat_pbmc_ref.h5")
rna_so <- H5weaver::read_h5_seurat("seurat_pbmc_ref.h5")

proj <- addGeneIntegrationMatrix(ArchRProj = proj,
                                 useMatrix = "GeneScoreMatrix",
                                 matrixName = "GeneIntegrationMatrix",
                                 reducedDims = "IterativeLSI",
                                 seRNA = rna_so,
                                 groupRNA = "celltype",
                                 nameCell = "un_pred_type",
                                 nameGroup = "un_pred_group",
                                 nameScore = "un_pred_score",
                                 transferParams = list(dims = 1:10,
                                                       k.weight = 20),
                                 nGenes = 4000,
                                 plotUMAP = FALSE,
                                 force = TRUE
)

plotEmbedding(proj,
              name = "un_pred_group")

saveArchRProject(proj,
                 paste0(batch_ids[1],"_filtered_ArchRProject/"))


## Generate a river plot to compare the ADT and ArchR labels
adt_anno <- read.table("../../../common/x044_adt_cell_type_anno.tsv", 
                       sep = "\t", header = TRUE, comment.char = "|")
pbmc_anno <- read.table("../../../common/seurat_pbmc_type_anno.tsv", 
                        sep = "\t", header = TRUE, comment.char = "|")

pbmc_anno <- pbmc_anno %>%
  filter(seurat_pbmc_type %in% unique(proj$un_pred_group)) %>%
  mutate(seurat_pbmc_type_id = 1:n()) %>%
  rename(seurat_pbmc_type_label = seurat_pbmc_type)

proj_anno <- data.frame(sample_name = proj$cellNames,
                        peaks_frac = proj$FRIP,
                        label_score = proj$un_pred_score,
                        seurat_pbmc_type_label = proj$un_pred_group,
                        cluster_label = proj$cluster_label) %>%
  left_join(pbmc_anno) %>%
  left_join(adt_anno)

type_river <- scrattch.vis::build_river_plot(proj_anno,
                                             grouping = c("seurat_pbmc_type","cluster"),
                                             fill_group = "cluster",
                                             min_link_size = 0.05,
                                             pad = 0.2)

type_river

ggsave("Figure_3_type_river.pdf",
       type_river,
       width = 2, height = 4.5,
       useDingbats = FALSE
       )

proj_anno$barcodes <- sub(".+#","",proj_anno$sample_name)

write.csv(proj_anno,
          "X044_archr_proj_anno.csv")

proj_anno <- read.csv("X044_archr_proj_anno.csv")

atac_umap_coords <- read.csv("../../X044_analysis_results.csv")[,c("barcodes","umap_1","umap_2")]
adt_umap_coords <- read.csv("../../X044_adt_umap_results.csv")[,c("barcodes","adt_umap_1","adt_umap_2","cluster_label")]

proj_anno <- proj_anno %>% 
  left_join(atac_umap_coords) %>%
  left_join(adt_umap_coords)

# Output for supplementary data
adt_counts <- fread("../../X044_adt_count_df_full_depth.csv")
fwrite(adt_counts,
       "Figure-3_Supp-Data-2_ADT_UMI_Counts.csv.gz")

umap_type_df <- proj_anno %>%
  left_join(filter(adt_counts, target == "total")) %>%
  dplyr::rename(atac_type_label_score = label_score,
                atac_pbmc_type_label = seurat_pbmc_type_label,
                atac_pbmc_type_color = seurat_pbmc_type_color,
                adt_manual_type_label = cluster_label,
                adt_manual_type_color = cluster_color,
                atac_umap_1 = umap_1,
                atac_umap_2 = umap_2,
                adt_umis = count) %>%
  mutate(well_id = sub("#.+","",sample_name)) %>%
  select(well_id, barcodes,
         peaks_frac, atac_pbmc_type_label, atac_type_label_score, atac_pbmc_type_color,
         atac_umap_1, atac_umap_2,
         adt_umis, adt_manual_type_label, adt_manual_type_color,
         adt_umap_1, adt_umap_2)
fwrite(umap_type_df, "Figure-3_Supp-Data-1_type_labels_UMAP_coords.csv.gz")
  
## Plot Type UMAPs
atac_umap_coords <- atac_umap_coords %>%
  left_join(proj_anno)
adt_umap_coords <- adt_umap_coords %>%
  left_join(adt_anno)

atac_type_umap <- ggplot() +
  geom_point_rast(data = atac_umap_coords %>%
                    filter(!is.na(seurat_pbmc_type_color)),
                  aes(x = umap_1,
                      y = umap_2,
                      color = seurat_pbmc_type_color),
                  size = 0.1,
                  alpha = 0.2,
                  raster.dpi = 600) + 
  theme_bw(base_size = 7) +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.ticks = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_color_identity()

ggsave("Figure_3_atac_type_umap.pdf",
       atac_type_umap,
       width = 1.875, height = 1.875,
       useDingbats = FALSE)

adt_type_umap <- ggplot() +
  geom_point_rast(data = adt_umap_coords,
                  aes(x = adt_umap_1,
                      y = adt_umap_2,
                      color = cluster_color),
                  size = 0.1,
                  alpha = 0.2,
                  raster.dpi = 600) + 
  # geom_label(data = adt_umap_coords %>%
  #              group_by(cluster_label) %>%
  #              summarize(x = median(adt_umap_1),
  #                        y = median(adt_umap_2)),
  #            aes(x = x,
  #                y = y,
  #                label = cluster_label),
  #            alpha = 0.5) +
  theme_bw(base_size = 7) +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.ticks = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_color_identity()

ggsave("Figure_3_adt_type_umap.pdf",
       adt_type_umap,
       width = 1.875, height = 1.875,
       useDingbats = FALSE)

peak_files <- list.files("X044_filtered_ArchRProject/PeakCalls/",
                         pattern = ".rds",
                         full.names = TRUE)

peak_results <- map(peak_files,
                    readRDS)

unique_peak_results <- map(1:length(peak_results),
                           function(i) {
                             fg_peaks <- peak_results[[i]]
                             bg_peaks <- peak_results[-i]
                             bg_peaks <- do.call(c,bg_peaks)
                             bg_peaks <- GenomicRanges::reduce(bg_peaks)
                             
                             ol <- findOverlaps(fg_peaks, bg_peaks)
                             
                             qh <- unique(queryHits(ol))
                             
                             if(length(qh) < length(fg_peaks)) {
                               fg_peaks[-unique(queryHits(ol))]
                             } else {
                               NULL
                             }
                           })

n_unique_df <- data.frame(cell_type = sub(".+/","",sub("-.+","",peak_files)),
                          n_peaks = map_int(unique_peak_results, length))

unique_peaks_gr <- do.call(c, unique_peak_results)
unique_peaks_gr <- GenomicRanges::reduce(unique_peaks_gr)
unique_peaks_gr <- GenomicRanges::sort(unique_peaks_gr)

saveRDS(unique_peaks_gr,
        "X044_unique_peaks_gr.rds")


### Get all fragments for each cell type

all_proj_fragment_list <- map(getArrowFiles(proj),
                              getFragmentsFromArrow)

cell_types <- unique(archr_meta$cluster_label)
type_fragment_list <- map(cell_types,
                          function(cell_type) {
                            cell_names <- archr_meta$archr_name[archr_meta$cluster_label == cell_type]
                            
                            res <- map(all_proj_fragment_list,
                                       function(fragments) {
                                         fragments[fragments$RG %in% cell_names]
                                       })
                            
                            gr_res <- res[[1]]
                            for(i in 2:length(res)) {
                              gr_res <- c(gr_res, res[[i]])
                            }
                            GenomicRanges::sort(gr_res)
                          })

type_bed_list <- map(type_fragment_list,
                     gr_to_bed)

dir.create("X044_type_bed")
walk2(type_bed_list,
      cell_types,
     function(bed, cell_type) {
       write.table(bed,
                   paste0("X044_type_bed/",cell_type,".bed"),
                   row.names = FALSE,
                   col.names = FALSE,
                   quote = FALSE,
                   sep = "\t")
     })

dir.create("X044_type_peaks")
walk(cell_types,
     function(cell_type) {
       macs2_command <- paste("/shared/apps/miniconda3/bin/macs2 callpeak",
                              "-f BED",
                              "-g hs",
                              "-t",paste0("X044_type_bed/",cell_type,".bed"),
                              "--outdir X044_type_peaks",
                              "-n",cell_type,
                              "--nomodel"
                              )
       system(macs2_command)
     })

# Make bedgraphs and bigwigs for display
# Generate bedgraphs
dir.create("X044_type_bedgraph")

gr_to_bedgraph <- function(gr) {
  gr <- sort(gr)
  gr_cov <- as.list(coverage(gr))
  
  gr_ends <- lapply(gr_cov, function(x) { cumsum(runLength(x)) } )
  gr_starts <- lapply(gr_ends, function(x) { c(1,(x[1:(length(x)-1)] + 1)) } )
  gr_scores <- lapply(gr_cov, function(x) { runValue(x) } )
  gr_lens <- sapply(gr_ends, length)
  
  bg_df <- data.frame(chr = rep(names(gr_lens), gr_lens),
                      start = as.integer(unlist(gr_starts)),
                      end = as.integer(unlist(gr_ends)),
                      score = as.integer(unlist(gr_scores)))
  
  bg_df[bg_df$score > 0,]
  
}

walk2(type_fragment_list, cell_types,
     function(type_fragments, cell_type) {
       cluster_bg <- gr_to_bedgraph(type_fragments)
       
       bg_file <- paste0(cell_type, ".bedgraph")
       
       write.table(cluster_bg,
                   file = file.path("X044_type_bedgraph",bg_file),
                   quote = FALSE,
                   sep = "\t",
                   row.names = FALSE,
                   col.names = FALSE)
     })

# Convert bedgraph to bigwig
dir.create("X044_type_bigwig")
bg_to_bw <- "/shared/apps/ucsctools/bedGraphToBigWig"
file.copy(system.file("reference/hg38.chrom.sizes", package = "ATAComb"),
          "hg38.chrom.sizes")

chrom_sizes <- "hg38.chrom.sizes"

walk(cell_types,
     function(cell_type) {
       bg_file <- paste0(cell_type, ".bedgraph")
       bg_file <- file.path("X044_type_bedgraph",bg_file)
       bw_file <- paste0(cell_type, ".bw")
       bw_file <- file.path("X044_type_bigwig",bw_file)
       
       sort_command <- paste("cat", bg_file, "| sort -k1,1 -k2,2n > temp.bedgraph; mv temp.bedgraph", bg_file)
       system(sort_command)
       conv_command <- paste(bg_to_bw, bg_file, chrom_sizes, bw_file)
       system(conv_command)
     })




## Assemble a pseudobulk matrix for projection onto singlecells
archr_se <- getMatrixFromProject(proj,
                                  useMatrix = "TileMatrix",
                                  binarize = TRUE)

archr_mat <- assays(archr_se)@listData$TileMatrix

cluster_labels <- as.data.frame(getCellColData(proj, select = "cluster_label"))
unique_labels <- unique(cluster_labels$cluster_label)

pseudobulk_mat_list <- map(unique_labels,
                           function(cluster_label) {
                             print(cluster_label)
                             cluster_barcodes <- rownames(cluster_labels)[cluster_labels$cluster_label == cluster_label]
                             cluster_mat <- archr_mat[,cluster_barcodes, drop = FALSE]
                             cluster_sums <- rowSums(cluster_mat)
                             mat <- Matrix::sparseMatrix(x = cluster_sums[cluster_sums > 0],
                                                         i = which(cluster_sums > 0),
                                                         j = rep(1, length(cluster_sums[cluster_sums > 0])),
                                                         dims = c(length(cluster_sums), 1))
                             colnames(mat) <- cluster_label
                             as(mat, "dgTMatrix")
                           })

pseudobulk_mat <- do.call(cbind, pseudobulk_mat_list)

pseudobulk_se <- SummarizedExperiment(assays = list(counts = pseudobulk_mat))

pseudobulk_gr <- GRanges(rep(archr_se@elementMetadata@listData$seqnames@values,
                             archr_se@elementMetadata@listData$seqnames@lengths),
                         IRanges(start = archr_se@elementMetadata@listData$start + 1,
                                 end = archr_se@elementMetadata@listData$start + 2e4))
rowRanges(pseudobulk_se) <- pseudobulk_gr

saveRDS(pseudobulk_se,
        "X044_pseudobulk_sce.rds")


#Pairwise tests
all_types <- unique(proj$cluster_label)
pw_types <- combn(all_types, 2)

pw_results <- map2(pw_types[1,],
                   pw_types[2,],
                   function(fg, bg) {
                     list(
                       fg = getMarkerFeatures(
                         ArchRProj = proj,
                         useMatrix = "PeakMatrix",
                         groupBy = "cluster_label",
                         testMethod = "wilcoxon",
                         bias = c("TSSEnrichment","log10(nFrags)"),
                         useGroups = fg,
                         bgdGroups = bg
                       ),
                       bg = getMarkerFeatures(
                         ArchRProj = proj,
                         useMatrix = "PeakMatrix",
                         groupBy = "cluster_label",
                         testMethod = "wilcoxon",
                         bias = c("TSSEnrichment","log10(nFrags)"),
                         useGroups = bg,
                         bgdGroups = fg
                       )
                     )
                   })

marker_results_unfiltered <- map(pw_results,
                      function(pw_result) {
                        list(
                          fg = getMarkers(pw_result[["fg"]]),
                          bg = getMarkers(pw_result[["bg"]])
                        )
                      })

saveRDS(marker_results_unfiltered,
        "all_types_pairwise_unfiltered.rds")

cutoff <- "FDR <= 0.1"

marker_results <- map(pw_results,
                      function(pw_result) {
                        list(
                          fg = getMarkers(pw_result[["fg"]],
                                          cutOff = cutoff),
                          bg = getMarkers(pw_result[["bg"]],
                                          cutOff = cutoff)
                        )
                      })

saveRDS(marker_results,
        "all_types_pairwise.rds")

pw_summary <- map2_dfr(1:ncol(pw_types),
                       marker_results,
                       function(pw_col, pw_res) {
                         
                         fg_df <- as.data.frame(pw_res[["fg"]]@listData[[1]])
                         
                         fg_res <- data.frame(fg = pw_types[1, pw_col],
                                              bg = pw_types[2, pw_col],
                                              n_DAP = nrow(fg_df),
                                              mean_log2fc = mean(fg_df$Log2FC))
                         
                         bg_df <- as.data.frame(pw_res[["bg"]]@listData[[1]])
                         bg_res <- data.frame(fg = pw_types[2, pw_col],
                                              bg = pw_types[1, pw_col],
                                              n_DAP = nrow(bg_df),
                                              mean_log2fc = mean(bg_df$Log2FC))
                         rbind(fg_res, bg_res)
                       })

ggplot(pw_summary) +
  geom_point(aes(x = log10(n_DAP),
                 y = mean_log2fc))

pw_plot_df <- pw_summary

type_pos <- data.frame(type=c("B.Naive",
                              "B.Active",
                              "T.CD4.Naive",
                              "T.CD4.Memory",
                              "T.CD4.CD45RO.CD45RA",
                              "T.CD4.Exhausted",
                              "T.CD8.Naive",
                              "T.CD8.Effector",
                              "T.CD8.Exhausted",
                              "T.DoublePositive",
                              "T.DoubleNegative",
                              "NK",
                              "DC.Myeloid",
                              "Mono.CD14",
                              "Mono.CD16",
                              "DC.Plasmacytoid"))
type_pos$xpos <- 1:nrow(type_pos)
type_pos$ypos <- nrow(type_pos):1

pw_plot_df$ypos <- type_pos$ypos[match(pw_summary$fg, type_pos$type)]
pw_plot_df$xpos <- type_pos$xpos[match(pw_summary$bg, type_pos$type)]

ggplot(pw_plot_df) +
  geom_tile(aes(x = xpos,
                y = ypos,
                fill = log10(n_DAP)),
            color = "gray40",
            size = 0.1) +
  scale_x_continuous("Background Type",
                     breaks = type_pos$xpos,
                     labels = type_pos$type ,
                     expand = c(0,0)) +
  scale_y_continuous("Foreground Type",
                     breaks = type_pos$ypos,
                     labels = type_pos$type, expand = c(0,0)) +
  scale_fill_viridis_c() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.3),
        axis.text = element_text(color = "#000000"),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "#000000")) +
  ggtitle("Perm: No Neutrophils")

# Generate bed and bam files to try diffbind
dir.create("X044_cluster_label_bed")

gr_to_bed <- function(gr) {
  gr <- sort(gr)
  data.frame(chr = seqnames(gr),
             start = as.integer(start(gr)),
             end = as.integer(end(gr)),
             name = "",
             score = 1,
             strand = "+")
}

coverage_files <- list.files("X044_filtered_ArchRProject/GroupCoverages/cluster_label/")

walk(coverage_files,
     function(coverage_file) {
       h5_file <- file.path("X044_filtered_ArchRProject/GroupCoverages/cluster_label",coverage_file)
       coverage_cells <- as.vector(rhdf5::h5read(h5_file, "/Coverage/Info/CellNames"))
       coverage_samples <- unique(sub("#.+","",coverage_cells))
       coverage_arrows <- paste0(coverage_samples, ".arrow")

       cluster_label <- sub("\\._.+","",coverage_file)
       cluster_meta <- as.data.frame(archr_meta[archr_meta$cluster_label == cluster_label,])
       
       cluster_grs <- map(coverage_arrows,
                          function(arrowFile) {
                            arrow_sample <- sub(".arrow","",arrowFile)
                            arrow_cells <- cluster_meta$archr_name[cluster_meta$Sample == arrow_sample]
                            getFragmentsFromArrow(
                              ArrowFile = arrowFile,
                              cellNames = arrow_cells
                            )
                            
                          })
       cluster_gr <- do.call(c, cluster_grs)
       
       cluster_bed <- gr_to_bed(cluster_gr)
       
       bed_file <- sub(".insertions.coverage.h5",".bed",coverage_file)
       
       write.table(cluster_bed,
                   file = file.path("X044_cluster_label_bed",bed_file),
                   quote = FALSE,
                   sep = "\t",
                   row.names = FALSE,
                   col.names = FALSE)
       
     })

dir.create("X044_cluster_label_bam")
file.copy(system.file("reference/hg38.chrom.sizes", package = "ATAComb"),
          "hg38.chrom.sizes")

bed_files <- list.files("X044_cluster_label_bed")

walk(bed_files,
     function(bed_file) {
       bam_file <- sub(".bed", ".bam", bed_file)
       
       command <- paste("/shared/apps/bedtools2/bin/bedtools bedtobam",
                        "-i", file.path("X044_cluster_label_bed", bed_file),
                        "-g hg38.chrom.sizes",
                        ">", file.path("X044_cluster_label_bam", bam_file))
       
       system(command)
       
     })


dir.create("X044_cluster_label_peak_bed")

walk2(peak_results,peak_files,
     function(peak_result,peak_file) {
       peak_bed <- gr_to_bed(peak_result)
       bed_file <- sub(".+/","",sub("-.+","_peaks.bed",peak_file))
       write.table(peak_bed,
                   file = file.path("X044_cluster_label_peak_bed",bed_file),
                   quote = FALSE,
                   sep = "\t",
                   row.names = FALSE,
                   col.names = FALSE)
     })

# Generate bedgraphs
dir.create("X044_cluster_label_bedgraph")

gr_to_bedgraph <- function(gr) {
  gr <- sort(gr)
  gr_cov <- as.list(coverage(gr))
  
  gr_ends <- lapply(gr_cov, function(x) { cumsum(runLength(x)) } )
  gr_starts <- lapply(gr_ends, function(x) { c(1,(x[1:(length(x)-1)] + 1)) } )
  gr_scores <- lapply(gr_cov, function(x) { runValue(x) } )
  gr_lens <- sapply(gr_ends, length)
  
  bg_df <- data.frame(chr = rep(names(gr_lens), gr_lens),
                      start = as.integer(unlist(gr_starts)),
                      end = as.integer(unlist(gr_ends)),
                      score = as.integer(unlist(gr_scores)))
  
  bg_df[bg_df$score > 0,]
  
}

walk(unique_labels,
     function(cluster_label) {
       cluster_meta <- archr_meta[archr_meta$cluster_label == cluster_label,]
       cluster_grs <- map(newArrowFile,
                         function(arrowFile) {
                           arrow_sample <- sub(".arrow","",arrowFile)
                           arrow_cells <- cluster_meta$archr_name[cluster_meta$Sample == arrow_sample]
                           getFragmentsFromArrow(
                             ArrowFile = arrowFile,
                             cellNames = arrow_cells
                           )
                         })
       cluster_gr <- do.call(c, cluster_grs)
       
       cluster_bg <- gr_to_bedgraph(cluster_gr)
       
       bg_file <- paste0(cluster_label, ".bedgraph")
       write.table(cluster_bg,
                   file = file.path("X044_cluster_label_bedgraph",bg_file),
                   quote = FALSE,
                   sep = "\t",
                   row.names = FALSE,
                   col.names = FALSE)
     })

# Convert bedgraph to bigwig
dir.create("X044_cluster_label_bigwig")
bg_to_bw <- "/shared/apps/ucsctools/bedGraphToBigWig"
file.copy(system.file("reference/hg38.chrom.sizes", package = "ATAComb"),
          "hg38.chrom.sizes")

chrom_sizes <- "hg38.chrom.sizes"

walk(unique_labels,
     function(cluster_label) {
       bg_file <- paste0(cluster_label, ".bedgraph")
       bg_file <- file.path("X044_cluster_label_bedgraph",bg_file)
       bw_file <- paste0(cluster_label, ".bw")
       bw_file <- file.path("X044_cluster_label_bigwig",bw_file)
       
       sort_command <- paste("cat", bg_file, "| sort -k1,1 -k2,2n > temp.bedgraph; mv temp.bedgraph", bg_file)
       system(sort_command)
       conv_command <- paste(bg_to_bw, bg_file, chrom_sizes, bw_file)
       system(conv_command)
     })

# all_marker_peaks <- map_dfr(marker_results,
#                             function(marker_result) {
#                               
#                               fg_df <- as.data.frame(marker_result[[1]][[1]])
#                               bg_df <- as.data.frame(marker_result[[2]][[1]])
#                               
#                               if(nrow(fg_df) & nrow(bg_df) > 0) {
#                                 joint_df <- rbind(fg_df[,c("seqnames","idx","start","end")],
#                                                   bg_df[,c("seqnames","idx","start","end")])
#                               } else if(nrow(fg_df) > 0) {
#                                 fg_df
#                               } else {
#                                 bg_df
#                               }
#                             })
# 
# all_marker_peaks <- unique(all_marker_peaks)
# all_marker_gr <- GRanges(seqnames = all_marker_peaks$seqnames,
#                          IRanges(start = all_marker_peaks$start,
#                                  end = all_marker_peaks$end))
#   
# proj2 <- proj
# 
# 
# proj2 <- addPeakSet(proj2,
#                     peakSet = all_marker_gr,
#                     force = TRUE)
# 
# proj2 <- addPeakMatrix(proj2,
#                        force = TRUE)
# 
# proj2 <- addIterativeLSI(proj,
#                          useMatrix = "PeakMatrix",
#                          iterations = 1,
#                          force = TRUE)
# 
# proj2 <- addClusters(input = proj2,
#                      reducedDims = "IterativeLSI",
#                      knnAssign = 7,
#                      force = TRUE)
# 
# proj2 <- addUMAP(proj2, 
#                  reducedDims = "IterativeLSI",
#                  nNeighbors = 40,
#                  force = TRUE)
# 
# umap_clusters_plot2 <- plotEmbedding(ArchRProj = proj2, 
#                                      colorBy = "cellColData", 
#                                      name = "Clusters", 
#                                      embedding = "UMAP")
# umap_clusters_plot2
# 
# umap_clusters_plot3 <- plotEmbedding(ArchRProj = proj2, 
#                                      colorBy = "cellColData", 
#                                      name = "cluster_label", 
#                                      embedding = "UMAP")
# umap_clusters_plot3
# 
# umap_clusters_plot4 <- plotEmbedding(ArchRProj = proj2, 
#                                      colorBy = "cellColData", 
#                                      name = "nFrags", 
#                                      embedding = "UMAP")
# umap_clusters_plot4
# 
# outDir <- "X037_ArchRProject/"
# proj <- addGeneIntegrationMatrix(ArchRProj = proj,
#                                  addToArrow = FALSE,
#                                  useMatrix = "GeneScoreMatrix",
#                                  matrixName = "GeneIntegrationMatrix",
#                                  reducedDims = "IterativeLSI",
#                                  seRNA = rna_so,
#                                  groupRNA = "celltype",
#                                  nameCell = "un_pred_type",
#                                  nameGroup = "un_pred_group",
#                                  nameScore = "un_pred_score",
#                                  force = TRUE
# )
# 
# un_conf_mat <- as.matrix(confusionMatrix(proj$Clusters, proj$un_pred_group))
# 
# png("un_conf_mat.png",width = 72*8, height = 72*8)
# heatmap(un_conf_mat)
# dev.off()
# 
# rna_types <- unique(rna_so@meta.data$celltype)
# 
# rna_patterns <- list(b_pdc = "^B|Plasmacytoid",
#                      mono_dc = "^Mono|Myeloid",
#                      tnk_cell = "^T|NK")
# 
# rna_clusters <- lapply(rna_patterns,
#                        function(x) {
#                          rna_types[grepl(x, rna_types)]
#                        })
# 
# rna_cells <- lapply(rna_clusters,
#                     function(x) {
#                       colnames(rna_so)[rna_so@meta.data$celltype %in% x]
#                     })
# 
# preClust <- colnames(un_conf_mat)[apply(un_conf_mat, 1 , which.max)]
# 
# atac_clusters <- lapply(rna_clusters,
#                         function(x) {
#                           rownames(un_conf_mat)[preClust %in% x]
#                         })
# 
# atac_cells <- lapply(atac_clusters,
#                      function(x) {
#                        proj$cellNames[proj$Clusters %in% x]
#                      })
# 
# groupList <- lapply(names(atac_cells),
#                     function(x) {
#                       SimpleList(ATAC = atac_cells[[x]],
#                                  RNA = rna_cells[[x]])
#                     })
# groupList <- SimpleList(groupList)
# 
# proj <- addGeneIntegrationMatrix(ArchRProj = proj,
#                                  useMatrix = "GeneScoreMatrix",
#                                  matrixName = "GeneIntegrationMatrix",
#                                  reducedDims = "IterativeLSI",
#                                  seRNA = rna_so,
#                                  groupList = groupList,
#                                  groupRNA = "celltype",
#                                  nameCell = "co_pred_type",
#                                  nameGroup = "co_pred_group",
#                                  nameScore = "co_pred_score",
#                                  force = TRUE
# )
# 
# co_conf_mat <- as.matrix(confusionMatrix(proj$Clusters, proj$co_pred_group))
# 
# png("co_conf_mat.png",width = 72*8, height = 72*8)
# heatmap(co_conf_mat)
# dev.off()
# 
# old_labels <- rownames(co_conf_mat)
# new_labels <- colnames(co_conf_mat)[apply(co_conf_mat, 1, which.max)]
# 
# label_map <- mapLabels(old_labels,
#                        oldLabels = old_labels,
#                        newLabels = new_labels)
# label_map
# 
# proj$seurat_pbmc_type <- mapLabels(proj$Clusters, 
#                                    newLabels = label_map, 
#                                    oldLabels = old_labels)
# 
# results <- data.frame(barcodes = sub(".+#","",rownames(proj@cellColData)),
#                       doublet_score = proj@cellColData$DoubletScore,
#                       doublet_enrichment = proj@cellColData$DoubletEnrichment,
#                       seurat_pbmc_type = proj$seurat_pbmc_type,
#                       seurat_pbmc_type_score = proj$co_pred_score)
# 
# out_table_path <- paste0(batch_ids[1], "_archr_results.csv.gz")
# fwrite(results,
#        out_table_path)
# 
# saveArchRProject(proj, "X037_ArchRProject")
# 
# umap_types_plot <- plotEmbedding(ArchRProj = proj, 
#                                     colorBy = "cellColData", 
#                                     name = "seurat_pbmc_type ", 
#                                     embedding = "UMAP")
# 
# ggsave("umap_types_plot.png",
#        umap_types_plot,
#        width = 8, height = 8)
# 
# 

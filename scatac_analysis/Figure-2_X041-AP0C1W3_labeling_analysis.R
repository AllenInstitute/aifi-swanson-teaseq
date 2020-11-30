library(ArchR)
library(purrr)
options(stringsAsFactors = FALSE)
set.seed(3030)

addArchRThreads(30)
addArchRGenome("hg38")
macs2 <- "/home/lucas.graybuck/miniconda3/bin/macs2"

proj <- loadArchRProject("X041-AP0C1W3_ArchRProject/")

proj <- addGroupCoverages(ArchRProj = proj,
                          groupBy = "seurat_pbmc_type")

proj <- addReproduciblePeakSet(ArchRProj = proj,
                               groupBy = "seurat_pbmc_type",
                               pathToMacs2 = macs2)

proj <- addPeakMatrix(proj)

proj <- addMotifAnnotations(ArchRProj = proj,
                            motifSet = "cisbp",
                            name = "Motif")

markersPeaks <- getMarkerFeatures(proj,
                                  useMatrix = "PeakMatrix",
                                  groupBy = "seurat_pbmc_type",
                                  bias = c("TSSEnrichment","log10(nFrags)"),
                                  testMethod = "wilcoxon")

enrichMotifs <- peakAnnoEnrichment(seMarker = markersPeaks,
                                   ArchRProj = proj,
                                   peakAnnotation = "Motif",
                                   cutOff = "FDR <= 0.01 & Log2FC >= log2(1.5)")

saveRDS(enrichMotifs,
        "X041-AP0C1W3_enrichMotifs.rds")

enrichMotifs <- readRDS("X041-AP0C1W3_enrichMotifs.rds")

motifs_df <- assay(enrichMotifs) %>%
  mutate(motif = rownames(enrichMotifs)) %>%
  reshape2::melt() %>%
  dplyr::rename(cell_type = variable,
                neg_log10_padj = value) %>%
  mutate(well_label = "Nuclei: No Neutrophils",
         well_id = "X041-AP0C1W3") %>%
  select(well_label, well_id, cell_type, motif, neg_log10_padj)

fwrite(motifs_df,
       "X041-AP0C1W1_motif_enrichment.csv.gz")

heatmapEM <- plotEnrichHeatmap(enrichMotifs, 
                               n = 10, 
                               transpose = FALSE)

ComplexHeatmap::draw(heatmapEM, 
                     heatmap_legend_side = "bot", 
                     annotation_legend_side = "bot")

proj <- saveArchRProject(proj, "X024-AP0C1W1_ArchRProject")

all_types <- unique(proj$seurat_pbmc_type)
pw_types <- combn(all_types, 2)

pw_results <- map2(pw_types[1,],
                   pw_types[2,],
                   function(fg, bg) {
                     list(
                       fg = getMarkerFeatures(
                         ArchRProj = proj,
                         useMatrix = "PeakMatrix",
                         groupBy = "seurat_pbmc_type",
                         testMethod = "wilcoxon",
                         bias = c("TSSEnrichment","log10(nFrags)"),
                         useGroups = fg,
                         bgdGroups = bg
                       ),
                       bg = getMarkerFeatures(
                         ArchRProj = proj,
                         useMatrix = "PeakMatrix",
                         groupBy = "seurat_pbmc_type",
                         testMethod = "wilcoxon",
                         bias = c("TSSEnrichment","log10(nFrags)"),
                         useGroups = bg,
                         bgdGroups = fg
                       )
                     )
                   })

cutoff <- "FDR <= 0.05 & Log2FC >= log2(1.5)"

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

marker_results <- readRDS("all_types_pairwise.rds")

pw_df <- map2_dfr(1:ncol(pw_types),
                  marker_results,
                  function(pw_col, pw_res) {
                    
                    fg_df <- as.data.frame(pw_res[["fg"]]@listData[[1]])
                    
                    if(nrow(fg_df) > 0) {
                      fg_df$fg_type <- pw_types[1, pw_col]
                      fg_df$bg_type <- pw_types[2, pw_col]
                    }
                    
                    bg_df <- as.data.frame(pw_res[["bg"]]@listData[[1]])
                    
                    if(nrow(bg_df) > 0) {
                      bg_df$fg_type <- pw_types[2, pw_col]
                      bg_df$bg_type <- pw_types[1, pw_col]
                    }
                    
                    res <- rbind(fg_df, bg_df)
                    if(nrow(res) > 0) {
                      res$well_label <- "Nuclei: No Neutrophils"
                      res$well_id <- "X041-AP0C1W3"
                      res <- res[,c("well_label","well_id",
                                    "fg_type","bg_type",
                                    "seqnames","start","end",
                                    "Log2FC","FDR","MeanDiff")]
                      names(res)[names(res) == "seqnames"] <- "chr"
                      res
                    }
                    
                  })

fwrite(pw_df,
       "X041-AP0C1W43_pw_peak_df.csv.gz")


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

type_pos <- data.frame(type=c("Mono.CD14",
                              "NK",
                              "T.CD4.Naive",
                              "B.Naive",
                              "Mono.CD16",
                              #"DC.Plasmacytoid",
                              "T.CD4.Memory"))
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
  ggtitle("Nuclei: 10x")

saveRDS(pw_plot_df,
        "X041-AP0C1W3_pw_plot_df.rds")
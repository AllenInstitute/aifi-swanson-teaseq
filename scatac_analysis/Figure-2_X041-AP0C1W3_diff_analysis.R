library(ArchR)
options(stringsAsFactors = FALSE)
set.seed(3030)

addArchRThreads(30)
addArchRGenome("hg38")

reference_source <- system.file("reference/seurat_rna_pbmc_ref.h5",
                                package = "ATAComb")
file.copy(reference_source, "seurat_pbmc_ref.h5")
seurat_reference <- "seurat_pbmc_ref.h5"
rna_so <- H5weaver::read_h5_seurat(seurat_reference)

drive_mount <- "/mnt/x041-atac-in"
downsample_level <- "200M"

well_id <- "X041-AP0C1W3"
well_label <- "X041_aifi-nuclei-facs-sort-neu-dd"

batch_id <- sub("_.+","",well_label)

qc_dir <- paste(drive_mount,
                paste0("downsampled_",
                       downsample_level),
                batch_id, well_id, 
                "atac_qc", sep = "/")

frag_file <- list.files(qc_dir, 
                        pattern = "filtered_fragments.tsv.gz$",
                        full.names = TRUE)

newArrowFile <- createArrowFiles(inputFiles = frag_file,
                                 sampleName = well_label,
                                 filterTSS = 4,
                                 filterFrags = 1000,
                                 addTileMat = TRUE,
                                 addGeneScoreMat = TRUE)

doubScores <- addDoubletScores(input = newArrowFile,
                               k = 20,
                               knnMethod = "UMAP")

proj <- ArchRProject(ArrowFiles = newArrowFile,
                     outputDirectory = paste0(well_id,"_ArchRProject/"),
                     copyArrows = TRUE)

proj <- filterDoublets(ArchRProj = proj,
                       cutEnrich = 1)

proj <- addIterativeLSI(ArchRProj = proj,
                        name = "IterativeLSI",
                        force = TRUE)
proj <- addImputeWeights(proj)

proj <- addUMAP(ArchRProj = proj, 
                reducedDims = "IterativeLSI",
                force = TRUE,
                saveModel = FALSE)

proj <- addClusters(input = proj,
                    reducedDims = "IterativeLSI",
                    resolution = 2,
                    force = TRUE)

umap_clusters_plot <- plotEmbedding(ArchRProj = proj, 
                    colorBy = "cellColData", 
                    name = "Clusters", 
                    embedding = "UMAP")

ggsave("umap_clusters_plot.png",
       umap_clusters_plot,
       width = 8, height = 8)

outDir <- "X0241-AP0C1W3_ArchRProject/"
proj <- addGeneIntegrationMatrix(ArchRProj = proj,
                                 addToArrow = FALSE,
                                 useMatrix = "GeneScoreMatrix",
                                 matrixName = "GeneIntegrationMatrix",
                                 reducedDims = "IterativeLSI",
                                 seRNA = rna_so,
                                 groupRNA = "celltype",
                                 nameCell = "un_pred_type",
                                 nameGroup = "un_pred_group",
                                 nameScore = "un_pred_score",
                                 force = TRUE
)

un_conf_mat <- as.matrix(confusionMatrix(proj$Clusters, proj$un_pred_group))

png("un_conf_mat.png",width = 72*8, height = 72*8)
heatmap(un_conf_mat)
dev.off()

rna_types <- unique(rna_so@meta.data$celltype)

rna_patterns <- list(b_pdc = "^B|Plasmacytoid",
                     mono_dc = "^Mono|Myeloid",
                     tnk_cell = "^T|NK")

rna_clusters <- lapply(rna_patterns,
                       function(x) {
                         rna_types[grepl(x, rna_types)]
                       })

rna_cells <- lapply(rna_clusters,
                    function(x) {
                      colnames(rna_so)[rna_so@meta.data$celltype %in% x]
                    })

preClust <- colnames(un_conf_mat)[apply(un_conf_mat, 1 , which.max)]

atac_clusters <- lapply(rna_clusters,
                        function(x) {
                          rownames(un_conf_mat)[preClust %in% x]
                        })

atac_cells <- lapply(atac_clusters,
                     function(x) {
                       proj$cellNames[proj$Clusters %in% x]
                     })

groupList <- lapply(names(atac_cells),
                    function(x) {
                      SimpleList(ATAC = atac_cells[[x]],
                                 RNA = rna_cells[[x]])
                    })
groupList <- SimpleList(groupList)

proj <- addGeneIntegrationMatrix(ArchRProj = proj,
                                 useMatrix = "GeneScoreMatrix",
                                 matrixName = "GeneIntegrationMatrix",
                                 reducedDims = "IterativeLSI",
                                 seRNA = rna_so,
                                 groupList = groupList,
                                 groupRNA = "celltype",
                                 nameCell = "co_pred_type",
                                 nameGroup = "co_pred_group",
                                 nameScore = "co_pred_score",
                                 force = TRUE
)

co_conf_mat <- as.matrix(confusionMatrix(proj$Clusters, proj$co_pred_group))

png("co_conf_mat.png",width = 72*8, height = 72*8)
heatmap(co_conf_mat)
dev.off()

old_labels <- rownames(co_conf_mat)
new_labels <- colnames(co_conf_mat)[apply(co_conf_mat, 1, which.max)]

label_map <- mapLabels(old_labels,
                       oldLabels = old_labels,
                       newLabels = new_labels)
label_map

proj$seurat_pbmc_type <- mapLabels(proj$Clusters, 
                                   newLabels = label_map, 
                                   oldLabels = old_labels)

results <- data.frame(barcodes = sub(".+#","",rownames(proj@cellColData)),
                      doublet_score = proj@cellColData$DoubletScore,
                      doublet_enrichment = proj@cellColData$DoubletEnrichment,
                      seurat_pbmc_type = proj$seurat_pbmc_type,
                      seurat_pbmc_type_score = proj$co_pred_score)

out_table_path <- paste0(well_label, "_archr_results.csv.gz")
fwrite(results,
       out_table_path)

saveArchRProject(proj, "X041-AP0C1W3_ArchRProject")

umap_types_plot <- plotEmbedding(ArchRProj = proj, 
                                    colorBy = "cellColData", 
                                    name = "seurat_pbmc_type ", 
                                    embedding = "UMAP")

ggsave("umap_types_plot.png",
       umap_types_plot,
       width = 8, height = 8)



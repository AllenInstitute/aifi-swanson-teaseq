library(ArchR)
library(data.table)

addArchRGenome("hg38")
addArchRThreads(30)

if(!dir.exists("ArchROutput")) {
  # Trying ArchR's defaults
  cellranger_barcodes <- fread("/mnt/x060-multiome/X061-Well1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", header = FALSE)
  keep_barcodes <- SimpleList(list("X061-AP0C1W1" = cellranger_barcodes[["V1"]]))
  
  arrowFile <- createArrowFiles("/mnt/x060-multiome/X061-Well1/outs/atac_fragments.tsv.gz",
                                sampleNames = "X061-AP0C1W1",
                                validBarcodes = keep_barcodes)
  
  proj <- ArchRProject(arrowFile)
} else {
  proj <- loadArchRProject("ArchROutput")
}

# scATAC-seq analysis

proj <- addIterativeLSI(proj)

proj <- addClusters(proj)

proj <- addUMAP(proj,
                saveModel = FALSE)

plotEmbedding(proj, name = "Clusters")

macs2 <- "/shared/apps/miniconda3/bin/macs2"

proj <- addGroupCoverages(proj,
                          groupBy = "Clusters")

proj <- addReproduciblePeakSet(proj,
                               groupBy = "Clusters",
                               pathToMacs2 = macs2)

proj <- addPeakMatrix(proj)

proj <- saveArchRProject(proj)

peak_res <- getMatrixFromProject(proj, "PeakMatrix")
peak_gr <- getPeakSet(proj)
peak_names <- paste0(seqnames(peak_gr), ":" ,start(peak_gr), "-" ,end(peak_gr))

peak_mat <- assays(peak_res)@listData$PeakMatrix
rownames(peak_mat) <- peak_names
colnames(peak_mat) <- sub(".+#", "", colnames(peak_mat))

saveRDS(peak_mat,
        "atac_peak_mat.rds")

proj <- addIterativeLSI(proj,
                        useMatrix = "PeakMatrix",
                        name = "PeakLSI")

proj <- addUMAP(proj,
                reducedDims = "PeakLSI",
                name = "PeakUMAP",
                saveModel = FALSE)

plotEmbedding(proj, "PeakUMAP",
              name = "Clusters")

archr_meta <- getCellColData(proj)

peak_umap <- getEmbedding(proj, "PeakUMAP", returnDF = TRUE) %>%
  as.data.frame() %>%
  mutate(archr_name = rownames(.)) %>%
  mutate(original_barcode = sub(".+#","",archr_name)) %>%
  mutate(archr_cluster = archr_meta$Clusters)
names(peak_umap)[1:2] <- c("peak_umap_1","peak_umap_2")
peak_umap <- select(peak_umap, 
                    archr_name, original_barcode, 
                    archr_cluster,
                    peak_umap_1, peak_umap_2)

write.csv(peak_umap,
          "peak_umap_coords.csv",
          row.names = FALSE)

tile_umap <- getEmbedding(proj, "UMAP", returnDF = TRUE) %>%
  as.data.frame() %>%
  mutate(archr_name = rownames(.)) %>%
  mutate(original_barcode = sub(".+#","",archr_name)) %>%
  mutate(archr_cluster = archr_meta$Clusters)
names(tile_umap)[1:2] <- c("tile_umap_1","tile_umap_2")
tile_umap <- select(tile_umap, 
                    archr_name, original_barcode, 
                    archr_cluster,
                    tile_umap_1, tile_umap_2)

write.csv(tile_umap,
          "tile_umap_coords.csv",
          row.names = FALSE)

### Integrative analysis

seRNA <- import10xFeatureMatrix(
  input = "/mnt/x060-multiome/X061-pipeline-analysis/rna_preprocessed/X061-P0_X061-Well1_labeled.h5",
  names = "X061-AP0C1W1"
)

rna_meta <- H5weaver::read_h5_cell_meta("/mnt/x060-multiome/X061-pipeline-analysis/rna_preprocessed/X061-P0_X061-Well1_labeled.h5")
rownames(rna_meta) <- paste0("X061-AP0C1W1#",paste0(rna_meta$original_barcodes,"-1"))

colnames(seRNA) <- rownames(rna_meta)
colData(seRNA) <- as(rna_meta, "DataFrame")

common_bcs <- intersect(getCellNames(proj), rownames(rna_meta))

seRNA <- seRNA[,common_bcs]

proj <- addGeneIntegrationMatrix(proj,
                                 seRNA = seRNA,
                                 groupRNA = "seurat_pbmc_type",
                                 force = TRUE)

proj <- addGeneExpressionMatrix(input = proj,
                                seRNA = seRNA,
                                force = TRUE)

## transfer some info to fool ArchR into using GeneExpressionMatrix
library(rhdf5)
arrowFile <- getArrowFiles(proj)
predscore <- h5read(arrowFile, "/GeneIntegrationMatrix/Info/predictionScore")
predscore2 <- rep(1, length(predscore))
h5write(predscore2, arrowFile, "/GeneExpressionMatrix/Info/predictionScore")

proj <- addCoAccessibility(
  ArchRProj = proj,
  reducedDims = "IterativeLSI"
)

cA <- getCoAccessibility(
  proj,
  corCutOff = 0.7,
  resolution = 1,
  returnLoops = FALSE
)

peaks_gr <- getPeakSet(proj)
peakMatrix <- assay(getMatrixFromProject(proj, "PeakMatrix"))

get_linked_peak_scores <- function(peaks_gr,
                                   peakMatrix,
                                   geneSymbol) {
  
  prom_idx <- which(mcols(peaks_gr)$nearestGene == geneSymbol & mcols(peaks_gr)$peakType == "Promoter")
  prom_cA <- cA[cA$queryHits %in% cd8_prom_idx | cA$subjectHits %in% cd8_prom_idx,]
  linked_peak_idx <- unique(c(prom_cA$queryHits, prom_cA$subjectHits))
  
  linked_peak_mat <- peakMatrix[linked_peak_idx,]
  colSums(linked_peak_mat)
}

cd8_scores <- get_linked_peak_scores(peaks_gr,
                                     peakMatrix,
                                     "CD8A")

proj <- addPeak2GeneLinks(
  proj,
  reducedDims = "IterativeLSI",
  useMatrix = "GeneExpressionMatrix"
  )

p2g <- metadata(proj@peakSet)$Peak2GeneLinks

proj <- addCellColData(proj,
                       cd8_scores,
                       "CD8_score",
                       cells = names(cd8_scores),
                       force = TRUE)

plotEmbedding(proj,
              colorBy = "cellColData",
              name = "CD8_score")

plotEmbedding(proj,
              colorBy = "GeneExpressionMatrix",
              name = "CD8A")

plotEmbedding(proj,
              colorBy = "GeneScoreMatrix",
              name = "CD8A")

gs_mat <- getMatrixFromProject(proj, "GeneScoreMatrix")
saveRDS(gs_mat, "atac_GeneScore_matrix.rds")

library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(patchwork)
library(immutils)

ref_so <- readRDS("reference/pbmc_multimodal.rds")

peak_mat <- readRDS("atac_peak_mat.rds")

tenx_res <- Read10X_h5("/mnt/x060-multiome/X061-Well1/outs/filtered_feature_bc_matrix.h5")
tenx_mat <- tenx_res$`Gene Expression`

adt_mat <- readRDS("adt_count_mat.rds")
adt_mat <- Matrix::t(as(adt_mat, "dgCMatrix"))
adt_mat <- adt_mat[!grepl("Control", rownames(adt_mat)),]
adt_mat <- adt_mat[!rownames(adt_mat) %in% c("CD10", "CD127", "CD197", "CD25", "CD278", "CD279", "CD304", "CD319", "CD80", "TCR-Va24-Ja18", "TCR-Va7.2"),]

common_cells <- intersect(colnames(peak_mat), colnames(tenx_res$`Gene Expression`))
common_cells <- intersect(common_cells, colnames(adt_mat))

so <- CreateSeuratObject(counts = tenx_mat[,common_cells],
                         project = "X061-Multiome")
so[["ADT"]] <- CreateAssayObject(counts = adt_mat[,common_cells])
so[["ATAC"]] <- CreateAssayObject(counts = peak_mat[,common_cells])

# Remove low-content cells
ggplot() +
  geom_point(data = so@meta.data,
             aes(x = nFeature_RNA,
                 y = nFeature_ATAC),
             size = 0.2) +
  geom_hline(aes(yintercept = 5000),
             color = "red", linetype = "dashed") +
  geom_vline(aes(xintercept = 750),
             color = "red", linetype = "dashed")

so <- so[,so$nFeature_RNA > 750]

DefaultAssay(so) <- 'RNA'
so <- NormalizeData(so)
so <- FindVariableFeatures(so)
so <- ScaleData(so)
so <- RunPCA(so)

so <- RunUMAP(so,
              assay = 'RNA',
              reduction.name = "umap",
              dims = 1:30)

DefaultAssay(so) <- 'ADT'
VariableFeatures(so) <- rownames(so[["ADT"]])
so <- NormalizeData(so, 
                    normalization.method = 'LogNormalize', 
                    scale.factor = 1000)
so <- ScaleData(so) 
so <- RunPCA(so,
             npcs = 20,
             reduction.name = 'adt_pca')
so <- RunUMAP(so,
              dims = 1:12,
              reduction = 'adt_pca',
              assay = 'ADT',
              reduction.name = 'adt_umap')

DefaultAssay(so) <- 'ATAC'
so <- RunTFIDF(so)
so <- FindTopFeatures(so, min.cutoff = 'q75')
so <- RunSVD(so,
             reduction.name = "atac_lsi")
so <- RunUMAP(so,
              dims = 1:30,
              reduction = 'atac_lsi',
              assay = 'ATAC',
              reduction.name = 'atac_umap')

# Label transfer via RNA
DefaultAssay(so) <- 'RNA'
so <- SCTransform(so)
anchors <- FindTransferAnchors(
  reference = ref_so,
  query = so,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)
so <- MapQuery(
  anchorset = anchors,
  query = so,
  reference = ref_so,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

saveRDS(so, "joint_seurat_object.rds")

ref_umap_df <- as.data.frame(so@reductions$ref.umap@cell.embeddings)
ref_umap_df$celltype.l2 <- so$predicted.celltype.l2

rna_umap_df <- as.data.frame(so@reductions$umap@cell.embeddings)
rna_umap_df$celltype.l2 <- so$predicted.celltype.l2

atac_umap_df <- as.data.frame(so@reductions$atac_umap@cell.embeddings)
atac_umap_df$celltype.l2 <- so$predicted.celltype.l2

adt_umap_df <- as.data.frame(so@reductions$adt_umap@cell.embeddings)
adt_umap_df$celltype.l2 <- so$predicted.celltype.l2

p1 <- ggplot(rna_umap_df) +
  geom_point(aes(x = UMAP_1,
                 y = UMAP_2,
                 color = celltype.l2),
             size = 0.2) +
  scale_color_varibow() +
  large_guides()
p2 <- ggplot(ref_umap_df) +
  geom_point(aes(x = refUMAP_1,
                 y = refUMAP_2,
                 color = celltype.l2),
             size = 0.2) +
  scale_color_varibow() + 
  large_guides()
p3 <- ggplot(atac_umap_df) +
  geom_point(aes(x = atac_umap_1,
                 y = atac_umap_2,
                 color = celltype.l2),
             size = 0.2) +
  scale_color_varibow() +
  large_guides()
p4 <- ggplot(adt_umap_df) +
  geom_point(aes(x = adt_umap_1,
                 y = adt_umap_2,
                 color = celltype.l2),
             size = 0.2) +
  scale_color_varibow() + 
  large_guides()

(p1 + p2) / (p3 + p4)


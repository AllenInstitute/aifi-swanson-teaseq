library(Seurat)

mats <- Read10X_h5("/mnt/x059-multiome/X059_cellranger-arc/X059-W1/outs/filtered_feature_bc_matrix.h5")

so <- CreateSeuratObject(mats$`Gene Expression`)

so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")

VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

so <- NormalizeData(so)
so <- FindVariableFeatures(so)
so <- ScaleData(so)

so <- RunPCA(so, features = VariableFeatures(object = so))

DimPlot(so, reduction = "pca")
ElbowPlot(so)

so <- FindNeighbors(so, dims = 1:15)
so <- FindClusters(so, resolution = 0.5)


so <- RunUMAP(so, dims = 1:15)

DimPlot(so, reduction = "umap")

VlnPlot(so, features = c("MS4A1", "CD79A"))

FeaturePlot(so, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "IL7R", 
                               "CD8A"))

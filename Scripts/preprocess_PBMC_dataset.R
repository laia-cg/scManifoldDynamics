# ===================
# PBMC 8k Preprocessing Script
#
# Dependencies:
# -------------------------------------------------------------------------------
# install.packages("BiocManager")  # (if it's not already installed)
# BiocManager::install(version = "3.16")

# CRAN packages:
# install.packages("dplyr")
# install.packages("Seurat")
# install.packages("patchwork")

# Bioconductor packages:
# BiocManager::install(c("celldex", "SingleR", "SingleCellExperiment"))

# In case you want to run UMAP, you also need:
# reticulate::py_install("umap-learn") 

# ===================
# Load Libraries
# ===================
library(dplyr)
library(Seurat)
library(patchwork)
library(celldex)
library(SingleR)
library(SingleCellExperiment)

# ===================
# 1. Load the PBMC dataset
# ===================
# Please replace the path accordingly
pbmc.data <- Read10X(
  data.dir = "/pbmc8k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/GRCh38/"
)

# ===================
# 2. Create a Seurat object
# ===================
pbmc <- CreateSeuratObject(
  counts = pbmc.data,
  project = "pbmc8k",
  min.cells = 3,
  min.features = 200
)
pbmc

# ===================
# 3. Calculate mitochondrial percentage and visualize QC metrics
# ===================
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# ===================
# 4. Subset cells based on QC thresholds, normalize, find variable features
# ===================
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 200)

# Scale all genes
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# ===================
# 5. Perform PCA and find neighbors/clusters
# ===================
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:16)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Uncomment if you have UMAP installed (via 'umap-learn' in Python):
# pbmc <- RunUMAP(pbmc, dims = 1:10, umap.method = 'umap-learn')

# ===================
# 6. Marker Identification
# ===================
pbmc.markers <- FindAllMarkers(pbmc,
                               only.pos = TRUE,
                               min.pct = 0.25,
                               logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)

# ===================
# 7. Rename cluster IDs
# ===================
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T",
                     "FCGR3A+ Mono", "NK", "DC", "Platelet", rep("Unknown", 7))
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

pbmc8k <- pbmc
pbmc8k_hvg <- VariableFeatures(pbmc8k)

# ===================
# 8. Annotation with SingleR
# ===================
ref <- BlueprintEncodeData()
pbmc8k_sce <- as.SingleCellExperiment(pbmc8k)
pred <- SingleR(test = pbmc8k_sce, ref = ref, labels = ref$label.main)
table(pred$labels)

colLabels(pbmc8k_sce) <- pred$pruned.labels
pbmc8k_clust_annotation <- colLabels(pbmc8k_sce)

# ===================
# 9. Save data 
# ===================
pbmc8k_counts <- counts(pbmc8k_sce[pbmc8k_hvg, ])
pbmc8k_log_counts <- logcounts(pbmc8k_sce[pbmc8k_hvg, ])

write.csv(pbmc8k_counts, file = "check_PBMC8k_counts.csv")
write.csv(pbmc8k_log_counts, file = "check_PBMC8k_log_counts.csv")
write.csv(pbmc8k_clust_annotation, file = "check_PBMC8k_clust_annotation.csv")



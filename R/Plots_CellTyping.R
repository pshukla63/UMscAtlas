# this script was used to for displaying of final typed dataset 
# Figures 1d-g, 2 a-b, 3b,
# Supplementary Figures 1e, 3a, 7e-f

# load libraries
library(pals)
library(grDevices)
library(paletteer)
library(colorspace)
library(brew)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplotify)
library(ggplot2)
library(qs)
library(rlang)
library(harmony)

# color palettes used colors saved in colors.R

# load datasets (full datset)
seu <- qread("../final_Dataset.qs")

# prepare subsets
# melanoma/melanocyte subset (not integrated)
Mela <- subset(seu, subset = (majority_celltype %in% c("Melanocytes", "Melanoma cells")))

# subsetting all T&NK Cells
tnk_cells <- c("CD4_EM", "CD4_EX", "CD4_N", "CD4_REG", "CD8_EM", "CD8_EMRA", "CD8_EX", "CD8_N", "CD8_RM", "gdT", "MAIT", "NK_CYTO", "NK_REST", "Proliferating cells")
TNK <- subset(seu, subset = SubTyping %in% tnk_cells)

# integration of T/NK cell subset
# Harmony integration
TNK <- RunHarmony(
  object = TNK,
  group.by.vars = c("Condition", "orig.ident"),
  reduction.use = "pca",
  dims.use = 1:50,
  assay.use = "RNA",
  theta = c(2, 2)
)

# re-compute clustering and UMAP after integration
TNK <- RunUMAP(TNK, dims = 1:30, reduction = "harmony", verbose = FALSE)
TNK <- FindNeighbors(TNK, dims = 1:30, reduction = "harmony", verbose = FALSE)
# use range of clustering resolutions
TNK <- FindClusters(TNK, resolution = c(0, 0.01, 0.2, 0.4, 0.6, 0.8, 1), verbose = FALSE)

# Figure 1d
Seurat::DimPlot(seu, reduction = "umap", group.by = "majority_celltype", raster = FALSE, cols = majority_celltype_col, shuffle = TRUE, label = FALSE)

# Differential expression of T/NK cells for Figure 1e dotplot
order_celltypes <- c("B/Plasma cells", "Melanocytes", "Melanoma cells", "Endothelial cells", "Pericytes", "Glial cells", "Myeloid cells", "T/NK cells",   "Smooth muscle cells", "Mesenchymal cells (healthy)", "Epithelial cells", "Fibroblasts", "Photoreceptor cells")

de_genes <- Seurat::FindAllMarkers(seu,  min.pct = 0.25,
                                   only.pos = TRUE, group.by = "majority_celltype")

de_genes <- subset(de_genes, de_genes$p_val_adj < 0.05)

top_specific_markers <- de_genes %>%
  group_by(cluster) %>%
  top_n(2, avg_log2FC)

# set factor levels in the Seurat object metadata
seu$majority_celltype <- factor(seu$majority_celltype, levels = order_celltypes)

markers_final <- c("IGHM", "MS4A1", "MLANA", "TYR", "PLVAP", "VWF", "PECAM1", "HIGD1B", "FHL5", "MBP", "MAG", "IL1B", "LYZ", "GZMK", "ICOS", "CD2", "FCGR3A", "SFRP4", "MCAM", "DES", "ACTA2", "TAGLN",  "NKX2-1", "OTOS", "KRT7", "COL1A1", "COL1A2", "DCN", "GNAT1", "RAX2", "PRAME", "PMEL")

# Figure 1e
dittoSeq::dittoDotPlot(seu, vars = markers_final, group.by = "majority_celltype", min = -3, max =  3)

# Figure 1f
Seurat::DimPlot(seu, reduction = "umap", group.by = "Malignant", raster = FALSE, cols = iCNV_colors, shuffle = TRUE)

# Figure 1g
Seurat::DimPlot(seu, reduction = "umap", group.by = "secondary_mutation", raster = FALSE, cols = second_mutations_col, shuffle = TRUE)


# Supplementary Figure 1e
Seurat::DimPlot(seu, reduction = "umap", group.by = "orig.ident", raster = FALSE, cols = sample_col, shuffle = TRUE)

# Figure 2a
Seurat::DimPlot(seu, reduction = "umap", group.by = "location", raster = FALSE, cols = location_colors, shuffle = TRUE)

# Figure 2b
Seurat::FeaturePlot(seu, features = "MLANA", cols = c("lightgrey", "turquoise4"), reduction = "umap", raster = FALSE)

# Figure 3b
#UMAP colored by origin
desired_order <- c("healthy", "primary", "metastasis")
Seurat::DimPlot(Mela, reduction = "umap", group.by = "origin", raster = FALSE, cols = origin_colors, shuffle = TRUE)

# Supplementary Figure 3a
Seurat::DimPlot(seu, reduction = "umap", group.by = "primary_mutation", raster = FALSE, cols = primary_mutations_col, shuffle = TRUE)

# Supplementary Figure 7e
Seurat::DimPlot(TNK, reduction = "umap", group.by = "SubTyping", raster = FALSE, cols = TME_SubTyp_col, shuffle = TRUE)

# Differential expression of T/NK cells for # Supplementary Figure 7f dotplot
de_genesTNK <- Seurat::FindAllMarkers(TNK,  min.pct = 0.25,
                                      only.pos = TRUE, group.by = "SubTyping")

de_genesTNK <- subset(de_genesTNK, de_genesTNK$p_val_adj < 0.05)

top_specific_markersTNK <- de_genesTNK %>%
  group_by(cluster) %>%
  top_n(2, avg_log2FC)

# order for dotplot
order_TNKcelltypes <- c("CD8_EX",  "CD8_EM", "NK_CYTO", "CD8_EMRA", "CD4_N", "CD8_N","NK_REST", "CD4_EX", "CD4_REG", "Proliferating cells", "CD4_EM",  "CD8_RM", "gdT", "MAIT")

markers_TNK <- c("LAG3", "GZMK", "GZMB",  "NKG7", "CCL5", "TRDC", "TRGC1", "TRGV9", "KLRD1", "CD69", "ITGAE", "CXCR6", "ZNF683",  "FGFBP2", "PRF1", "GZMH", "KLRG1", "CCR7", "IL7R", "TCF7", "LEF1", "FCGR3A", "KLRC1", "XCL1", "XCL2", "PDCD1", "CTLA4", "HAVCR2", "TOX", "CCR6", "CXCR3", "SELL",  "FOXP3", "IL2RA", "TIGIT", "MKI67", "TOP2A", "CCNB1", "TRAV1-2", "KLRB1")

# set factor levels in the Seurat/SCE object metadata
TNK$SubTyping <- factor(TNK$SubTyping, levels = order_TNKcelltypes)

# Supplementary Figure 7f
dittoSeq::dittoDotPlot(TNK, vars = markers_TNK, group.by = "SubTyping", min = -2, max =  2)

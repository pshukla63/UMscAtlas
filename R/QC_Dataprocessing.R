# this script was used to for raw processing including normalizing, scaling and clustering of the data
# Supplementary Figures 1a-d

# packages
library(qs)
library(Seurat)

# define colors (see "color_palettes.Rmd")
# primary_samp_col defines colors for each entry in orig.ident for primary tumors(samples with same orig.ident were obtained from same tumor biopsy and processed in duplicates if there is more than one samplename entry in dataset)
pUM_samples <- c("pUM01", "pUM02", "pUM03", "pUM04", "pUM05", "pUM06", "pUM09", "pUM10", "pUM11", "pUM12", "pUM14", "pUM15", "pUM16", "pUM17", "pUM19", "pUM20", "pUM21", "pUM22", "pUM23", "pUM24", "pUM25")
mint_colors <- as.character(paletteer::paletteer_c("grDevices::Mint", n = 24))
cats <- unique(pUM_samples)
if (length(cats) > 24) warning("More than 21 categories: colors will recycle")

# Map categories to colors (up to 21)
offset <- 3
col_map <- setNames(mint_colors[ (seq_along(cats) - 1 + offset) %% 24 + 1 ], cats)

# Assign colors to each element of x
primary_samp_col <- col_map[pUM_samples]


# mets_samp_col defines colors for each entry in orig.ident for metastatic tumors(samples with same orig.ident were obtained from same tumor biopsy and processed in duplicates if there is more than one samplename entry in dataset)
mUM_samples <- c("mUM01", "mUM02", "mUM04", "mUM05", "mUM06", "mUM07", "mUM08", "mUM09", "mUM11", "mUM13", "mUM14", "mUM15", "mUM16", "mUM18", "mUM19", "mUM20")

# Get all original colors (the full discrete palette)
acton_all <- as.character(paletteer::paletteer_d("khroma::acton"))

# Build an interpolator over the full palette and sample 16 evenly spaced colors
acton_16 <- grDevices::colorRampPalette(acton_all)(16)

# Print each color with its index
#for (i in seq_along(acton_16)) {
#  cat(sprintf("%2d: %s\n", i, acton_16[i]))}

# Choose a stable order for unique samples (customize if you need a specific order)
sample_levels <- unique(mUM_samples)  # or: sort(unique(sample_names))
if (length(sample_levels) > 16) {
  stop("You have more than 16 unique samples; increase n or reduce categories.")
}

# Map sample -> color (1:1, in the chosen order)
col_map_M <- setNames(acton_16[seq_along(sample_levels)], sample_levels)

# Per-element colors aligned to sample_names
mets_samp_col <- col_map_M[mUM_samples]


# sample_col combines healthy_samp_col, primary_samp_col and mets_samp_col to color all different tumor biopsies sequenced
sample_col <- c(healthy_samp_col, primary_samp_col, mets_samp_col)
# load full data set
seu <- qread("../dataset_merged_seurat_metadata.rds")

# check scatter plot and violine plot to visualize QC per cell and gene, save as PDF (first QC)
DefaultAssay(seu) <- "RNA"
seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
Seurat::FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", raster = FALSE)


seu <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mito")
seu <- Seurat::PercentageFeatureSet(seu, pattern = "^MRP[SL]", col.name = "percent.ribo")

Seurat::VlnPlot(seu, features = c("nCount_RNA", "nFeature_RNA"), group.by = "orig.ident", pt.size = 0.1)
Seurat::VlnPlot(seu, features = c("nFeature_RNA", "percent.mito"), group.by = "orig.ident", pt.size = -1)

# remove cells
remove cells with very low features (keep only features >300 and <8000, % mito <50% and Counts < 100000)

seu <- subset(seu, subset = nFeature_RNA > 300 & 
                  nFeature_RNA < 8000 &
                  percent.mito < 50 &
                  nCount_RNA < 100000)

# Supplementary Figure 1a-d
# Supplementary Figure 1a
Seurat::VlnPlot(seu, features = "nCount_RNA", group.by = "orig.ident", assay = NULL, raster = FALSE, cols = sample_col, pt.size = 0.1, alpha = 0.03) + 
  theme(text = element_text(size = 25), axis.text  = element_text(size = 20), axis.title = element_text(size = 20)
  )
# Supplementary Figure 1b
Seurat::VlnPlot(seu, features = "nFeature_RNA", group.by = "orig.ident", raster = FALSE, cols = sample_col, pt.size = 0.1, alpha = 0.03) + 
  theme(text = element_text(size = 25), axis.text  = element_text(size = 20), axis.title = element_text(size = 20)
  )

# Supplementary Figure 1c
Seurat::VlnPlot(seu, features = "percent.mito", group.by = "orig.ident", raster = FALSE, cols = sample_col, pt.size = 0.1, alpha = 0.03) + 
  theme(text = element_text(size = 25), axis.text  = element_text(size = 20), axis.title = element_text(size = 20)
  )

# Supplementary Figure 1d
Seurat::VlnPlot(seu, features = "percent.ribo", group.by = "orig.ident", raster = FALSE, cols = sample_col, pt.size = 0.1, alpha = 0.03) + 
  theme(text = element_text(size = 25), axis.text  = element_text(size = 20), axis.title = element_text(size = 20)
  )


# Normalizing and find Variable Features
Normalize dataset & find variable features & check features (add few extra genes of interest for UM biology: PTEN, CCR4, PMEL, MLANA, MITF, SOX10, TYR)

seu <- Seurat::NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)


seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)

# extra genes of interest
my_genes <- c("PTEN", "CCR4", "PMEL", "MLANA", "MITF", "SOX10", "TYR")

# Combine with existing HVGs
current_hvgs <- Seurat::VariableFeatures(seu)
combined_hvgs <- unique(c(current_hvgs, my_genes))

# Set the new list of HVGs
Seurat::VariableFeatures(seu) <- combined_hvgs

# Get HVGs list in hvgs
hvgs <- Seurat::VariableFeatures(seu)

# Check which of the genes are in the HVGs
my_genes_in_hvgs <- my_genes %in% hvgs

# Print results
data.frame(Gene = my_genes, In_HVGs = my_genes_in_hvgs)

length(VariableFeatures(seu))

# Scale Data
seu <- Seurat::ScaleData(seu)

# run PCA dimensionality reduction and chose with cumu PCA where over 90% of variability is described

DefaultAssay(seu) <- "RNA"

seu <- Seurat::RunPCA(seu)

Seurat::DimPlot(seu, reduction = "pca")

Seurat::DimHeatmap(seu, dims = 1:12, cells = 500, balanced = TRUE)

Seurat::ElbowPlot(seu, ndims = 40)

cumu <- cumsum((seu[["pca"]]@stdev)^2 / sum((seu[["pca"]]@stdev)^2) * 100)
npc <- which(cumu > 90)[1]
seu <- RunUMAP(seu, dims = 1:npc)


# UMAP
run UMAP dimensionality reduction and save new object

Seurat:: DimPlot(seu, reduction = "umap", group.by = "orig.ident")


# Clustering
clustering, find neigbors using same dimensions as for UMAP generation run clustering in 07_01_Clustering (already removed outliers etc. before in first round of analysis!)

seu <- FindNeighbors(seu, dims = 1:npc)

seu <- Seurat::FindClusters(seu, resolution = 0.4)

qsave(seu, "../dataset_merged_seurat_Clustering_seq_04.qs")

# subset data set depending on organ to run reference-based cell typing for each organ (eye, brain, liver, skin and thyroid)
eyes <- subset(seu, subset = organ == "eye")
brain <- subset(seu, subset = organ == "brain")
liver <- subset(seu, subset = organ == "liver")
skin <- subset(seu, subset = organ == "skin")
thyroid <- subset(seu, subset = organ == "thyroid")


# save data
qsave(eyes, "../dataset_Eyes_seurat_forCelltyping.qs")
qsave(brain, "../dataset_Brain_seurat_forCelltyping.qs")
qsave(liver, "../dataset_Liver_seurat_forCelltyping.qs")
qsave(skin, "../dataset_Skin_seurat_forCelltyping.qs")
qsave(thyroid, "../dataset_Thyroid_seurat_forCelltyping.qs")


sessioninfo::session_info()
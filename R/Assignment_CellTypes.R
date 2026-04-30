# this script was used to assing cell types
# It shows eye subset as example for all 5 different organ subsets

# packages
library(qs)
library(Seurat)
library(dplyr)
library(readr)
library(tibble)
library(Azimuth)
library(tidyverse)

# define colors (see "color_palettes.Rmd")

# load data sets for cell typing
qsave(eyes, "../dataset_Eyes_seurat_forCelltyping.qs")
qsave(brain, "../dataset_Brain_seurat_forCelltyping.qs")
qsave(liver, "../dataset_Liver_seurat_forCelltyping.qs")
qsave(skin, "../dataset_Skin_seurat_forCelltyping.qs")
qsave(thyroid, "../dataset_Thyroid_seurat_forCelltyping.qs")


ref_seurat <- qs_read("../eye_ref_downsampled.qs2")


# cell typing using Seurat

DefaultAssay(ref_seurat) <- "RNA"  # or the appropriate assay name
DefaultAssay(eyes) <- "RNA"  # or the appropriate assay name

ref_seurat <- FindVariableFeatures(ref_seurat, selection.method = "vst", nfeatures = 3000)

# extra genes of interest
my_genes <- c("PTEN", "CCR4", "PMEL", "MLANA", "MITF", "SOX10", "TYR")

# Combine with existing HVGs
current_hvgs <- Seurat::VariableFeatures(ref_seurat)
combined_hvgs <- unique(c(current_hvgs, my_genes))

# Set the new list of HVGs
Seurat::VariableFeatures(ref_seurat) <- combined_hvgs

# Get HVGs list in hvgs
hvgs <- Seurat::VariableFeatures(ref_seurat)

# Check which of the genes are in the HVGs
my_genes_in_hvgs <- my_genes %in% hvgs

common_features <- intersect(VariableFeatures(ref_seurat), VariableFeatures(eyes))

# Find anchors between reference and query datasets
anchors <- FindTransferAnchors(reference = ref_seurat, query = eyes, dims = 1:50, features = common_features)

# Transfer labels from reference to query dataset
predictions <- TransferData(anchorset = anchors, refdata = ref_seurat$author_cell_type, dims = 1:50)
eyes <- AddMetaData(eyes, metadata = predictions)

# Visualize results
DimPlot(eyes, reduction = "pca", group.by = "predicted.id")
# Visualize results
DimPlot(eyes, reduction = "umap", group.by = "predicted.id", raster = FALSE, label = TRUE, repel = TRUE, cols = 'polychrome')

# save data
qsave(eyes, "../dataset_Eyes_seurat_RefBasedCelltyping.qs")


# extract metadata for import in full Dataset
# Define new columns to export (make sure they exist in your Seurat object)
new_columns <- c("predicted.id", "prediction.score.oligodendrocyte.precursor.cell", "prediction.score.astrocyte", "prediction.score.oligodendrocyte", "prediction.score.fibroblast", "prediction.score.microglial.cell", "prediction.score.vascular.associated.smooth.muscle.cell", "prediction.score.macrophage", "prediction.score.pericyte", "prediction.score.endothelium", "prediction.score.natural.killer.cell", "prediction.score.T.cell", "prediction.score.melanocyte", "prediction.score.myelinating.Schwann.cell", "prediction.score.B.cell", "prediction.score.Schwann.cell", "prediction.score.pigmented.epithelial.cell", "prediction.score.Melanoma", "prediction.score.max")  # example new columns

# Extract metadata and add cell IDs
metadata_subset <- merge@meta.data[, new_columns] %>%
  rownames_to_column("cell_id")

# Save to CSV
write_csv(metadata_subset, "../metadataEyeSubset.csv")

# load full dataset and metadate from cell type assignment
seu <- qread("../dataset_merged_seurat_Clustering_seq_04.qs")

eyes_meta <- read.csv("../metadataEyeSubset.csv")
skin_meta <- read.csv("../metadataSkinSubset.csv")
liver_meta <- read.csv("../metadataLiverSubset.csv")
thyroid_meta <- read.csv("../metadataThyroidSubset.csv")
brain_meta <- read.csv("../metadataBrainSubset.csv")


# Load the Excel file with how to homogenize cellnames
cell_info <- read_excel("../Homogenizing_cellnames.xlsx")

# merge dataset into one dataset for later clusterbased celltyping and merging of celltypes from reference based celltyping and save merged object
# Merge all metadata frames by cell_id
merged_meta <- purrr::reduce(
  list(eyes_meta, skin_meta, liver_meta, thyroid_meta, brain_meta),
  function(x, y) full_join(x, y, by = "cell_id")
)

# Set cell_id as rownames
merged_meta <- column_to_rownames(merged_meta, var = "cell_id")

# Collapse all predicted.id columns into one combined column
predicted_cols <- grep("predicted.id", colnames(merged_meta), value = TRUE)

merged_meta$predicted.id_combined <- apply(
  merged_meta[, predicted_cols], 1,
  function(x) paste(na.omit(x), collapse = "; ")
)

# Remove only the original predicted.id columns (not the new combined one)
merged_meta <- merged_meta[, !(colnames(merged_meta) %in% predicted_cols)]


# add simplified cell type names to dataset setup for renaming from excel dataset
# Create a named vector: names are Original Cell Types, values are Simplified Names
mapping_vector <- setNames(cell_info$`Simplified Name`, cell_info$`Original Cell Type`)

# Extract first matching cell type from predicted.id_combined
merged_meta$homogenized_cellnames <- sapply(merged_meta$predicted.id_combined, function(x) {
  types <- unlist(strsplit(x, "; "))
  match <- types[types %in% names(mapping_vector)]
  if (length(match) > 0) {
    return(mapping_vector[match[1]])  # Use first match
  } else {
    return(NA)
  }
})

seu <- AddMetaData(seu, metadata = merged_meta)

# cluster based cell type assignment add new cell type column (if 80% of all cells of a cluster are one cell type add this cell type to the new column)
# Ensure cluster identity is set
Idents(seu) <- seu$RNA_snn_res.0.4

# Create a data frame with cluster and cell type info
df <- data.frame(
  Cluster = Idents(seu),
  CellType = seu$homogenized_cellnames
)

# Calculate percentage of each cell type per cluster
celltype_distribution <- df %>%
  group_by(Cluster, CellType) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Cluster) %>%
  mutate(Fraction = Count / sum(Count)) %>%
  ungroup()

# Identify clusters where one cell type exceeds 80%
dominant_labels <- celltype_distribution %>%
  group_by(Cluster) %>%
  filter(Fraction == max(Fraction)) %>%
  filter(Fraction > 0.8) %>%
  ungroup()

# Create a named vector for mapping
cluster_to_label <- setNames(dominant_labels$CellType, dominant_labels$Cluster)

# Assign new labels to cells based on cluster
new_labels <- as.character(Idents(seu))
new_labels <- ifelse(new_labels %in% names(cluster_to_label),
                     cluster_to_label[new_labels],
                     "Mixed")

# Add new column to Seurat object
seu$majority_celltype <- new_labels

# assing mixed clusters manually 

# Ensure cluster identity is set correctly
Idents(seu) <- seu$RNA_snn_res.0.4

# Copy the current majority cell type column
new_labels <- seu$majority_celltype

# Manually update specific clusters
new_labels[Idents(seu) == "9"]  <- "Pericytes"
new_labels[Idents(seu) == "21"] <- "Mesenchymal cells (healthy)"
new_labels[Idents(seu) == "25"] <- "Melanocytes"
new_labels[Idents(seu) == "28"] <- "Astrocytes"

# Assign updated labels back to the Seurat object
seu$majority_celltype <- new_labels

qsave(seu, "../dataset_mergedAll_seurat_ClusterBasedCelltyping.qs")


# cell typing using Azimuth
# subset dataset by predicted celltype Immune cells
seu_Immunecells <- subset(seu, subset = majority_celltype %in% c("B/Plasma cells", "Myeloid cells", "T/NK cells"))

# Harmony integration
seu_Immunecells <- RunHarmony(
  object = seu_Immunecells,
  group.by.vars = c("Condition", "orig.ident"),
  reduction.use = "pca",
  dims.use = 1:50,
  assay.use = "RNA",
  theta = c(2, 2)
)

# re-compute clustering and UMAP after integration
seu_Immunecells <- RunUMAP(seu_Immunecells, dims = 1:30, reduction = "harmony", verbose = FALSE)
seu_Immunecells <- FindNeighbors(seu_Immunecells, dims = 1:30, reduction = "harmony", verbose = FALSE)
# use range of clustering resolutions
seu_Immunecells <- FindClusters(seu_Immunecells, resolution = c(0, 0.01, 0.2, 0.4, 0.6, 0.8, 1), verbose = FALSE)

# !!!!!! Ask Mauro for Azimuth script!!!!!!!!!!!
seu_Immunecells <- RunAzimuth(seu_Immunecells, ref = "pbmcref")
# check if it worked
DimPlot(seu_Immunecells, reduction = "umap", group.by = "predicted.celltype.l1") +
  coord_fixed()
DimPlot(seu_Immunecells, reduction = "umap", group.by = "predicted.celltype.l2") +
  coord_fixed()
# Load the Excel file with how to homogenize cell names
cell_info <- read_excel("../Synchronized_cellnames.xlsx")


# Create mapping vector
mapping_vector <- setNames(cell_info$"Simplified Name", cell_info$"Original Cell Type")


# Generate the synchronized cell names
synchronized_names <- sapply(seu_Immunecells$predicted.celltype.l2, function(x) {
  types <- unlist(strsplit(x, ";\\s*"))  # Split on semicolon and optional space
  match <- types[types %in% names(mapping_vector)]
  if (length(match) > 0) {
    return(mapping_vector[match[1]])  # Use first match
  } else {
    return(NA)
  }
})

# Add to metadata using AddMetaData
# Make sure names match the cell names in the Seurat object
names(synchronized_names) <- colnames(seu_Immunecells)

# Add as metadata
seu_Immunecells <- AddMetaData(seu_Immunecells, metadata = synchronized_names, col.name = "synchronized_cellnames")

qsave(seu_Immunecells, "../dataset_Immunesubset_seurat_Azimuth.qs")

# cell typing using external dataset (TILs) 
# Subset T/NK cells from Azimuth cell typing to type using TIL reference

# Define the cell types to keep
celltypes_to_keep <- c("CD8 T cells", "Proliferating cells", "CD4 T cells", "Treg", 
                       "NK", "ILC", "dnT", "gdT", "NK_CD56bright", "CD4 CTL")

# Get cells to keep based on synchronized_cellnames
cells_to_keep <- rownames(seu_Immunecells@meta.data)[
  seu_Immunecells@meta.data$synchronized_cellnames %in% celltypes_to_keep
]

# Subset the Seurat object
seu_TNK <- subset(seu_Immunecells, cells = cells_to_keep)
# Harmony integration
seu_TNK <- RunHarmony(
  object = seu_Immunecells,
  group.by.vars = c("Condition", "orig.ident"),
  reduction.use = "pca",
  dims.use = 1:50,
  assay.use = "RNA",
  theta = c(2, 2)
)

# re-compute clustering and UMAP after integration
seu_TNK <- RunUMAP(seu_TNK, dims = 1:30, reduction = "harmony", verbose = FALSE)
seu_TNK <- FindNeighbors(seu_TNK, dims = 1:30, reduction = "harmony", verbose = FALSE)
# use range of clustering resolutions
seu_TNK <- FindClusters(seu_TNK, resolution = c(0, 0.01, 0.2, 0.4, 0.6, 0.8, 1), verbose = FALSE)

# load dataset for T/NK cell subtyping
tcell_ref <- readRDS("../Tcell_Reference_Seurat.rds")

# Find variable features for the reference dataset
tcell_ref <- FindVariableFeatures(tcell_ref, selection.method = "vst", nfeatures = 2000)

# Identify common features between reference and query datasets
common_features <- intersect(VariableFeatures(tcell_ref), VariableFeatures(seu_TNK))

# Find transfer anchors
anchors <- FindTransferAnchors(reference = tcell_ref, query = seu_TNK, dims = 1:30)

# Transfer labels from reference to query dataset
predictions <- TransferData(anchorset = anchors, refdata = tcell_ref$cellSubType, dims = 1:30)

# Rename the prediction column to avoid overwriting existing metadata
colnames(predictions)[1] <- "predicted.cellSubType"

# Add new predictions to the Seurat object
seu_TNK <- AddMetaData(seu_TNK, metadata = predictions)

qsave(seu_TNK, "../dataset_TNKsubset_seurat_RefBasedCelltyping.qs")


## merge subtyping data into seu Seurat object with all cells

# Prepare Immune_fine for adding TNK subtypng to Immunecell object (copy of seu:Immune)
Immune_fine <- seu_Immune

# Add TNK metadata for cells in common with Immune cells
common_cells_tnk <- intersect(colnames(Immune_fine), colnames(seu_TNK))
cat("Number of common cells between Immune and TIL:", length(common_cells_tnk), "\n")

immune_cols <- colnames(Immune_fine@meta.data)
tnk_cols <- colnames(seu_TNK@meta.data)
new_cols_tnk <- setdiff(tnk_cols, immune_cols)
cat("Adding columns from TIL to Immune:", new_cols_tnk, "\n")

# Create a metadata frame with NA for all Immune cells
tnk_meta_to_add <- data.frame(matrix(NA, nrow = ncol(Immune_fine), ncol = length(new_cols_tnk)))
rownames(tnk_meta_to_add) <- colnames(Immune_fine)
colnames(tnk_meta_to_add) <- new_cols_tnk

# Fill in values for common cells
tnk_meta_to_add[common_cells_tnk, ] <- seu_TNK@meta.data[common_cells_tnk, new_cols_tnk, drop = FALSE]

# Add metadata to Immune_fine
Immune_fine <- AddMetaData(Immune_fine, metadata = tnk_meta_to_add)


# Create new column called "SubTyping"
final_subtyping_immune <- rep(NA, ncol(Immune_fine))
names(final_subtyping_immune) <- colnames(Immune_fine)

# Priority 1: TNK SubTyping for cells in common with immune cells
if ("SubTyping" %in% colnames(seu_TNK@meta.data)) {
  tnk_subtyping <- seu_TNK@meta.data$SubTyping
  names(tnk_subtyping) <- colnames(seu_TNK)
  common_tnk <- intersect(names(final_subtyping_immune), names(tnk_subtyping))
  final_subtyping_immune[common_tnk] <- tnk_subtyping[common_til]
}

# Priority 2: add Immune synchronized_cellnames for remaining cells into "Subtyping"
if ("synchronized_cellnames" %in% colnames(Immune@meta.data)) {
  immune_subtyping <- Immune@meta.data$synchronized_cellnames
} else {
  stop("Column 'synchronized_cellnames' not found in Immune@meta.data")
}
names(immune_subtyping) <- colnames(Immune)
missing_idx <- is.na(final_subtyping_immune)
final_subtyping_immune[missing_idx] <- immune_subtyping[missing_idx]

# Add SubTyping to Immune_fine and save
Immune_fine@meta.data$SubTyping <- final_subtyping_immune
qsave(Immune_fine, "../dataset_Immunesubset_seurat_RefBasedCelltyping.qs")

# add SubTyping columnt to full dataset "seu" and add majority cell types to SubTyping of all cells not included in Immune seurat object
# extract colnames from TME and seu
cols_Immune <- colnames(Immune_fine@meta.data)
cols_seu <- colnames(seu@meta.data)
#check which columns are only in TME dataset and not full dataset
Immune_not_seu <- setdiff(cols_Immune, cols_seu)
Immune_not_seu

# Intersect cell barcodes
overlap_cells <- intersect(colnames(seu), rownames(Immune_fine@meta.data))

# Build transfer table for overlapping cells
transfer_tbl <- Immune_fine@meta.data[overlap_cells, Immune_fine_not_seu, drop = FALSE]

# Export Metadata incl. cell barcodes of TME in Data TME subset folder
export_tbl <- transfer_tbl
export_tbl$barcode <- rownames(transfer_tbl)
export_tbl <- export_tbl[, c("barcode", setdiff(colnames(export_tbl), "barcode"))]
write.csv(export_tbl, "../TME_only_metadata_with_barcodes_toMerge.csv", row.names = FALSE)

# Merge metadata of overlapping cells into new data set
seu_SubTyping <- AddMetaData(object = seu, metadata = transfer_tbl)
cols_newseu <- colnames(seu_SubTyping@meta.data)
cols_newseu

# save final dataset with all majority_celltype and SubTyping
qsave(merge_filtered, "../final_Dataset.qs")

sessioninfo::session_info()
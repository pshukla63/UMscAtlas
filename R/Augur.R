
# This script was used to create the following figures:
# Figures 4a-g
# Supplementary Figures 7a-d,g

library(tidyverse)
library(Seurat)
library(harmony)
library(qs)
library(pals)
library(Augur)
library(magrittr)
library(clustree)
library(viridis)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
set.seed(1234)

# define colors (see "color_palettes.Rmd")
majority_celltype_col <- c(
  "Glial cells" = "#5A5156", 
  "B/Plasma cells" = "lightslateblue", 
  "Myeloid cells" = "lightblue", 
  "T/NK cells" = "#325A9B", 
  "Endothelial cells" = "chartreuse4", 
  "Epithelial cells" = "#1CBE4F", 
  "Fibroblasts" = "yellowgreen", 
  "Mesenchymal cells (healthy)" = "darkgreen", 
  "Pericytes" = "olivedrab4", 
  "Smooth muscle cells" = "#1CFFCE", 
  "Proliferating cells" = "#001889FF"
  )

second_mutations_col <- c(
  "BAP1" = "#1F867BFF", 
  "SF3B1" = "#EAC541FF", 
  "EIF1AX" = "#7B1C5DFF"
  )

# load full data set
seu <- qread("../full_dataset.qs")

# select only primaries
Idents(seu) <- "origin"
seu <- subset(seu, idents = "primary")

# select primary tumors with known secondary mutations
Idents(seu) <- "secondary_mutation"
seu <- subset(seu, idents = c("none", "unknown"), invert = TRUE)

# run standard Seurat processing for subset
DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, verbose = FALSE)
# define number of PCs that explain at least 90% of variance
cumu <- cumsum((seu[["pca"]]@stdev)^2 / sum((seu[["pca"]]@stdev)^2) * 100)
npc <- which(cumu > 90)[1]
seu <- RunUMAP(seu, dims = 1:npc, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:npc, verbose = FALSE)
seu <- FindClusters(seu, verbose = FALSE)

# confirm that integration is needed
DimPlot(seu, group.by = "orig.ident") +
  coord_fixed()

# Harmony integration
seu <- RunHarmony(
  object = seu,
  group.by.vars = c("Condition", "orig.ident"),
  reduction.use = "pca",
  dims.use = 1:50,
  assay.use = "RNA",
  theta = c(2, 2)
  )

# re-compute clustering and UMAP after integration
seu <- RunUMAP(seu, dims = 1:30, reduction = "harmony", verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:30, reduction = "harmony", verbose = FALSE)
# use range of clustering resolutions
seu <- FindClusters(seu, resolution = c(0, 0.01, 0.2, 0.4, 0.6, 0.8, 1), verbose = FALSE)

# re-label cells from cluster 15 as proliferating
seu$majority_celltype[seu$RNA_snn_res.0.8 == 15] <- "Proliferating cells"

# Fig. 4a
DimPlot(seu, group.by = "majority_celltype", raster = FALSE, shuffle = TRUE) +
  scale_color_manual(values = majority_celltype_col) +
  labs(title = NULL) +
  coord_fixed()

# Suppl. Fig. 7a
DimPlot(seu, group.by = "RNA_snn_res.0.8", label = TRUE) +
  NoLegend() +
  coord_fixed()

# Suppl. Fig. 7b
FeaturePlot(seu, features = "MKI67") +
  coord_fixed()

# Fig. 4b
DimPlot(seu, group.by = "secondary_mutation", raster = FALSE, shuffle = TRUE) +
  coord_fixed() +
  labs(title = NULL) +
  scale_color_manual(values = second_mutations_col)


# run Augur at different clustering resolutions
auc_list <- c()

res_list <- c(0, 0.01, 0.2, 0.4, 0.6, 0.8, 1)

for (resolution in res_list){
  
  temp_auc <- calculate_auc(input = seu,
                            label_col = "secondary_mutation",
                            cell_type_col = paste0("RNA_snn_res.", resolution),
                            n_threads = 16,
                            show_progress = TRUE,
                            min_cells = 100)
  
  auc_list[[paste0("resolution_", resolution)]] <- temp_auc
  
}

# prepare results data frame
auc_list2 <- c()

for (i in 1:length(auc_list)){
  
  auc_list2[[names(auc_list[i])]] <- auc_list[[i]]$AUC
}

cluster_aucs <- bind_rows(auc_list2, .id = "resolution")
cluster_aucs$resolution <- gsub("resolution_", "RNA_snn_res.", cluster_aucs$resolution)

# prepare input for clustering tree visualization
# prepare layout
layout <- clustree(seu, prefix = "RNA_snn_res.", return = "layout", prop_filter = 0.01, show.axis = TRUE)

# add a 'node' column to index this to the tree
cluster_aucs %<>%
  mutate(node = paste0(resolution, "C", cell_type)) %>%
  arrange(node)

# add the AUC in a separate column
stat <- layout %>%
  left_join(cluster_aucs %>%
              dplyr::select(node, auc)) %>%
  mutate(auc = ifelse(is.na(auc), 0.5, auc))

layout$auc = stat$auc[match(layout$node, stat$node)]

# Fig. 4c
ggraph(layout) +
  geom_edge_link(arrow = arrow(length = unit(2, "points"), ends = "last"), 
                 end_cap = circle(5.5 * 1.5, "points"), 
                 start_cap = circle(5.5 * 1.5, "points"), 
                 aes_(alpha = ~ in_prop)) +
  scale_edge_colour_manual(values = "black", guide = "none") +
  scale_edge_alpha(name = "(%)", labels = function(x) x * 100, limits = c(0, 1)) +
  clustree:::add_node_points(node_colour = "auc", node_size = 6, node_alpha = 1, colnames(layout)) +
  scale_color_viridis(name = "AUC", option = "turbo") +
  geom_node_text(aes_(label = ~cluster), size = 4, colour = "black") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.y = element_line(colour = "grey92"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_y_continuous(breaks = unique(layout$y), labels = unique(layout$RNA_snn_res.), name = "Clustering resolution")

# add AUC values to Seurat object's meta data
toadd <- stat %>%
  dplyr::select(RNA_snn_res., cluster, auc) %>%
  dplyr::filter(RNA_snn_res. == 0.8) %>%
  dplyr::select(-RNA_snn_res.) %>%
  dplyr::rename(RNA_snn_res.0.8 = cluster)

reallyadd <- seu@meta.data %>%
  dplyr::select(RNA_snn_res.0.8) %>%
  rownames_to_column("barcode") %>%
  left_join(toadd, by = "RNA_snn_res.0.8") %>%
  column_to_rownames("barcode") %>%
  dplyr::select(auc)

seu <- AddMetaData(seu, metadata = reallyadd)

# Fig. 4d
FeaturePlot(seu, features = "auc") +
  scale_color_viridis(option = "turbo", name = "AUC") +
  coord_fixed() +
  ggtitle("AUC at clustering resolution 0.8")

# Suppl. Fig. 7c
ggplot(seu@meta.data, aes(x = RNA_snn_res.0.8, fill = orig.ident)) +
  geom_bar(position = "fill") +
  labs(x = "Cluster", y = "Fraction") +
  scale_fill_manual(values = primary_samp_col, name = "Sample") +
  theme_bw(base_size = 14)

# Suppl. Fig. 7d
ggplot(seu@meta.data, aes(x = RNA_snn_res.0.8, fill = secondary_mutation)) +
  geom_bar(position = "fill") +
  labs(x = "Cluster", y = "Fraction") +
  scale_fill_manual(values = second_mutations_col, name = "Secondary Mutation") +
  theme_bw(base_size = 14)

# re-run Augur on cell types with fine-grained T/NK cell types added
# for this we can use the column "detailed_celltype"
multiclass_augur <- calculate_auc(
  input = seu,
  label_col = "secondary_mutation",
  cell_type_col = "detailed_celltype",
  n_threads = 16,
  show_progress = TRUE,
  min_cells = 100
)

# make custom lollipop function to have more control over plot aesthetics
custom_lollipop <- function(aucs){
  
  size_sm <- 14
  size_lg <- 16
  range <- range(aucs$auc)
  expand <- abs(diff(range)) * 0.1
  p <- aucs %>%
    ggplot(aes(x = reorder(cell_type, auc), y = auc)) +
    geom_hline(aes(yintercept = 0.5), linetype = 'dotted', size = 0.5) +
    geom_point(size = 1.5) +
    geom_text(aes(label = format(auc, digits = 3),
                  y = ifelse(auc < 0.5, 0.5, auc)), size = 4,
              nudge_y = expand,
              hjust = 0.5) +
    geom_segment(aes(xend = cell_type, yend = 0.5)) +
    scale_y_continuous('AUC', limits = c(min(range[1] - expand, 0.5), 
                                         range[2] + expand * 1.5)) +
    coord_flip() +
    theme_bw() + 
    theme(axis.text.x = element_text(size = size_sm),
          axis.text.y = element_text(size = size_sm),
          axis.title.x = element_text(size = size_lg),
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          strip.text = element_text(size = size_lg),
          strip.background = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_blank(),
          legend.position = "top",
          legend.text = element_text(size = size_sm),
          legend.title = element_text(size = size_sm),
          legend.key.size = unit(0.6, "lines"),
          legend.margin = margin(rep(0, 4)),
          legend.background = element_blank(),
          plot.title = element_text(size = size_lg, hjust = 0.5))
  
  p
  
}

# Fig. 4e
custom_lollipop(multiclass_augur$AUC) +
  theme(plot.margin = margin(5, 20, 5, 5, "pt")) +
  labs(title = "Multiclass comparison")

# create pseudobulk object for differential expression analysis
# remove samples that have too few CD8_EM cells
Idents(seu) <- "orig.ident"
sub_seu <- subset(seu, idents = c("pUM02", "pUM20"), invert = TRUE)

# aggregate expression by sample and celltype
pseudo_seu <- AggregateExpression(sub_seu, assays = "RNA", return.seurat = TRUE, group.by = c("orig.ident", "secondary_mutation", "SubTyping"))
Idents(pseudo_seu) <- "SubTyping"

# Suppl. Fig. 7g
DE_results <- FindMarkers(
  pseudo_seu,
  subset.ident = "CD8-EM",
  group.by = "secondary_mutation",
  ident.1 = "BAP1",
  ident.2 = "EIF1AX",
  assay = "RNA",
  logfc.threshold = -Inf,
  min.pct = -Inf,
  min.diff.pct = -Inf,
  only.pos = FALSE,
  test.use = "DESeq2"
  )

EnhancedVolcano(
  DE_results,
  lab = rownames(DE_results),
  x = "avg_log2FC",
  y = "p_val",
  title = "BAP1 vs. EIF1AX",
  subtitle = "CD8_EM T cells",
  caption = "positive foldchange = upregulated in BAP1",
  drawConnectors = TRUE,
  pCutoff = 0.05,
  xlim = c(-3.0, 3.25),
  ylim = c(0, 10)
  )

DE_results <- FindMarkers(
  pseudo_seu,
  subset.ident = "CD8-EM",
  group.by = "secondary_mutation",
  ident.1 = "BAP1",
  ident.2 = "SF3B1",
  assay = "RNA",
  logfc.threshold = -Inf,
  min.pct = -Inf,
  min.diff.pct = -Inf,
  only.pos = FALSE,
  test.use = "DESeq2"
  )

EnhancedVolcano(
  DE_results,
  lab = rownames(DE_results),
  x = "avg_log2FC",
  y = "p_val",
  title = "BAP1 vs. SF3B1",
  subtitle = "CD8_EM T cells",
  caption = "positive foldchange = upregulated in BAP1",
  drawConnectors = TRUE,
  pCutoff = 0.05,
  xlim = c(-2.25, 4.0),
  ylim = c(0, 7)
  )

DE_results <- FindMarkers(
  pseudo_seu,
  subset.ident = "CD8-EM",
  group.by = "secondary_mutation",
  ident.1 = "EIF1AX",
  ident.2 = "SF3B1",
  assay = "RNA",
  logfc.threshold = -Inf,
  min.pct = -Inf,
  min.diff.pct = -Inf,
  only.pos = FALSE,
  test.use = "DESeq2"
  )

EnhancedVolcano(
  DE_results,
  lab = rownames(DE_results),
  x = "avg_log2FC",
  y = "p_val",
  title = "EIF1AX vs. SF3B1",
  subtitle = "CD8_EM T cells",
  caption = "positive foldchange = upregulated in EIF1AX",
  drawConnectors = TRUE,
  pCutoff = 0.05,
  xlim = c(-2.25, 3.5),
  ylim = c(0, 7)
  )

# also run DE between BAP1 and the other two for heatmap visualization and GSEA
DE_results <- FindMarkers(
  pseudo_seu,
  subset.ident = "CD8-EM",
  group.by = "secondary_mutation",
  ident.1 = "BAP1",
  assay = "RNA",
  logfc.threshold = -Inf,
  min.pct = -Inf,
  min.diff.pct = -Inf,
  only.pos = FALSE,
  test.use = "DESeq2"
  )

DE_results$gene <- rownames(DE_results)

BvOthers_markers <- DE_results

# prepare for heatmap
# subset only CD8_EM T cells 
Idents(pseudo_seu) <- "SubTyping"
pseudo_cd8 <- subset(pseudo_seu, idents = "CD8-EM")
# scale data
pseudo_cd8 <- ScaleData(pseudo_cd8)
# extract expression matrix
mat <- pseudo_cd8@assays$RNA$scale.data
# preapare annotation data frame
anno_df <- pseudo_cd8@meta.data %>%
  dplyr::rename(`Secondary Mutation` = secondary_mutation)
# define which genes should be plotted
genes_to_plot <- rownames(dplyr::filter(BvOthers_markers, p_val < 0.05))

# Fig. 4f
pheatmap(
  mat[genes_to_plot,],
  scale = "row", 
  cellwidth = 10,
  show_rownames = FALSE,
  annotation_col = dplyr::select(anno_df, `Secondary Mutation`),
  annotation_colors = list(`Secondary Mutation` = second_mutations_col[1:3]),
  drop_levels = TRUE,
  annotation_names_col = FALSE,
  legend_labels = c("Secondary Mutation"),
  color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(100)
  )

# run GSEA
# prepare ranked gene list
ranks <- BvOthers_markers %>%
  dplyr::filter(!is.na(p_val)) %>%
  dplyr::mutate(metric = -log(p_val)*sign(avg_log2FC)) %>%
  dplyr::arrange(desc(metric)) %>%
  dplyr::pull(metric, name = gene)

# run GSEA using GO terms
GO_results <- gseGO(
  geneList = ranks,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "SYMBOL",
  pvalueCutoff = 0.01,
  seed = TRUE
  )

# Fig. 4g
ggplot(GO_results@result, aes(x = NES, y = fct_reorder(Description, NES), fill = p.adjust)) +
  geom_col() +
  labs(y = NULL) +
  scale_fill_viridis(direction = -1) +
  theme_bw(base_size = 14)

# save metadata for compositional analysis in "Compositions.R" script
qsave(seu@metadata, file = "../pUM_secondarymutation_metadata.qs")

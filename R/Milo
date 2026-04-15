# This script was used to create the following figures:
# Figures 5a-i
# Supplementary Figures 9a-g

library(tidyverse)
library(Seurat)
library(qs)
library(pals)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(patchwork)
library(harmony)
set.seed(1234)

# define colors (see "color_palettes.Rmd")
second_mutations_col <- c(
  "BAP1" = "#1F867BFF", 
  "SF3B1" = "#EAC541FF", 
  "EIF1AX" = "#7B1C5DFF"
)

pUM_samples <- c("pUM01", "pUM02", "pUM03", "pUM04", "pUM05", "pUM06", "pUM09", "pUM10", "pUM11", "pUM12", "pUM14", "pUM15", "pUM16", "pUM17", "pUM19", "pUM20", "pUM21", "pUM22", "pUM23", "pUM24", "pUM25")
mint_colors <- as.character(paletteer::paletteer_c("grDevices::Mint", n = 24))
cats <- unique(pUM_samples)
if (length(cats) > 24) warning("More than 21 categories: colors will recycle")

# Map categories to colors (up to 21)
offset <- 3
col_map <- setNames(mint_colors[ (seq_along(cats) - 1 + offset) %% 24 + 1 ], cats)

# Assign colors to each element of x
primary_samp_col <- col_map[pUM_samples]

# load full data set
seu <- qread("../full_dataset.qs")

# select only melanoma cells
Idents(seu) <- "majority_celltype"
seu <- subset(seu, idents = "Melanoma cells")

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
seu <- FindClusters(seu, verbose = FALSE)

# Fig. 5a
DimPlot(seu, group.by = "secondary_mutation", cols = second_mutations_col)  +
  coord_fixed() +
  labs(title = NULL)

# Suppl. Fig. 9a
DimPlot(seu, group.by = "seurat_clusters", label = TRUE)  +
  coord_fixed()

# convert Seurat object to SingleCellExperiment object because this is needed as input for Milo
SCE <- as.SingleCellExperiment(seu)

# create Milo object
milo_obj <- Milo(SCE)

# generate neighborhoods
milo_obj <- buildGraph(milo_obj, k = 30, d = 30, reduced.dim = "HARMONY")
milo_obj <- makeNhoods(milo_obj, prop = 0.1, k = 30, d = 30, refined = TRUE, reduced_dims = "HARMONY")

# Suppl. Fig. 9b
plotNhoodSizeHist(milo_obj) +
  geom_vline(xintercept = mean(colSums(nhoods(milo_obj))), color = "red", linetype = "dashed") +
  labs(y = "Count")

# count cells in neighborhoods
milo_obj <- countCells(milo_obj, meta.data = data.frame(colData(milo_obj)), samples = "orig.ident")

# prepare design table for differential abundance testing
milo_design <- data.frame(colData(milo_obj))[,c("orig.ident", "secondary_mutation")]
milo_design <- distinct(milo_design)
rownames(milo_design) <- milo_design$orig.ident

# run differential abundance testing 
# perform all pairwise comparisons using contrasts
# use the "graph-overlap" method here for FDR weighting to speed up processing time

contrast.use <- c("secondary_mutationBAP1 - secondary_mutationEIF1AX")
milo_results <- testNhoods(
  milo_obj, 
  design = ~ 0 + secondary_mutation, 
  design.df = milo_design, 
  model.contrasts = contrast.use,
  fdr.weighting = "graph-overlap",
  reduced.dim = "HARMONY"
  )

# Fig. 5b
ggplot(milo_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "Differential abundance of cell neighborhoods\nBAP1 vs. EIF1AX") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(size = 16))

# Fig. 5d
plotNhoodGraphDA(milo_obj, milo_results, alpha = 0.05, is.da = TRUE, highlight.da = TRUE) +
  coord_fixed() +
  theme_classic() +
  scale_x_continuous(breaks = c(-4, 0, 4)) +
  labs(x = "umap_1", y = "umap_2") +
  guides(fill = guide_colorbar(order = 3))

contrast.use <- c("secondary_mutationBAP1 - secondary_mutationSF3B1")
milo_results <- testNhoods(
  milo_obj, 
  design = ~ 0 + secondary_mutation, 
  design.df = milo_design, 
  model.contrasts = contrast.use,
  fdr.weighting = "graph-overlap",
  reduced.dim = "HARMONY"
  )

# Fig. 5c
ggplot(milo_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "Differential abundance of cell neighborhoods\nBAP1 vs. SF3B1") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(size = 16))

# Fig. 5e
plotNhoodGraphDA(milo_obj, milo_results, alpha = 0.05, node_stroke = 0.4) +
  coord_fixed() +
  theme_classic() +
  scale_x_continuous(breaks = c(-4, 0, 4)) +
  labs(x = "umap_1", y = "umap_2") +
  guides(fill = guide_colorbar(order = 3))

contrast.use <- c("secondary_mutationEIF1AX - secondary_mutationSF3B1")
milo_results <- testNhoods(
  milo_obj, 
  design = ~ 0 + secondary_mutation, 
  design.df = milo_design, 
  model.contrasts = contrast.use,
  fdr.weighting = "graph-overlap",
  reduced.dim = "HARMONY"
  )

# Suppl. Fig. 9c
ggplot(milo_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "Differential abundance of cell neighborhoods\nEIF1AX vs. SF3B1") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(size = 16))

# perform differential abundance testing of BAP1 against the other two mutations
contrast.use <- c("secondary_mutationBAP1 - (secondary_mutationEIF1AX + secondary_mutationSF3B1)/2")
milo_results <- testNhoods(
  milo_obj, 
  design = ~ 0 + secondary_mutation, 
  design.df = milo_design, 
  model.contrasts = contrast.use,
  fdr.weighting = "graph-overlap",
  reduced.dim = "HARMONY"
  )

# find NhoodGroups
milo_obj <- buildNhoodGraph(milo_obj)
milo_results <- groupNhoods(
  milo_obj, 
  milo_results, 
  da.fdr = 0.05,
  max.lfc.delta = 3, 
  overlap = 2, 
  compute.new = TRUE
  )

# Fig. 5f
plotNhoodGroups(milo_obj, milo_results) + 
  coord_fixed() +
  plotNhoodGraphDA(milo_obj, milo_results, alpha = 0.05) +
  coord_fixed() +
  plot_layout(guides = "collect")

# Fig. 5g
plotDAbeeswarm(milo_results, group.by = "NhoodGroup", alpha = 0.5)

# add log normalized counts to Milo object
milo_obj <- logNormCounts(milo_obj)

# run DE for all NhoodGroups enriched exclusively in BAP1 on top 2000 variable features
dec <- modelGeneVar(milo_obj)
hvgs <- getTopHVGs(dec, n=2000)

nhood_markers <- findNhoodGroupMarkers(
  milo_obj, 
  milo_results,
  subset.row = hvgs,
  subset.groups = c("1", "10", "13", "14")
)

rownames(nhood_markers) <- nhood_markers$GeneID

# define genes of interest to be highlighted on heatmap
genes_of_interest <- c("MET", "ERBB4", "NRG3", "AKT3", "ZEB2", "S100B", "TYR", "DCT", "CAV1", "TRPM1", "ADD3", "AKAP12", "UACA", "NAV3", "SEMA5A", "SEMA6D")

# define facet panel labels
facet_labels <- as_labeller(
  c(
    "1" = "Nhood Group 1",
    "2" = "Nhood Group 2",
    "9" = "Nhood Group 9",
    "10" = "Nhood Group 10",
    "13" = "Nhood Group 13",
    "14" = "Nhood Group 14"
    )
  )

# define colors
group_colors <- c(
  "2" = "#832424", 
  "9" = "#832424", 
  "1" = "#3A3A98", 
  "10" = "#3A3A98",
  "14" = "#3A3A98",
  "13" = "#3A3A98"
  )

# Create a list of element_rect objects for strip backgrounds
strip_fills <- lapply(group_colors, function(col) element_rect(fill = col))

# prepare expression matrix as input for the heatmap
x <- calcNhoodExpression(milo_obj, subset.row = markers13)
expr_mat <- nhoodExpression(x)[markers13, ]
colnames(expr_mat) <- seq_len(ncol(nhoods(x)))

expr_mat <- expr_mat[,milo_results$NhoodGroup %in% c("13", "1", "10", "14", "2", "9"), drop=FALSE]
expr_mat <- t(apply(expr_mat, 1, function(X) (X - min(X))/(max(X)- min(X))))

rownames(expr_mat) <- sub(pattern = "-", replacement = ".", rownames(expr_mat))

pl_df <- data.frame(t(expr_mat)) %>%
  rownames_to_column("Nhood") %>%
  mutate(Nhood=as.double(Nhood)) %>%
  left_join(milo_results, by="Nhood") %>%
  group_by(NhoodGroup) %>%
  mutate(logFC_rank=rank(logFC, ties.method="random")) %>%
  ungroup()

row.order <- hclust(dist(expr_mat))$order
ordered_features <- rownames(expr_mat)[row.order]

rownames(expr_mat) <- str_replace(rownames(expr_mat), pattern="(^[0-9]+)", replacement="X\\1")

pl_df <- pl_df %>%
  pivot_longer(cols=rownames(expr_mat), names_to='feature', values_to="avg_expr") %>%
  mutate(feature=factor(feature, levels=ordered_features))

pl_df <- pl_df %>%
  mutate(label=ifelse(feature %in% genes_of_interest, as.character(feature), NA)) %>%
  mutate(da = ifelse(NhoodGroup %in% c(2, 9), "depleted", "enriched"))

pl_bottom <- pl_df %>%
  ggplot(aes(logFC_rank, feature, fill=avg_expr)) +
  geom_tile() +
  scale_fill_viridis_c(
    option = "magma", 
    name = "Avg.Expr.", 
    guide = guide_colorbar(
      position = "bottom", 
      label.theme = element_text(angle = 90, vjust = 0.5), title.position = "left", title.vjust = 1)) +
  xlab("Neighbourhoods") + 
  ylab("Genes") +
  scale_x_continuous(expand = c(0.01, 0)) +
  theme_classic(base_size = 16) +
  coord_cartesian(clip = "off") +
  ggh4x::facet_grid2(
    .~factor(NhoodGroup, levels = c(2, 9, 1, 10, 14, 13)), 
    scales = "free", 
    space = "fixed", 
    labeller = facet_labels,
    strip = strip_themed(background_x = strip_fills,
                         text_x = element_text(color = "white"))) +
  theme(
    axis.text.x = element_blank(), 
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.line.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    axis.text.y = element_blank())

xmax <- pl_df %>%
  dplyr::select(Nhood, NhoodGroup) %>%
  dplyr::filter(NhoodGroup == 13) %>%
  unique() %>%
  nrow()

# Fig. 5h
pl_bottom +
  ggrepel::geom_text_repel(data=. %>%
                             filter(!is.na(label)) %>%
                             group_by(label) %>%
                             summarise(logFC_rank=max(logFC_rank), avg_expr=mean(avg_expr), feature=dplyr::first(feature)) %>%
                             dplyr::mutate(NhoodGroup = 13),
                           aes(label=label, x = xmax+1),
                           size = 4,
                           xlim = c(NA, Inf),
                           direction = "y",
                           hjust = -0.25,
                           min.segment.length = 0,
                           point.padding = 0,
                           label.padding = 0,
                           seed =  1234) +
  theme(plot.margin = margin(5, 80, 5, 5, "pt"))

# Suppl. Fig. 9d
ggplot(seu@meta.data, aes(x = seurat_clusters, fill = orig.ident)) +
  geom_bar(position = "fill") +
  labs(x = "Seurat Clusters", y = "Fraction of Cells") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = primary_samp_col, name = "Sample")

# run ORA on NhoodGroup marker genes
# get marker genes
markers1 <- rownames(nhood_markers)[nhood_markers$adj.P.Val_1 < 0.05 & nhood_markers$logFC_1 > 1]
markers10 <- rownames(nhood_markers)[nhood_markers$adj.P.Val_10 < 0.05 & nhood_markers$logFC_10 > 1]
markers13 <- rownames(nhood_markers)[nhood_markers$adj.P.Val_13 < 0.05 & nhood_markers$logFC_13 > 1]
markers14 <- rownames(nhood_markers)[nhood_markers$adj.P.Val_14 < 0.05 & nhood_markers$logFC_14 > 1]

allmarkers <- list(markers1, markers10, markers13, markers14)
names(allmarkers) <- c("NhoodGroup 1", "NhoodGroup 10", "NhoodGroup 13", "NhoodGroup 14")

# Fig. 5i
# Suppl. Fig. 9e-g
for (i in 1:length(allmarkers)){

  GO_enrichment <- enrichGO(
    gene = allmarkers[[i]],
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    readable = TRUE
  )
  
  p1 <- dotplot(GO_enrichment) +
    ggtitle(paste0("GO BP for ", names(allmarkers[i])))
  
  gene_df <- bitr(
    allmarkers[[i]], 
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  Reactome_enrichment <- enrichPathway(
    gene = unique(gene_df$ENTREZID),
    organism = "human", 
    pvalueCutoff = 0.05, 
    readable = TRUE)
  
  p2 <- dotplot(Reactome_enrichment) +
    ggtitle(paste0("Reactome for ", names(allmarkers[i])))
  
  print(p1)
  print(p2)
  
}

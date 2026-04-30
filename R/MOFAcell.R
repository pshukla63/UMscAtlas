
library(Seurat)
library(qs)
library(ggplot2)
library(Matrix)
library(ezRun)
library(RColorBrewer)
library(qs)
library(skimr)
library(edgeR)
library(tidyverse)
library(MOFAcellulaR)
library(MOFA2)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)

# Load integrated Seurat object
scIntegrated <- qread("../Seurat.qs")

DefaultAssay(scIntegrated) <- "RNA"

# Visualise the major cell types
DimPlot(scIntegrated, reduction="umap", label=TRUE, label.size=3, group.by = "majority_celltype")
DimPlot(scIntegrated, reduction="umap", label=TRUE, label.size=3, group.by = "majority_celltype", split.by = "organ", ncol=3)

################################################################################################################################################################

# Explore the metadata
meta <- scIntegrated@meta.data
skim_output <- skim(meta)

# Extract Names of Character/Factor Columns

# Column names from the character summary block:
char_cols <- skim_output %>% 
  filter(skim_type == "character") %>% 
  pull(skim_variable)

# Column names from the factor summary block (if any):
factor_cols <- skim_output %>% 
  filter(skim_type == "factor") %>% 
  pull(skim_variable)

# Combine and ensure your primary sample ID is included
sample_id_col <- "Sample" # Or "study_cohort"
char_factor_cols_to_check <- unique(c(sample_id_col, char_cols))

# Select only the identified character/factor columns
metadata_for_constancy <- meta %>%
  select(all_of(char_factor_cols_to_check))

# Group by the sample ID and calculate the number of unique values (n_distinct) for every other column.
constancy_check <- metadata_for_constancy %>%
  group_by(!!sym(sample_id_col)) %>%
  summarise(across(everything(), n_distinct, .names = "n_distinct_{.col}"))

# Find columns that are constant (n_distinct == 1) across ALL samples
constant_cols <- colnames(constancy_check)[
  sapply(constancy_check, function(x) all(x == 1))
]

# Clean up column names and ensure the sample ID column is present
final_sample_cols <- gsub("n_distinct_", "", constant_cols)

# Ensure 'orig.ident' is at the start and is the only sample identifier used
final_sample_cols <- final_sample_cols[final_sample_cols != sample_id_col]
final_sample_cols <- unique(c(sample_id_col, final_sample_cols))

# Create the final clean sample metadata table
sample_summary_table <- meta %>%
  select(all_of(final_sample_cols)) %>%
  distinct() %>% # Remove redundant cell-level rows
  arrange(!!sym(sample_id_col))

rownames(sample_summary_table) <- NULL

################################################################################################################################################################

# Prepare pseudobulk count matrix
pseudoBulk <- AggregateExpression(scIntegrated, group.by = c("majority_celltype", "Sample"), assays = "RNA", slot = "counts")
pseudoBulk <- pseudoBulk$RNA

colnames(pseudoBulk) <- gsub("-", "_", colnames(pseudoBulk))
dim(pseudoBulk)

# Remove mito and ribo genes
mt.idx <- grep("^MT-", rownames(pseudoBulk))
ribo.idx <- grep("^RP[S|L]", rownames(pseudoBulk))

pseudoBulk <- pseudoBulk[-c(mt.idx, ribo.idx),]
dim(pseudoBulk)

saveRDS(pseudoBulk, "pbRawCounts_perSample_noMitoRibo.rds")

# Prepare information for each donor ~ cell type combination
meta <- data.frame(table(scIntegrated$Sample, scIntegrated$majority_celltype))
colnames(meta) <- c("donor_id", "cell_type", "cell_counts")
rownames(meta) <- paste0(meta$cell_type, "_", meta$donor_id)
rownames(meta) <- gsub("-", "_", rownames(meta))

idx <- which(rownames(meta) %in% colnames(pseudoBulk))
meta <- meta[idx,]
saveRDS(meta, "colData_perSample.rds")

################################################################################################################################################################

# Load pseudobulk count matrix and major cell type information
coldat <- readRDS("colData_perSample.rds")
pb_data <- readRDS("pbRawCounts_perSample_noMitoRibo.rds")
colnames(pb_data) <- gsub("-", "_", colnames(pb_data))

pb_data <- pb_data[,rownames(coldat)]

# Define major cell types
cts <- coldat$cell_type %>% 
  unique() %>%
  set_names()

# Pipeline for differential expression
de_res <- map(cts, function(ct) {
  print(ct)
  ct_meta_data <- coldat %>%
    mutate(test_column = ifelse(cell_type == ct, ct, "rest"))
  
  dat <- DGEList(pb_data, samples = DataFrame(ct_meta_data))
  
  keep <- filterByExpr(dat, group = ct_meta_data$test_column)
  
  dat <- dat[keep,]
  
  dat <- calcNormFactors(dat)
  
  design <- model.matrix(~factor(test_column,
                                 levels = c("rest",ct)), dat$samples)
  
  colnames(design) <- c("int", ct)
  
  dat <- estimateDisp(dat, design)
  
  fit <- glmQLFit(dat, design, robust=TRUE)
  
  res <- glmQLFTest(fit, coef=ncol(design))
  
  de_res <- topTags(res, n = Inf) %>%
    as.data.frame() %>%
    rownames_to_column("gene")
  
  return(de_res)
  
})

de_res <- de_res %>% 
  enframe() %>%
  unnest()

# Keep the significant major cell type markers for background estimation
de_res <- de_res %>%
  dplyr::filter(logFC > 1 & FDR < 0.01) %>%
  arrange(name, FDR, - logFC) 

table(de_res$name)

write_csv(de_res, file = "./bgCellMarkers_perSample.csv")

####################################################################################################################################

# Load metadata
meta <- readRDS("colData_perSample.rds")
rownames(meta) <- gsub("/", "_", rownames(meta))
meta$cell_type <- gsub("/", "_", meta$cell_type)
meta$donor_id <- gsub("-", "_", meta$donor_id)

# Load pseudo bulk count matrix
pseudoBulk <- readRDS("pbRawCounts_perSample_noMitoRibo.rds")
colnames(pseudoBulk) <- gsub("-", "_", colnames(pseudoBulk))
colnames(pseudoBulk) <- gsub("/", "_", colnames(pseudoBulk))

pseudoBulk <- pseudoBulk[,rownames(meta)]

# Load cell type markers based on pseudobulk
mrkr_genes <- read_csv("bgCellMarkers_perSample.csv")
mrkr_genes$name <- gsub("/", "_", mrkr_genes$name)

mrkr_genes <- mrkr_genes %>% #dplyr::filter(!name %in% exclude_ct) %>%
  dplyr::filter(FDR < 0.01, logFC > 1) %>%
  dplyr::select(name, gene) %>%
  dplyr::rename("lineage" = name) %>%
  group_by(lineage) %>%
  nest() %>%
  dplyr::mutate(data = map(data, ~.x[[1]])) %>%
  deframe()

####################################################################################################################################

# Create a SummrizedExp object for pre-processing
pb_obj <- MOFAcellulaR::create_init_exp(counts = pseudoBulk, coldata = meta)

# Filter pseudo bulk samples coming from < 10 cells (Assumption: unreliable gene count estimate)
ct_list <- MOFAcellulaR::filt_profiles(pb_dat = pb_obj,
                                       cts = c("Endothelial cells", "Melanoma cells", "Myeloid cells", "Pericytes", "T_NK cells"),
                                       ncells = 50, 
                                       counts_col = "cell_counts", 
                                       ct_col = "cell_type") 

# Possible to filter cell types with very few samples
ct_list <- MOFAcellulaR::filt_views_bysamples(pb_dat_list = ct_list,
                                              nsamples = 15)

# Identify low and highly expressed genes per cell type
ct_list <- MOFAcellulaR::filt_gex_byexpr(pb_dat_list = ct_list,
                                         min.count = 10, # Modify!!
                                         min.prop = 0.25) # Modify!!

# Check if any sample has very few genes
# ct_list <- MOFAcellulaR::filt_views_bygenes(pb_dat_list = ct_list,
#                                             ngenes = 15)

# Normalize pseudobulk expression profiles
ct_list <- MOFAcellulaR::tmm_trns(pb_dat_list = ct_list,
                                  scale_factor = 1000000)

# Identify highly variable genes per cell type (optional if the number of genes is very low)
ct_list <- MOFAcellulaR::filt_gex_byhvg(pb_dat_list = ct_list,
                                        prior_hvg = NULL,
                                        var.threshold = 0)

# Background gene correction
ct_list <- MOFAcellulaR::filt_gex_bybckgrnd(pb_dat_list = ct_list,
                                            prior_mrks = mrkr_genes)

# Convert into a MOFA object
multiview_dat <- pb_dat2MOFA(pb_dat_list = ct_list, 
                             sample_column = "donor_id")

####################################################################################################################################

# Create MOFA
MOFAobject <- MOFA2::create_mofa(multiview_dat)

# Run MOFA
data_opts <- MOFA2::get_default_data_options(MOFAobject)
train_opts <- MOFA2::get_default_training_options(MOFAobject)
model_opts <- MOFA2::get_default_model_options(MOFAobject)


# This avoids the regularization of multicellular programs per cell type.
# This avoids less sparse gene weights
model_opts$spikeslab_weights <- FALSE 

# Define the number of factors needed
model_opts$num_factors <- 6

# Prepare MOFA model:
MOFAobject <- MOFA2::prepare_mofa(object = MOFAobject,
                                  data_options = data_opts,
                                  model_options = model_opts,
                                  training_options = train_opts)

outfile <- file.path("./mofa_perSample_noMitoRibo.hdf5")
model <- MOFA2::run_mofa(MOFAobject, outfile, use_basilisk = TRUE)

MOFA2::plot_data_overview(model)

####################################################################################################################################

# Load sample metadata
samples <- samples_metadata(model)[[2]]

meta2 <- read_tsv("sample_metadata.tsv")
meta2 <- subset(meta2, sample %in% samples)

samples_metadata(model) <- meta2

####################################################################################################################################

# Visualize sample variability
UMAP_embedding <- MOFAcellulaR::plot_sample_2D(model = model,
                                               method = "UMAP",
                                               metadata = meta,
                                               sample_id_column = "sample",
                                               color_by = "Condition")

manual_colors <- brewer.pal(n = 11, name = "Set3")

# Generate the UMAP coloured by location
ggscatter(
  data = UMAP_embedding, 
  x = "UMAP_1", 
  y = "UMAP_2", 
  color = "location", 
  shape = "origin",
  size = "origin",             # Map size to origin to allow variation
  alpha = 0.5
) + 
  # Update Legend Titles
  labs(color = "Location", shape = "Disease Progression") +
  
  # Color Scale
  scale_color_manual(values = location_colors) +
  
  # Shape Scale
  scale_shape_manual(values = c(
    "healthy" = 16,     
    "primary" = 17,     
    "metastasis" = 18   
  )) +
  
  # SIZE SCALE: Increase Metastasis (18) to ~8 to match the visual weight of the Circle (16)
  scale_size_manual(values = c(
    "healthy" = 6, 
    "primary" = 6.5, 
    "metastasis" = 8.5
  )) +
  
  # THE FIX: Override the legend so all icons look the same size (6)
  guides(
    color = guide_legend(order = 1, override.aes = list(shape = 15, size = 5)), 
    shape = guide_legend(order = 2, override.aes = list(size = 5)),
    size = "none"              # Hide the extra 'size' legend
  ) +
  
  # Theme adjustments
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"), 
    axis.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 10, face = "italic"),
    legend.position = "right",
    legend.box = "vertical"
  )

# Generate the UMAP coloured by secondary mutation
ggscatter(
  data = UMAP_embedding, 
  x = "UMAP_1", 
  y = "UMAP_2", 
  color = "secondary_mutation", 
  shape = "origin",
  size = "origin",             # Map size to origin to allow variation
  alpha = 0.8
) + 
  # Update Legend Titles
  labs(color = "Secondary Mutation", shape = "Disease Progression") +
  geom_hline(yintercept = 0.1, alpha = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.1, alpha = 0.5, linetype = "dashed") +
  # Color Scale
  scale_color_manual(values = second_mutations_col) +
  
  # Shape Scale
  scale_shape_manual(values = c(
    "healthy" = 16,     
    "primary" = 17,     
    "metastasis" = 18   
  )) +
  
  # SIZE SCALE: Increase Metastasis (18) to ~8 to match the visual weight of the Circle (16)
  scale_size_manual(values = c(
    "healthy" = 6, 
    "primary" = 6.5, 
    "metastasis" = 8.5
  )) +
  
  # THE FIX: Override the legend so all icons look the same size (6)
  guides(
    color = guide_legend(order = 1, override.aes = list(shape = 15, size = 5)), 
    shape = guide_legend(order = 2, override.aes = list(size = 5)),
    size = "none"              # Hide the extra 'size' legend
  ) +
  
  # Theme adjustments
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"), 
    axis.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 10, face = "italic"),
    legend.position = "right",
    legend.box = "vertical"
  )

# Make scatterplot between Factor2 and Factor3
ggscatter(factors_sub, x = "Factor2", y = "Factor3", 
          color = "secondary_mutation", 
          shape = "origin",          # Adds shape based on the same variable
          size = 5,                    # Makes dots bigger (default is usually 1.5)
          palette = second_mutations_col,   # Define your custom colors here
          alpha = 0.8) +               # Increased alpha so you can actually see the bigger dots
  geom_hline(yintercept = 0, alpha = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0, alpha = 0.5, linetype = "dashed") +
  scale_shape_manual(values = c("healthy" = 16, "primary" = 17, "metastasis" = 18)) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"), 
        axis.text = element_text(size = 10, face = "bold"),
        strip.text.x = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        legend.title = element_text(size = 10, face = "italic"),
        legend.position = "right")


# Extract explained variance R2
r2_list <- model@cache$variance_explained$r2_per_factor

# If it's a grouped model, prefix the view names with the group names
group <- FALSE 
if(group) {
  r2_list <- purrr::map2(r2_list, names(r2_list), function(dat, nam) {
    colnames(dat) <- paste0(nam, "_", colnames(dat))
    return(dat)
  })
}

r2_matrix <- do.call(cbind, r2_list)
colnames(r2_matrix) <- gsub("T_NK", "T/NK", colnames(r2_matrix))

col_fun_r2 <- colorRamp2(c(0, max(r2_matrix, na.rm = TRUE)), c("white", "orange"))


# Convert list to long format
assoc_df <- assoc_list %>%
  tibble::enframe(name = "test") %>%
  tidyr::unnest(value)

# Matrix for colors (-log10 adj. p-value)
assoc_pvals <- assoc_df %>%
  dplyr::mutate(log_adjpval = -log10(adj_pvalue)) %>%
  dplyr::select(test, Factor, log_adjpval) %>%
  tidyr::pivot_wider(names_from = test, values_from = log_adjpval) %>%
  tibble::column_to_rownames("Factor") %>%
  as.matrix()

# Matrix for significance stars
star_matrix <- assoc_df %>%
  dplyr::mutate(stars = cut(p.value, 
                            breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                            labels = c("***", "**", "*", ""))) %>%
  dplyr::select(test, Factor, stars) %>%
  tidyr::pivot_wider(names_from = test, values_from = stars) %>%
  tibble::column_to_rownames("Factor") %>%
  as.matrix()

# Align rows with the R2 matrix
assoc_pvals <- assoc_pvals[rownames(r2_matrix), , drop = FALSE]
star_matrix <- star_matrix[rownames(r2_matrix), , drop = FALSE]

col_fun_p <- colorRamp2(c(0, max(assoc_pvals, na.rm = TRUE)), c("white", "darkblue"))
# This looks like a 'heat' map and makes high -log10 values (significance) pop in dark red
# col_fun_p <- colorRamp2(
#   seq(0, max(assoc_pvals, na.rm = TRUE), length = 100),
#   hcl.colors(100, "Greens 3", rev = TRUE)
# )

# Heatmap 1: Variance Explained
hb_r2 <- Heatmap(r2_matrix, 
                 name = "R2 (%)", 
                 col = col_fun_r2, 
                 row_names_side = "left",
                 cluster_columns = FALSE, 
                 cluster_rows = FALSE,
                 column_title = "Variance Explained",
                 rect_gp = gpar(col = "white", lwd = 1),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(sprintf("%.1f", r2_matrix[i, j]), x, y, 
                             gp = gpar(fontsize = 9))
                 })

# Heatmap 2: Associations
hb_p <- Heatmap(assoc_pvals, 
                name = "-log10(adj.p)", 
                col = col_fun_p, 
                cluster_columns = FALSE, 
                cluster_rows = FALSE,
                show_row_names = FALSE,
                column_title = "Covariate Associations",
                rect_gp = gpar(col = "white", lwd = 1),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(star_matrix[i, j], x, y, 
                            gp = gpar(fontsize = 14, fontface = "bold"))
                })

# Draw them side-by-side
draw(hb_r2 + hb_p, heatmap_legend_side = "right")

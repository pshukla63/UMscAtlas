
# This script was used to create the following figures:
# Supplementary Figures 8a-d

library(tidyverse)
library(Seurat)
library(qs)
library(vegan)
library(RColorBrewer)
library(scales)
library(janitor)
library(sccomp)
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

tcell_col <- c(
  "CD4_EM" = "#00366CFF", 
  "CD4_EX" = "#004B84FF", 
  "CD4_N" = "#00619FFF", 
  "CD4_REG" = "#0077BCFF", 
  "CD8_EM" = "#408DCEFF", 
  "CD8_EMRA" = "#68A1DCFF", 
  "CD8_EX" = "#88B4E8FF", 
  "CD8_N" = "#A4C7F3FF", 
  "CD8_RM" = "#BED7FBFF", 
  "gdT" = "#D4E6FFFF", 
  "NK_CYTO" = "#E8F2FFFF", 
  "NK_REST" = "#F9F9F9FF", 
  "Proliferating cells" = "#001889FF"
  )

second_mutations_col <- c(
  "BAP1" = "#1F867BFF", 
  "SF3B1" = "#EAC541FF", 
  "EIF1AX" = "#7B1C5DFF"
)

# load the metadata saved in the script "Augur.R"
metadata_raw <- qread("../pUM_secondarymutation_metadata.qs")

# first run compositional analysis for broad cell types
# prepare input
metadata <- metadata_raw %>%
  dplyr::select(orig.ident, secondary_mutation, primary_location) %>%
  rownames_to_column(var = "rownames") %>%
  dplyr::select(-rownames) %>%
  unique()

rownames(metadata) <- NULL

df <- metadata_raw %>%
  dplyr::select(orig.ident, majority_celltype) %>%
  rownames_to_column(var = "rownames") %>%
  dplyr::select(-rownames) %>%
  group_by(majority_celltype, orig.ident) %>%
  summarise(count = n(), .groups = 'drop') %>%
  left_join(metadata, by = "orig.ident")

# run sccomp
sccomp <- df %>%
  dplyr::rename("cell_group" = "majority_celltype") %>%
  sccomp_estimate( 
    formula_composition = ~ secondary_mutation + primary_location, 
    .sample = "orig.ident",
    .cell_group = cell_group,
    .count = count,
    cores = 30, verbose = FALSE, 
    variational_inference = FALSE)

sccomp_adjusted <- sccomp %>%
  sccomp_test() %>%
  sccomp_remove_unwanted_variation(formula_composition = ~ secondary_mutation + primary_location) %>%
  select(orig.ident, cell_group, adjusted_counts) %>%
  merge(metadata)

# Suppl. Fig. 8a
sccomp_adjusted %>%
  ggplot(aes(x = orig.ident, y = sqrt(adjusted_counts), fill = cell_group)) + 
  geom_col(position = "fill") + 
  scale_fill_manual(values = majority_celltype_col, name = "Cell Type") + 
  facet_wrap(~ secondary_mutation, scales = "free_x", space = "free_x") + 
  labs(y = expression(sqrt(Fraction)), x = NULL) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# prepare input for NMDS
df_comp <- sccomp_adjusted %>%
  pivot_wider(id_cols = orig.ident, names_from = cell_group, values_from = adjusted_counts) %>%
  mutate(across(everything(), ~ replace(., is.na(.), 0))) %>%
  inner_join(metadata, by = "orig.ident")

# calculate distance matrix
md_raw <- df_comp %>%
  column_to_rownames("orig.ident") %>%
  dplyr::select(-primary_location, -secondary_mutation) %>%
  metaMDS(., distance = "bray", autotransform = F)

# get points for plotting and add metadata
raw_data <- md_raw$points %>%
  as.data.frame() %>%
  rownames_to_column("orig.ident") %>%
  left_join(df_comp, by = "orig.ident") %>%
  group_by(secondary_mutation) %>%
  mutate(cent.1 = mean(MDS1), 
         cent.2 = mean(MDS2))

# prepare input for PERMANOVA
mat <- df_comp %>% 
  dplyr::select(-orig.ident, -primary_location, -secondary_mutation)

# run PERMANOVA
ad <- adonis2(mat ~ secondary_mutation , data = df_comp , method = "bray", sqrt.dist = T)

# get p-value for printing onto figure
p <- data.frame(p = paste0("PERMANOVA p = ", ad$`Pr(>F)`)[1])

# Suppl. Fig. 8b
ggplot(raw_data, aes(x = MDS1, y = MDS2, colour = secondary_mutation)) + 
  geom_point(size = 2) + 
  scale_colour_manual(values = second_mutations_col, name = "Secondary Mutation") + 
  geom_segment(aes(xend = cent.1, yend = cent.2)) +
  geom_text(data = p, x = -Inf, y = Inf, label = p, show.legend = F, 
            hjust = 0, vjust = 1, color = "black") + 
  ggtitle("Broad Cell Types") +
  stat_ellipse(level = 0.8) + 
  theme_bw(base_size = 14)

# do the same thing for the fine-grained T/NK cell subtypes
# define subtypes to make subset
TNK_subtypes <- c("CD4_EM", 
                  "NK_REST", 
                  "CD8_EM", 
                  "CD4_N", 
                  "CD4_REG", 
                  "CD4_EX", 
                  "Proliferating cells", 
                  "CD8_N", 
                  "CD8_RM", 
                  "gdT", 
                  "MAIT", 
                  "NK_CYTO", 
                  "CD8_EMRA", 
                  "CD8_EX"
                  )

metadata_TNK_raw <- metadata_raw %>%
  dplyr::filter(detailed_celltype %in% TNK_subtypes)

# prepare input
metadata_TNK <- metadata_TNK_raw %>%
  dplyr::select(orig.ident, secondary_mutation, primary_location) %>%
  rownames_to_column(var = "rownames") %>%
  dplyr::select(-rownames) %>%
  unique()

rownames(metadata_TNK) <- NULL

df_TNK <- metadata_TNK_raw %>%
  dplyr::select(orig.ident, detailed_celltype) %>%
  rownames_to_column(var = "rownames") %>%
  dplyr::select(-rownames) %>%
  group_by(detailed_celltype, orig.ident) %>%
  summarise(count = n(), .groups = 'drop') %>%
  left_join(metadata_TNK, by = "orig.ident")

# run sccomp
sccomp_TNK <- df_TNK %>%
  dplyr::rename("cell_group" = "detailed_celltype") %>%
  sccomp_estimate( 
    formula_composition = ~ secondary_mutation + primary_location, 
    .sample = "orig.ident",
    .cell_group = cell_group,
    .count = count,
    cores = 30, verbose = FALSE, 
    variational_inference = FALSE)

sccomp_adjusted_TNK <- sccomp_TNK %>%
  sccomp_test() %>%
  sccomp_remove_unwanted_variation(formula_composition = ~ secondary_mutation + primary_location) %>%
  select(orig.ident, cell_group, adjusted_counts) %>%
  merge(metadata_TNK)

# Suppl. Fig. 8c
sccomp_adjusted_TNK %>%
  ggplot(aes(x = orig.ident, y = sqrt(adjusted_counts), fill = cell_group)) + 
  geom_col(position = "fill") + 
  scale_fill_manual(values = tcell_col, name = "Cell Type") + 
  facet_wrap(~ secondary_mutation, scales = "free_x", space = "free_x") + 
  labs(y = expression(sqrt(Fraction)), x = NULL) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# prepare input for NMDS
df_comp_TNK <- sccomp_adjusted_TNK %>%
  pivot_wider(id_cols = orig.ident, names_from = cell_group, values_from = adjusted_counts) %>%
  mutate(across(everything(), ~ replace(., is.na(.), 0))) %>%
  inner_join(metadata_TNK, by = "orig.ident")

# calculate distance matrix
md_raw_TNK <- df_comp_TNK %>%
  column_to_rownames("orig.ident") %>%
  dplyr::select(-primary_location, -secondary_mutation) %>%
  metaMDS(., distance = "bray", autotransform = F)

# get points for plotting and add metadata
raw_data_TNK <- md_raw_TNK$points %>%
  as.data.frame() %>%
  rownames_to_column("orig.ident") %>%
  left_join(df_comp_TNK, by = "orig.ident") %>%
  group_by(secondary_mutation) %>%
  mutate(cent.1 = mean(MDS1), 
         cent.2 = mean(MDS2))

# prepare input for PERMANOVA
mat_TNK <- df_comp_TNK %>% 
  dplyr::select(-orig.ident, -primary_location, -secondary_mutation)

# run PERMANOVA
ad_TNK <- adonis2(mat_TNK ~ secondary_mutation , data = df_comp_TNK , method = "bray", sqrt.dist = T)

# get p-value for printing onto figure
p_TNK <- data.frame(p = paste0("PERMANOVA p = ", ad_TNK$`Pr(>F)`)[1])

# Suppl. Fig. 8d
ggplot(raw_data_TNK, aes(x = MDS1, y = MDS2, colour = secondary_mutation)) + 
  geom_point(size = 2) + 
  scale_colour_manual(values = second_mutations_col, name = "Secondary Mutation") + 
  geom_segment(aes(xend = cent.1, yend = cent.2)) +
  geom_text(data = p_TNK, x = -Inf, y = Inf, label = p, show.legend = F, 
            hjust = 0, vjust = 1, color = "black") + 
  ggtitle("T/NK Cell Subtypes") +
  stat_ellipse(level = 0.8) + 
  theme_bw(base_size = 14)

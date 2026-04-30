# This script generates the revised CytoTRACE2 + metaprogram plots (Groups B/C/D).
# Outputs: B1-B7 CytoTRACE2 panels, C1-C3 metaprogram violins, D1-D14 figures (Figures 2-4 and Extended Data).

# ============================================================================
# p31662 UM Atlas - Revised CytoTRACE2 + Metaprogram Plots (Groups B, C, D)
# Updated: 2026-03-27 (Round 2)
# Round 2 fixes: newest dataset, CytoTRACE2 barcode transfer, B1 location_colors,
#   B4 theme_classic, B5 healthy samples + ordering, B7 per-sample PNG,
#   D1-D4 from newest object, MP summary heatmap
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(qs)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ggpubr)
  library(ggbeeswarm)
  library(cowplot)
  library(patchwork)
  library(pals)
  library(viridis)
  library(ComplexHeatmap)
  library(circlize)
})

# Source color palettes (colors.R lives alongside this script in R/)
source("colors.R")

# Paths
base_dir <- ".."
spatial_csv <- "../unified_spatial_data_boston.csv"

out_b <- file.path(base_dir, "Revised_Plots_20260325", "GroupB_CytoTRACE")
out_c <- file.path(base_dir, "Revised_Plots_20260325", "GroupC_Metaprogram")
out_d <- file.path(base_dir, "Revised_Plots_20260325", "GroupD_Unchanged")
dir.create(out_b, recursive = TRUE, showWarnings = FALSE)
dir.create(out_c, recursive = TRUE, showWarnings = FALSE)
dir.create(out_d, recursive = TRUE, showWarnings = FALSE)

cat("=== p31662 CytoTRACE2 + Metaprogram Plot Revisions (Round 2) ===\n\n")

# Constants
origin_order <- c("healthy", "primary", "metastasis")
potency_labels <- c("Differentiated", "Unipotent", "Oligopotent", "Multipotent", "Pluripotent", "Totipotent")
potency_colors <- setNames(c("#9E0142", "#F46D43", "#FEE08B", "#E6F598", "#66C2A5", "#5E4FA2"), potency_labels)
pairwise_origin <- list(c("healthy", "primary"), c("primary", "metastasis"), c("healthy", "metastasis"))

# ============================================================================
# PHASE 1: Load newest dataset + transfer CytoTRACE2 scores
# ============================================================================
cat("Loading newest dataset (13_final_Dataset_reclustered_28012025.qs)...\n")
um_atlas <- qs::qread(file.path(base_dir, "13_final_Dataset_reclustered_28012025.qs"), nthreads = 16)
cat("  Loaded:", ncol(um_atlas), "cells,", nrow(um_atlas), "genes\n")
cat("  Samplenames:", length(unique(um_atlas$samplename)), "\n")

# Transfer CytoTRACE2 scores from old object
cat("Loading old CytoTRACE2 object for score transfer...\n")
old_ct2 <- qs::qread(file.path(base_dir, "UM_atlas_cytotrace2.qs2"), nthreads = 16)
cat("  Old object:", ncol(old_ct2), "cells\n")

# Transfer by barcode matching
common_cells <- intersect(colnames(um_atlas), colnames(old_ct2))
cat("  Common barcodes:", length(common_cells), "/", ncol(um_atlas), "\n")

um_atlas$CytoTRACE2_Score <- NA_real_
um_atlas$CytoTRACE2_Potency <- NA_character_
um_atlas$CytoTRACE2_Score[match(common_cells, colnames(um_atlas))] <-
  old_ct2$CytoTRACE2_Score[match(common_cells, colnames(old_ct2))]
um_atlas$CytoTRACE2_Potency[match(common_cells, colnames(um_atlas))] <-
  old_ct2$CytoTRACE2_Potency[match(common_cells, colnames(old_ct2))]

n_scored <- sum(!is.na(um_atlas$CytoTRACE2_Score))
cat("  Transferred CytoTRACE2 scores to", n_scored, "cells (",
    round(100 * n_scored / ncol(um_atlas), 1), "%)\n")

rm(old_ct2)
gc()

# ============================================================================
# PHASE 2: CytoTRACE2 plots (B1-B3, B5, D1-D4)
# ============================================================================

# Compute patient-level stats for boxplots
malignant_subset <- subset(um_atlas, Malignant == "Malignant" | majority_celltype == "Melanocytes")
cat("  Melanocytic subset:", ncol(malignant_subset), "cells\n")

sample_stats <- malignant_subset@meta.data |>
  group_by(Sample, Patient_nr) |>
  summarise(
    mean_score = mean(CytoTRACE2_Score, na.rm = TRUE),
    n_cells = n(),
    origin = first(origin),
    organ = first(organ),
    secondary_mutation = first(secondary_mutation),
    primary_mutation = first(primary_mutation),
    Gender = first(Gender),
    location = first(LocationPrimaryTumor),
    Age_numeric = as.numeric(first(Age_atDiagnosisY)),
    samplename = first(samplename),
    .groups = "drop"
  ) |>
  filter(n_cells >= 10) |>
  mutate(
    origin = factor(origin, levels = origin_order),
    age_group_base = case_when(
      is.na(Age_numeric) ~ NA_character_,
      Age_numeric < 40 ~ "Young",
      Age_numeric < 60 ~ "Middle",
      TRUE ~ "Old"
    )
  )

# Compute actual age ranges
age_ranges <- sample_stats |>
  filter(!is.na(age_group_base)) |>
  group_by(age_group_base) |>
  summarise(min_age = floor(min(Age_numeric)), max_age = floor(max(Age_numeric)), .groups = "drop")

age_label_map <- setNames(
  paste0(age_ranges$age_group_base, " [", age_ranges$min_age, "-", age_ranges$max_age, "]"),
  age_ranges$age_group_base
)
cat("Age group labels:", paste(age_label_map, collapse = ", "), "\n")

sample_stats <- sample_stats |>
  mutate(
    age_group = ifelse(is.na(age_group_base), NA_character_, age_label_map[age_group_base]),
    age_group = factor(age_group, levels = age_label_map[c("Young", "Middle", "Old")])
  )

# --- B1: Patient-Level CytoTRACE Scores by Disease Stage (colored by Location) ---
cat("B1: Patient-Level CytoTRACE Scores by Disease Stage (location colors)...\n")

p_b1 <- ggplot(sample_stats, aes(x = origin, y = mean_score)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, fill = "gray80") +
  geom_quasirandom(aes(color = location, size = n_cells), alpha = 0.7, dodge.width = 0.8) +
  stat_compare_means(comparisons = pairwise_origin, method = "wilcox.test",
                     label = "p.signif", tip.length = 0.02) +
  scale_x_discrete(limits = origin_order) +
  scale_color_manual(values = location_colors, name = "Uveal Location", na.value = "#CCCCCC") +
  scale_size_continuous(name = "Cell count", range = c(2, 8),
                        breaks = c(100, 500, 1000, 5000, 10000)) +
  theme_cowplot() +
  ggtitle("Patient-Level CytoTRACE2 Scores by Disease Stage") +
  xlab("Disease Progression") +
  ylab("Mean CytoTRACE2 Score per Patient") +
  theme(legend.position = "right", plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

ggsave(file.path(out_b, "B1_patient_scores_disease_stage.pdf"), p_b1,
       device = cairo_pdf, width = 12, height = 8)
cat("  Saved B1\n")

# --- B2: CytoTRACE Scores by Secondary Mutation (Malignant cells) ---
cat("B2: CytoTRACE Scores by Secondary Mutation...\n")

mal_only <- subset(um_atlas, Malignant == "Malignant")
mut_data <- mal_only@meta.data |>
  filter(secondary_mutation %in% c("BAP1", "SF3B1", "EIF1AX")) |>
  mutate(secondary_mutation = factor(secondary_mutation, levels = c("BAP1", "SF3B1", "EIF1AX")))

mut_comparisons <- list(c("BAP1", "SF3B1"), c("SF3B1", "EIF1AX"), c("BAP1", "EIF1AX"))

p_b2 <- ggviolin(mut_data, x = "secondary_mutation", y = "CytoTRACE2_Score",
                  fill = "secondary_mutation",
                  palette = unname(second_mutations_col[c("BAP1", "SF3B1", "EIF1AX")]),
                  add = "boxplot", add.params = list(fill = "white", width = 0.1)) +
  stat_compare_means(comparisons = mut_comparisons, method = "wilcox.test",
                     label = "p.signif", tip.length = 0.02) +
  labs(title = "CytoTRACE2 Scores by Secondary Mutation (Malignant Cells)",
       x = "Secondary Mutation", y = "CytoTRACE2 Score") +
  theme(legend.position = "none")

ggsave(file.path(out_b, "B2_scores_secondary_mutation.pdf"), p_b2,
       device = cairo_pdf, width = 8, height = 7)
cat("  Saved B2\n")
rm(mal_only, mut_data)

# --- B3: Patient-Level CytoTRACE Scores by Metastasis Organ ---
cat("B3: Patient-Level CytoTRACE Scores by Metastasis Organ...\n")

met_stats <- sample_stats |> filter(origin == "metastasis", !is.na(organ))
met_organs <- sort(unique(met_stats$organ))
met_organ_pairs <- if (length(met_organs) >= 2) combn(met_organs, 2, simplify = FALSE) else list()

p_b3 <- ggplot(met_stats, aes(x = organ, y = mean_score, fill = organ)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_quasirandom(aes(color = organ, size = n_cells), alpha = 0.7, dodge.width = 0.8) +
  {if (length(met_organ_pairs) > 0) stat_compare_means(comparisons = met_organ_pairs,
                                                        method = "wilcox.test",
                                                        label = "p.signif", tip.length = 0.02)} +
  scale_fill_manual(values = organ_colors, name = "Organ") +
  scale_color_manual(values = organ_colors, name = "Organ") +
  scale_size_continuous(name = "Cell count", range = c(2, 8)) +
  theme_cowplot() +
  ggtitle("Patient-Level CytoTRACE2 Scores by Metastasis Organ") +
  xlab("Organ") + ylab("Mean CytoTRACE2 Score per Patient") +
  theme(legend.position = "right", plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

ggsave(file.path(out_b, "B3_patient_scores_metastasis_organ.pdf"), p_b3,
       device = cairo_pdf, width = 10, height = 8)
cat("  Saved B3\n")

# --- B5: Patient-Level CytoTRACE2 by Disease Stage (colored by Age Class) ---
cat("B5: Patient-Level CytoTRACE2 by Disease Stage (Age Class)...\n")

# Include ALL samples (healthy has age data too)
age_colors_relabeled <- setNames(
  unname(ageGroup_col[age_ranges$age_group_base]),
  age_label_map[age_ranges$age_group_base]
)

p_b5 <- ggplot(sample_stats |> filter(!is.na(age_group)),
               aes(x = origin, y = mean_score)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, fill = "gray80") +
  geom_quasirandom(aes(color = age_group, size = n_cells), alpha = 0.7, dodge.width = 0.8) +
  scale_x_discrete(limits = origin_order) +
  scale_color_manual(values = age_colors_relabeled, name = "Age Class", na.translate = FALSE,
                     limits = age_label_map[c("Young", "Middle", "Old")]) +
  scale_size_continuous(name = "Cell count", range = c(2, 8),
                        breaks = c(100, 500, 1000, 5000, 10000)) +
  theme_cowplot() +
  ggtitle("Patient-Level CytoTRACE2 Scores by Disease Stage",
          subtitle = "Colored by Age Class") +
  xlab("Disease Progression") + ylab("Mean CytoTRACE2 Score per Patient") +
  theme(legend.position = "right", plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

ggsave(file.path(out_b, "B5_patient_scores_age_class.pdf"), p_b5,
       device = cairo_pdf, width = 12, height = 8)
cat("  Saved B5\n")

# --- D1: CytoTRACE2 Potency Categories UMAP (from newest object) ---
cat("D1: CytoTRACE2 Potency Categories UMAP (newest dataset)...\n")

um_atlas$CytoTRACE2_Potency <- factor(um_atlas$CytoTRACE2_Potency, levels = potency_labels)
p_d1 <- DimPlot(um_atlas, group.by = "CytoTRACE2_Potency", cols = potency_colors,
                raster = FALSE, order = TRUE) +
  theme(aspect.ratio = 1) +
  ggtitle("CytoTRACE2 Potency Categories")

ggsave(file.path(out_d, "D01_potency_categories_umap.png"), p_d1,
       width = 10, height = 9, dpi = 300, bg = "white")
cat("  Saved D1\n")

# --- D2: Developmental Potential by Cell Type violins ---
cat("D2: Developmental Potential by Cell Type violins...\n")

p_d2 <- VlnPlot(um_atlas, features = "CytoTRACE2_Score", group.by = "majority_celltype", pt.size = 0) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black") +
  scale_fill_manual(values = majority_celltype_col) +
  ggtitle("Developmental Potential by Cell Type") +
  xlab("Cell Type") + ylab("CytoTRACE2 Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

ggsave(file.path(out_d, "D02_developmental_potential_violins.pdf"), p_d2,
       device = cairo_pdf, width = 12, height = 7)
cat("  Saved D2\n")

# --- D3: CytoTRACE2 Scores by pUM Location violins ---
cat("D3: CytoTRACE2 Scores by pUM Location...\n")

pum_mal <- subset(malignant_subset, origin == "primary" &
                    LocationPrimaryTumor %in% names(primary_location_colors))

p_d3 <- VlnPlot(pum_mal, features = "CytoTRACE2_Score",
                 group.by = "LocationPrimaryTumor", pt.size = 0) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black") +
  scale_fill_manual(values = primary_location_colors) +
  ggtitle("CytoTRACE2 Scores by Anatomical Location in pUM Malignant Cells") +
  xlab("Anatomical Location") + ylab("CytoTRACE2 Score") +
  theme(legend.position = "none")

ggsave(file.path(out_d, "D03_cytotrace2_pum_location_violins.pdf"), p_d3,
       device = cairo_pdf, width = 9, height = 7)
cat("  Saved D3\n")
rm(pum_mal)

# --- D4: Patient-Level by Disease Stage (Gender) ---
cat("D4: Patient-Level by Disease Stage (Gender)...\n")

p_d4 <- ggplot(sample_stats, aes(x = origin, y = mean_score)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, fill = "gray80") +
  geom_quasirandom(aes(color = Gender, size = n_cells), alpha = 0.7, dodge.width = 0.8) +
  scale_x_discrete(limits = origin_order) +
  scale_color_manual(values = gender_col, name = "Gender") +
  scale_size_continuous(name = "Cell count", range = c(2, 8)) +
  theme_cowplot() +
  ggtitle("Patient-Level CytoTRACE2 Scores by Disease Stage",
          subtitle = "Colored by Gender") +
  xlab("Disease Progression") + ylab("Mean CytoTRACE2 Score per Patient") +
  theme(legend.position = "right", plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

ggsave(file.path(out_d, "D04_patient_scores_gender.pdf"), p_d4,
       device = cairo_pdf, width = 10, height = 8)
cat("  Saved D4\n")

rm(um_atlas, malignant_subset, sample_stats)
gc()

# ============================================================================
# PHASE 3: Spatial CSV (B6, B7)
# ============================================================================
cat("\nLoading spatial CSV for B6/B7...\n")
spatial_df <- read_csv(spatial_csv, show_col_types = FALSE)
cat("  Loaded:", nrow(spatial_df), "cells\n")
cat("  CytoTRACE2_Score range:", range(spatial_df$CytoTRACE2_Score, na.rm = TRUE), "\n")
cat("  NOTE: CytoTRACE2_Score = Top - Bottom UCell signature difference (can be negative)\n")

# --- B6: Spatial CytoTRACE Scores ---
cat("B6: Spatial CytoTRACE Scores...\n")

p_b6 <- ggplot(spatial_df, aes(x = x, y = y, color = CytoTRACE2_Score)) +
  geom_point(size = 0.1, alpha = 0.5) +
  scale_color_viridis_c(option = "magma", name = "Stemness\nSignature\nScore") +
  facet_wrap(~Sample, scales = "free") +
  xlab(expression("X (" * mu * "m)")) + ylab(expression("Y (" * mu * "m)")) +
  ggtitle("Spatial CytoTRACE2 Stemness Signature Scores") +
  theme_minimal() +
  theme(aspect.ratio = 1, strip.text = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8))

ggsave(file.path(out_b, "B6_spatial_cytotrace2.pdf"), p_b6,
       device = cairo_pdf, width = 16, height = 12)
cat("  Saved B6\n")

# --- B7: Correlation CNV Score - Stemness Score (per-sample, PNG) ---
cat("B7: Correlation CNV - Stemness (per-sample)...\n")

p_b7 <- ggplot(spatial_df, aes(x = CytoTRACE2_Score, y = cnv_score)) +
  geom_point(color = "grey50", alpha = 0.2, size = 0.5) +
  geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 0.8) +
  stat_cor(method = "pearson", label.x.npc = 0.05, label.y.npc = 0.95, size = 3.5) +
  facet_wrap(~Sample, scales = "free", ncol = 3) +
  xlab("Stemness Signature Score (Top - Bottom)") + ylab("CNV Score") +
  ggtitle("Correlation: CNV Score vs Stemness Signature Score (per Sample)") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))

ggsave(file.path(out_b, "B7_cnv_stemness_correlation.png"), p_b7,
       width = 14, height = 10, dpi = 300, bg = "white")
cat("  Saved B7 (PNG)\n")

rm(spatial_df)
gc()

# ============================================================================
# PHASE 4: Metaprogram object (B4, C1-C3, D5-D13, MP heatmap)
# ============================================================================
cat("\nLoading melanocytic_metaprograms.qs2...\n")
mp_obj <- tryCatch(
  qs2::qs_read(file.path(base_dir, "melanocytic_metaprograms.qs2"), nthreads = 16),
  error = function(e) {
    cat("  qs2 failed, trying qs::qread...\n")
    qs::qread(file.path(base_dir, "melanocytic_metaprograms.qs2"), nthreads = 16)
  }
)
cat("  Loaded:", ncol(mp_obj), "cells\n")

# Discover metaprogram columns
all_cols <- colnames(mp_obj@meta.data)
mp_cols <- grep("^MP\\d+ ", all_cols, value = TRUE)
cat("  Found", length(mp_cols), "metaprogram columns\n")

mp45_col <- grep("^MP45 ", mp_cols, value = TRUE)[1]
mp17_col <- grep("^MP17 ", mp_cols, value = TRUE)[1]
mp9_col <- grep("^MP9 ", mp_cols, value = TRUE)[1]
mp64_col <- grep("^MP64 ", mp_cols, value = TRUE)[1]
mp4_col <- grep("^MP4 Chromatin", mp_cols, value = TRUE)[1]
mp41_col <- grep("^MP41 ", mp_cols, value = TRUE)[1]
mp66_col <- grep("^MP66 ", mp_cols, value = TRUE)[1]

cat("  MP45:", mp45_col, "\n  MP17:", mp17_col, "\n  MP9:", mp9_col, "\n")

mp_obj$origin <- factor(mp_obj$origin, levels = origin_order)

# --- B4: MP45 Skin-pigmentation by intraocular location (theme_classic) ---
cat("B4: MP45 Skin-pigmentation by intraocular location...\n")

pum_mal_mp <- subset(mp_obj, origin == "primary" & Malignant == "Malignant")
loc_data <- pum_mal_mp@meta.data |>
  filter(LocationPrimaryTumor %in% c("choroid", "ciliary body", "iris")) |>
  mutate(Location = factor(LocationPrimaryTumor, levels = c("iris", "choroid", "ciliary body")))

loc_comparisons <- list(c("iris", "choroid"), c("choroid", "ciliary body"), c("iris", "ciliary body"))

p_b4 <- ggplot(loc_data, aes(x = Location, y = .data[[mp45_col]], fill = Location)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  stat_compare_means(comparisons = loc_comparisons, method = "wilcox.test",
                     label = "p.signif", tip.length = 0.02) +
  scale_fill_manual(values = primary_location_colors) +
  labs(title = "MP45: Skin-pigmentation by Intraocular Location",
       x = "Location", y = "UCell Score") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

ggsave(file.path(out_b, "B4_MP45_skin_pigmentation_location.pdf"), p_b4,
       device = cairo_pdf, width = 8, height = 7)
cat("  Saved B4\n")
rm(pum_mal_mp, loc_data)

# --- C1-C3: Metaprogram violins with axis labels ---
cat("C1-C3: Metaprogram violins...\n")

mp_violin <- function(obj, mp_col, mp_name, filename) {
  plot_data <- obj@meta.data |>
    dplyr::select(origin, score = all_of(mp_col)) |>
    filter(!is.na(origin), !is.na(score)) |>
    mutate(origin = factor(origin, levels = origin_order))

  p <- ggviolin(plot_data, x = "origin", y = "score", fill = "origin",
                palette = unname(origin_colors[origin_order]),
                add = "boxplot", add.params = list(fill = "white", width = 0.1)) +
    stat_compare_means(comparisons = pairwise_origin, method = "wilcox.test",
                       label = "p.signif", hide.ns = FALSE, tip.length = 0.02) +
    stat_compare_means(method = "kruskal.test", label.y.npc = 0.95, label = "p.format") +
    scale_x_discrete(limits = origin_order) +
    labs(title = mp_name, x = "Origin", y = "UCell Score") +
    theme_classic(base_size = 14) +
    theme(legend.position = "none", plot.title = element_text(face = "bold", size = 14))

  ggsave(filename, p, device = cairo_pdf, width = 7, height = 7)
}

mp_violin(mp_obj, mp45_col, "MP45: Skin-pigmentation", file.path(out_c, "C1_MP45_skin_pigmentation.pdf"))
cat("  Saved C1\n")
mp_violin(mp_obj, mp17_col, "MP17: EMT-III", file.path(out_c, "C2_MP17_EMT_III.pdf"))
cat("  Saved C2\n")
mp_violin(mp_obj, mp9_col, "MP9: Stress 2", file.path(out_c, "C3_MP9_Stress2.pdf"))
cat("  Saved C3\n")

# --- D5-D11: MP UMAP plots (as PNG) ---
cat("D5-D11: MP UMAP plots...\n")

mp_umap <- function(obj, mp_col, mp_name, filename) {
  p <- FeaturePlot(obj, features = mp_col, reduction = "umap", raster = FALSE) +
    theme(aspect.ratio = 1) +
    scale_color_viridis_c(option = "magma") +
    ggtitle(mp_name)
  ggsave(filename, p, width = 8, height = 7, dpi = 300, bg = "white")
}

mp_umap(mp_obj, mp45_col, "MP45: Skin-pigmentation", file.path(out_d, "D05_MP45_umap.png"))
cat("  Saved D5\n")
mp_umap(mp_obj, mp17_col, "MP17: EMT-III", file.path(out_d, "D06_MP17_umap.png"))
cat("  Saved D6\n")
mp_umap(mp_obj, mp9_col, "MP9: Stress 2", file.path(out_d, "D07_MP9_umap.png"))
cat("  Saved D7\n")
mp_umap(mp_obj, mp64_col, "MP64: Adherens", file.path(out_d, "D08_MP64_umap.png"))
cat("  Saved D8\n")
mp_umap(mp_obj, mp4_col, "MP4: Chromatin", file.path(out_d, "D09_MP4_umap.png"))
cat("  Saved D9\n")
mp_umap(mp_obj, mp41_col, "MP41: NPC/OPC", file.path(out_d, "D10_MP41_umap.png"))
cat("  Saved D10\n")
mp_umap(mp_obj, mp66_col, "MP66: Unassigned 1", file.path(out_d, "D11_MP66_umap.png"))
cat("  Saved D11\n")

# --- D12: Interferon/MHC by Disease Origin ---
cat("D12: Interferon/MHC by Disease Origin...\n")

ifn_mps <- grep("^MP2[2-4] ", mp_cols, value = TRUE)

p_d12 <- DotPlot(mp_obj, features = ifn_mps, group.by = "origin") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Interferon/MHC Metaprograms by Disease Origin")

ggsave(file.path(out_d, "D12_interferon_mhc_origin.pdf"), p_d12,
       device = cairo_pdf, width = 8, height = 6)
cat("  Saved D12\n")

# --- D13: Interferon/MHC by Treatment (Metastases only) ---
cat("D13: Interferon/MHC by Treatment (Metastases)...\n")

met_cells <- subset(mp_obj, origin == "metastasis")

p_d13 <- DotPlot(met_cells, features = ifn_mps,
                 group.by = "TreatmentStage_ofProcessedSample") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Interferon/MHC Metaprograms by Treatment (Metastases Only)")

ggsave(file.path(out_d, "D13_interferon_mhc_treatment.pdf"), p_d13,
       device = cairo_pdf, width = 10, height = 6)
cat("  Saved D13\n")
rm(met_cells)

# --- NEW: Metaprogram Summary Heatmap (top 30 correlation) ---
cat("MP Summary Heatmap: Top 30 metaprogram correlations...\n")

# Get mean scores per metaprogram
mp_means <- colMeans(mp_obj@meta.data[, mp_cols], na.rm = TRUE)
top30_mps <- names(sort(mp_means, decreasing = TRUE))[1:30]

cor_mat <- cor(mp_obj@meta.data[, top30_mps], use = "complete.obs")
rownames(cor_mat) <- gsub("^(MP[0-9]+) .*", "\\1", rownames(cor_mat))
colnames(cor_mat) <- gsub("^(MP[0-9]+) .*", "\\1", colnames(cor_mat))

col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
ht_mp <- Heatmap(
  cor_mat,
  name = "Pearson\nCorrelation",
  col = col_fun,
  cluster_rows = TRUE, cluster_columns = TRUE,
  show_row_names = TRUE, show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_title = "Top 30 Metaprogram Correlations (All Melanocytic Cells)",
  heatmap_legend_param = list(direction = "vertical")
)

cairo_pdf(file.path(out_d, "D14_metaprogram_summary_heatmap.pdf"), width = 14, height = 13)
draw(ht_mp)
dev.off()
cat("  Saved D14 (MP summary heatmap)\n")

rm(mp_obj)
gc()

# ============================================================================
cat("\n=== Groups B, C, D complete! ===\n")
for (d in c(out_b, out_c, out_d)) {
  cat("\n", d, ":\n")
  cat(paste(" ", list.files(d), collapse = "\n"), "\n")
}

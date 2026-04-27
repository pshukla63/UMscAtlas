# This script generates the revised GloScope plots (Group A) for the manuscript revision round.
# Outputs: A1 ANOSIM heatmap, A2-A6 MDS plots, A7 divergence heatmap (Extended Data).

# ============================================================================
# p31662 UM Atlas - Revised GloScope Plots (Group A) - Round 2
# Updated: 2026-03-27
# Tasks A1-A7 from Victoria/Mauro/Aizhan request list
# Round 2 fixes: theme_classic, mUM15 removal, ANOSIM subtitles,
#   samplename labels, filled shapes, summary heatmap data, none color
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(qs2)
  library(qs)
  library(openxlsx)
})

# Source color palettes (colors.R lives alongside this script in R/)
source("colors.R")

# Paths
base_dir <- ".."
gloscope_dir <- file.path(base_dir, "GloScope_Full_Outputs")
output_dir <- file.path(base_dir, "Revised_Plots_20260325", "GroupA_GloScope")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== p31662 GloScope Plot Revisions (Group A) - Round 2 ===\n")
cat("Output:", output_dir, "\n\n")

# Load sample metadata
sample_metadata <- read_csv(file.path(gloscope_dir, "sample_metadata.csv"), show_col_types = FALSE)

# Samplename mapping (from newest 51-sample object)
samplename_map <- c(
  "395751_25-hCh18_GEX_G6" = "hCh18", "395751_30-hCh20_GEX_G11" = "hCh20",
  "395751_32-hCh21_GEX_H1" = "hCh21", "395751_34-hCh22_GEX_H3" = "hCh22",
  "395751_36-hCh23_GEX_H5" = "hCh23", "395751_38-hCh24_GEX_H7" = "hCh24",
  "395751_40-hCh25_GEX_H9" = "hCh25", "395751_27-hIr19_GEX_G8" = "hIr19",
  "395751_08-mUM01_1_GEX_E9" = "mUM01_1", "395751_10-mUM02_GEX_E11" = "mUM02",
  "395751_12-mUM04_GEX_F5" = "mUM04", "395751_13-mUM05_GEX_F6" = "mUM05",
  "395751_14-mUM06_GEX_F7" = "mUM06", "395751_15-mUM07_GEX_F8" = "mUM07",
  "395751_16-mUM08_1_GEX_F9" = "mUM08_1", "395751_17-mUM08_2_GEX_F10" = "mUM08_2",
  "395751_18-mUM09_GEX_F11" = "mUM09", "395751_20-mUM11_GEX_G1" = "mUM11",
  "o28554_1_14-CNA_liv_bl_GEX_F10" = "mUM13_1", "o28554_1_13-FNA_liv_bl_GEX_F9" = "mUM13_2",
  "395751_48-mUM14_BL_GEX_E1" = "mUM14", "395751_51-mUM16_BL_GEX_E4" = "mUM16",
  "395751_53-mUM18_GEX_E6" = "mUM18", "395751_54-mUM19_GEX_E7" = "mUM19",
  "395751_55-mUM20_GEX_E8" = "mUM20",
  "341701_01-pUM1_637_GEX_H6" = "pUM01", "341701_02-pUM2_801_GEX_H7" = "pUM02",
  "341701_03-pUM3_669_GEX_H8" = "pUM03", "341701_04-pUM4_696_GEX_H12" = "pUM04",
  "341701_05-pUM5_713_GEX_A3" = "pUM05", "341701_06-pUM6_717_GEX_A4" = "pUM06",
  "373191_02-pUM9_GEX_F12" = "pUM09_1", "373191_03-pUM9_16S_GEX_G1" = "pUM09_2",
  "395751_22-pUM10_GEX_G3" = "pUM10",
  "373191_05-pUM11_GEX_G3" = "pUM11_1", "373191_06-pUM11_16S_GEX_G4" = "pUM11_2",
  "373191_07-pUM12_GEX_G5" = "pUM12_1", "373191_08-pUM12_16S_GEX_G6" = "pUM12_2",
  "373191_09-pUM14_GEX_G7" = "pUM14",
  "373191_10-pUM15_GEX_G8" = "pUM15_1", "373191_11-pUM15_16S_GEX_G9" = "pUM15_2",
  "395751_24-pUM16_GEX_G5" = "pUM16",
  "373191_12-pUM17_GEX_G10" = "pUM17_1", "373191_13-pUM17_16S_GEX_G11" = "pUM17_2",
  "395751_26-pUM19_GEX_G7" = "pUM19", "395751_29-pUM20_GEX_G10" = "pUM20",
  "395751_31-pUM21_GEX_G12" = "pUM21", "395751_33-pUM22_GEX_H2" = "pUM22",
  "395751_35-pUM23_GEX_H4" = "pUM23", "395751_37-pUM24_GEX_H6" = "pUM24",
  "395751_39-pUM25_GEX_H8" = "pUM25"
)

# Add samplename to metadata
sample_metadata$samplename <- samplename_map[sample_metadata$Sample]

# Helper: filter mUM15 from divergence matrix
filter_mum15 <- function(div_matrix) {
  mum15 <- grep("mUM15", rownames(div_matrix), value = TRUE)
  if (length(mum15) > 0) {
    cat("  Removing", length(mum15), "mUM15 sample(s) from divergence matrix\n")
    keep <- !rownames(div_matrix) %in% mum15
    div_matrix <- div_matrix[keep, keep]
  }
  div_matrix
}

# Helper: format ANOSIM p-value for display
format_anosim <- function(anosim_df) {
  p_display <- ifelse(anosim_df$p_value < 0.001, "p < 0.001",
                      paste0("p = ", format(anosim_df$p_value, digits = 3)))
  paste0("ANOSIM: R = ", round(anosim_df$statistic, 3), ", ", p_display)
}

# Helper: build MDS ggplot from divergence matrix using cmdscale
build_mds_plot <- function(div_matrix, metadata, group_col, color_vals,
                           title, subtitle = NULL, point_size = 6) {
  div_matrix <- filter_mum15(div_matrix)
  mds_coords <- cmdscale(as.dist(div_matrix), k = 2)
  df <- data.frame(
    Sample = rownames(mds_coords),
    dim1 = mds_coords[, 1],
    dim2 = mds_coords[, 2],
    stringsAsFactors = FALSE
  ) |>
    left_join(metadata, by = "Sample")

  p <- ggplot(df, aes(x = dim1, y = dim2, color = .data[[group_col]])) +
    geom_point(size = point_size) +
    scale_color_manual(values = color_vals, na.value = "#CCCCCC") +
    coord_fixed() +
    xlab("MDS1") + ylab("MDS2") +
    labs(title = title, subtitle = subtitle) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )

  list(plot = p, df = df)
}

# ============================================================================
# A1: Export ANOSIM Summary heatmap data (NOT Divergence heatmap)
# ============================================================================
cat("A1: Exporting ANOSIM summary heatmap data...\n")

anosim_summary <- read_csv(file.path(gloscope_dir, "anosim_summary_full_cohort.csv"), show_col_types = FALSE)
write.csv(anosim_summary, file.path(output_dir, "A1_anosim_summary_data.csv"), row.names = FALSE)
write.xlsx(anosim_summary, file.path(output_dir, "A1_anosim_summary_data.xlsx"))
cat("  Saved:", nrow(anosim_summary), "ANOSIM comparisons\n")

# ============================================================================
# A7: GloScope Divergence Heatmap Overview
# ============================================================================
cat("A7: Rebuilding divergence heatmap with corrected annotations...\n")

heatmap_objects <- tryCatch(
  qs2::qs_read(file.path(gloscope_dir, "gloscope_heatmap_all_cells_objects.qs2"), nthreads = 8),
  error = function(e) {
    cat("  qs2 failed, trying qs::qread...\n")
    qs::qread(file.path(gloscope_dir, "gloscope_heatmap_all_cells_objects.qs2"), nthreads = 8)
  }
)

div_matrix <- heatmap_objects$div_matrix
anno_data <- heatmap_objects$anno_data

# Remove mUM15 from heatmap
mum15_idx <- grep("mUM15", rownames(div_matrix))
if (length(mum15_idx) > 0) {
  cat("  Removing mUM15 from heatmap\n")
  div_matrix <- div_matrix[-mum15_idx, -mum15_idx]
  anno_data <- anno_data[-mum15_idx, ]
}

# Compute age ranges for relabeling
age_ranges <- sample_metadata |>
  filter(!is.na(Age_numeric), Age_group != "Unknown") |>
  group_by(Age_group) |>
  summarise(min_age = floor(min(Age_numeric)), max_age = ceiling(max(Age_numeric)), .groups = "drop")

age_relabel <- setNames(
  paste0(age_ranges$Age_group, " [", age_ranges$min_age, "-", age_ranges$max_age, "]"),
  age_ranges$Age_group
)
age_relabel <- c(age_relabel, "Unknown" = "Unknown")

# Corrected secondary mutation per samplename (from Seurat object — splits "none" vs "unknown")
smut_corrected_by_samplename <- c(
  # BAP1
  "hIr19"   = "BAP1", "mUM01_1" = "BAP1", "mUM04"   = "BAP1", "mUM05"   = "BAP1",
  "mUM06"   = "BAP1", "mUM07"   = "BAP1", "mUM08_1" = "BAP1", "mUM08_2" = "BAP1",
  "mUM09"   = "BAP1", "mUM11"   = "BAP1", "mUM13_1" = "BAP1", "mUM13_2" = "BAP1",
  "mUM14"   = "BAP1", "mUM18"   = "BAP1", "mUM19"   = "BAP1", "mUM20"   = "BAP1",
  "pUM09_1" = "BAP1", "pUM09_2" = "BAP1", "pUM11_1" = "BAP1", "pUM11_2" = "BAP1",
  "pUM19"   = "BAP1",
  # EIF1AX
  "pUM02"   = "EIF1AX", "pUM16"   = "EIF1AX", "pUM23"   = "EIF1AX", "pUM24"   = "EIF1AX",
  # SF3B1
  "pUM05"   = "SF3B1", "pUM06"   = "SF3B1", "pUM10"   = "SF3B1", "pUM20"   = "SF3B1",
  "pUM21"   = "SF3B1", "pUM22"   = "SF3B1", "pUM25"   = "SF3B1",
  # none (no secondary mutation)
  "hCh18"   = "none", "hCh20" = "none", "hCh21" = "none", "hCh22" = "none",
  "hCh23"   = "none", "hCh24" = "none", "hCh25" = "none",
  "mUM02"   = "none", "mUM16" = "none", "pUM04" = "none",
  # unknown (not assessed)
  "pUM01"   = "unknown", "pUM03"   = "unknown",
  "pUM12_1" = "unknown", "pUM12_2" = "unknown",
  "pUM14"   = "unknown", "pUM15_1" = "unknown", "pUM15_2" = "unknown",
  "pUM17_1" = "unknown", "pUM17_2" = "unknown"
)

# Build Sample → secondary_mutation_corrected via samplename_map
smut_corrected <- setNames(
  smut_corrected_by_samplename[samplename_map],
  names(samplename_map)
)

# Clean annotation data
anno_clean <- anno_data |>
  mutate(
    Gender = ifelse(is.na(Gender) | Gender %in% c("NA", ""), "Unknown", as.character(Gender)),
    Age_group = ifelse(is.na(Age_group) | Age_group %in% c("NA", ""), "Unknown", as.character(Age_group)),
    Age_group_label = age_relabel[Age_group],
    # Restore split none/unknown using per-sample corrected values (use Sample column, not rownames)
    secondary_mutation = coalesce(smut_corrected[Sample], as.character(secondary_mutation)),
    primary_mutation = ifelse(
      is.na(primary_mutation) | primary_mutation %in% c("unknown", "Unknown", "NA", ""),
      "Unknown", as.character(primary_mutation)
    )
  )

unknown_grey <- "#E1E1E1FF"
gender_colors_full <- c(gender_col, "Unknown" = unknown_grey)
age_label_colors <- setNames(
  unname(ageGroup_col[age_ranges$Age_group]),
  age_relabel[age_ranges$Age_group]
)
age_label_colors <- c(age_label_colors, "Unknown" = unknown_grey)
pmut_colors_clean <- c(primary_mutations_col, "Unknown" = unknown_grey)
# Use proper colors from second_mutations_col (none and unknown are separate)
smut_colors_clean <- c(second_mutations_col, "Unknown" = unknown_grey)

all_patients <- unique(anno_clean$Patient_nr)
patient_colors_full <- patientID_col[all_patients]
missing_pts <- all_patients[is.na(patient_colors_full)]
if (length(missing_pts) > 0) {
  patient_colors_full[missing_pts] <- scales::hue_pal()(length(missing_pts))
}
patient_colors_full <- patient_colors_full[!is.na(names(patient_colors_full))]

col_anno <- HeatmapAnnotation(
  Patient = anno_clean$Patient_nr,
  Origin = anno_clean$origin,
  Organ = anno_clean$organ,
  Gender = anno_clean$Gender,
  `Age Group` = anno_clean$Age_group_label,
  `Primary Mutation` = anno_clean$primary_mutation,
  `Secondary Mutation` = anno_clean$secondary_mutation,
  col = list(
    Patient = patient_colors_full,
    Origin = origin_colors,
    Organ = organ_colors,
    Gender = gender_colors_full,
    `Age Group` = age_label_colors,
    `Primary Mutation` = pmut_colors_clean,
    `Secondary Mutation` = smut_colors_clean
  ),
  annotation_name_gp = gpar(fontsize = 11),
  show_legend = c(Patient = FALSE, Origin = TRUE, Organ = TRUE,
                  Gender = TRUE, `Age Group` = TRUE,
                  `Primary Mutation` = TRUE, `Secondary Mutation` = TRUE)
)

ht <- Heatmap(
  div_matrix,
  name = "KL Divergence",
  col = colorRamp2(c(0, median(div_matrix), max(div_matrix)), c("white", "orange", "red")),
  top_annotation = col_anno,
  show_row_names = FALSE, show_column_names = FALSE,
  column_title = "GloScope Divergence: Disease Origin (Full Cohort - All Cells)",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_param = list(title = "KL Divergence", legend_height = unit(4, "cm"))
)

cairo_pdf(file.path(output_dir, "A7_heatmap_overview.pdf"), width = 14, height = 12)
draw(ht)
dev.off()
cat("  Saved: A7_heatmap_overview.pdf\n")

# Export heatmap data for collaborators (updated: mUM15 removed, none/unknown split)
write.csv(div_matrix, file.path(output_dir, "A1_heatmap_divergence_matrix.csv"))
write.csv(anno_clean, file.path(output_dir, "A1_heatmap_annotation_data.csv"), row.names = FALSE)
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Divergence Matrix")
openxlsx::writeData(wb, "Divergence Matrix", div_matrix, rowNames = TRUE)
openxlsx::addWorksheet(wb, "Annotation")
openxlsx::writeData(wb, "Annotation", anno_clean)
openxlsx::saveWorkbook(wb, file.path(output_dir, "A1_heatmap_data.xlsx"), overwrite = TRUE)
cat("  Exported heatmap data: divergence matrix + annotation (51 samples, corrected)\n")

rm(heatmap_objects, div_matrix, anno_data, anno_clean, ht, col_anno)
gc()

# ============================================================================
# A2: MDS Anatomical Location in Primary UM
# ============================================================================
cat("A2: MDS Anatomical Location in Primary UM...\n")

# Use full cohort matrix and subset to primary samples to include pUM06 (ciliochoroidal)
# pUM06 was excluded from divergence_location_pUM_full.rds because its metadata
# had LocationPrimaryTumor = "other/NA" (incorrect). We fix this here.
div_full <- readRDS(file.path(gloscope_dir, "divergence_origin_full_cohort.rds"))
div_full <- filter_mum15(div_full)

# Get primary UM samples
pum_samples <- sample_metadata$Sample[sample_metadata$origin == "primary"]
pum_in_mat <- intersect(pum_samples, rownames(div_full))
div_loc <- div_full[pum_in_mat, pum_in_mat]

meta_loc <- data.frame(Sample = rownames(div_loc)) |>
  left_join(sample_metadata, by = "Sample") |>
  # Fix pUM06: its primary_location is "ciliochoroidal" but metadata has "other/NA"
  mutate(LocationPrimaryTumor = case_when(
    Sample == "341701_06-pUM6_717_GEX_A4" ~ "ciliochoroidal",
    # Fix typo "chiliochoroidal" → "ciliochoroidal" if present
    LocationPrimaryTumor == "chiliochoroidal" ~ "ciliochoroidal",
    TRUE ~ LocationPrimaryTumor
  ))

location_colors_all <- c(primary_location_colors, "other/NA" = "#CCCCCC")

anosim_loc <- read_csv(file.path(gloscope_dir, "anosim_location_pUM_full.csv"), show_col_types = FALSE)

result <- build_mds_plot(div_loc, meta_loc, "LocationPrimaryTumor", location_colors_all,
                         title = "GloScope MDS: Anatomical Location in Primary UM",
                         subtitle = format_anosim(anosim_loc))
p <- result$plot + labs(color = "Location")

ggsave(file.path(output_dir, "A2_MDS_location_pUM.pdf"), p,
       device = cairo_pdf, width = 10, height = 8)
cat("  Saved: A2_MDS_location_pUM.pdf\n")

# ============================================================================
# A3: MDS Secondary Mutations in Primary UM (Full cohort)
# ============================================================================
cat("A3: MDS Secondary Mutations in Primary UM...\n")

div_smut <- readRDS(file.path(gloscope_dir, "divergence_secondary_mutation_pUM_full.rds"))
meta_smut <- data.frame(Sample = rownames(div_smut)) |>
  left_join(sample_metadata, by = "Sample")

smut_colors <- c(second_mutations_col, "none/unknown" = "#CCCCCC", "NA" = "#CCCCCC")
anosim_smut <- read_csv(file.path(gloscope_dir, "anosim_secondary_mutation_pUM_full.csv"), show_col_types = FALSE)

result <- build_mds_plot(div_smut, meta_smut, "secondary_mutation", smut_colors,
                         title = "GloScope MDS: Secondary Mutations in Primary UM (Full Cohort)",
                         subtitle = format_anosim(anosim_smut))
p <- result$plot + labs(color = "Secondary Mutation")

ggsave(file.path(output_dir, "A3_MDS_secondary_mutation_pUM.pdf"), p,
       device = cairo_pdf, width = 10, height = 8)
cat("  Saved: A3_MDS_secondary_mutation_pUM.pdf\n")

# ============================================================================
# A4: MDS Disease Progression (Full cohort)
# ============================================================================
cat("A4: MDS Disease Progression (Full cohort)...\n")

div_origin <- readRDS(file.path(gloscope_dir, "divergence_origin_full_cohort.rds"))
anosim_origin <- read_csv(file.path(gloscope_dir, "anosim_origin_full.csv"), show_col_types = FALSE)

meta_origin <- data.frame(Sample = rownames(div_origin)) |>
  left_join(sample_metadata, by = "Sample")

result <- build_mds_plot(div_origin, meta_origin, "origin", origin_colors,
                         title = "GloScope MDS: Disease Progression (Full Cohort)",
                         subtitle = format_anosim(anosim_origin))
p <- result$plot + labs(color = "Origin")

ggsave(file.path(output_dir, "A4_MDS_disease_progression.pdf"), p,
       device = cairo_pdf, width = 10, height = 8)
cat("  Saved: A4_MDS_disease_progression.pdf\n")

# ============================================================================
# A5: MDS Primary Mutations (Full cohort)
# ============================================================================
cat("A5: MDS Primary Mutations (Full cohort)...\n")

div_pmut <- readRDS(file.path(gloscope_dir, "divergence_primary_mutation_full.rds"))
meta_pmut <- data.frame(Sample = rownames(div_pmut)) |>
  left_join(sample_metadata, by = "Sample")

pmut_colors <- c(primary_mutations_col, "NA" = "#CCCCCC")
anosim_pmut <- read_csv(file.path(gloscope_dir, "anosim_primary_mutation_full.csv"), show_col_types = FALSE)

result <- build_mds_plot(div_pmut, meta_pmut, "primary_mutation", pmut_colors,
                         title = "GloScope MDS: Primary Mutations (Full Cohort)",
                         subtitle = format_anosim(anosim_pmut))
p <- result$plot + labs(color = "Primary Mutation")

ggsave(file.path(output_dir, "A5_MDS_primary_mutations.pdf"), p,
       device = cairo_pdf, width = 10, height = 8)
cat("  Saved: A5_MDS_primary_mutations.pdf\n")

# ============================================================================
# A6: MDS Organ/Metastasis Site (All cells)
# ============================================================================
cat("A6: MDS Organ/Metastasis Site (All cells)...\n")

div_organ <- readRDS(file.path(gloscope_dir, "divergence_organ_full_cohort.rds"))
div_organ <- filter_mum15(div_organ)

meta_organ <- data.frame(Sample = rownames(div_organ)) |>
  left_join(sample_metadata, by = "Sample")

anosim_organ <- read_csv(file.path(gloscope_dir, "anosim_organ_full.csv"), show_col_types = FALSE)

mds_coords_organ <- cmdscale(as.dist(div_organ), k = 2)
df_organ <- data.frame(
  Sample = rownames(mds_coords_organ),
  dim1 = mds_coords_organ[, 1],
  dim2 = mds_coords_organ[, 2],
  stringsAsFactors = FALSE
) |>
  left_join(meta_organ, by = "Sample") |>
  mutate(
    samplename = samplename_map[Sample],
    origin = factor(origin, levels = c("healthy", "primary", "metastasis"))
  )

# Filled shapes per UMscAtlas README
progression_shapes <- c("healthy" = 16, "primary" = 17, "metastasis" = 18)

p <- ggplot(df_organ, aes(x = dim1, y = dim2, color = organ, shape = origin)) +
  geom_point(size = 6) +
  geom_text_repel(aes(label = samplename), size = 3, max.overlaps = 25,
                  show.legend = FALSE) +
  scale_color_manual(values = organ_colors, name = "Organ") +
  scale_shape_manual(values = progression_shapes, name = "Disease Progression") +
  coord_fixed() +
  xlab("MDS1") + ylab("MDS2") +
  labs(title = "GloScope MDS: Organ/Metastasis Site (All Cells)",
       subtitle = format_anosim(anosim_organ)) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 5))) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

ggsave(file.path(output_dir, "A6_MDS_organ_all_cells.pdf"), p,
       device = cairo_pdf, width = 12, height = 9)
cat("  Saved: A6_MDS_organ_all_cells.pdf\n")

# ============================================================================
cat("\n=== Group A complete! ===\n")
cat("Output directory:", output_dir, "\n")
cat("Files:\n")
cat(paste(" ", list.files(output_dir), collapse = "\n"), "\n")

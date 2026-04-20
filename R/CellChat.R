# This script was used to create the following figures:
# Figures 4h-i

library(tidyverse)
library(Seurat)
library(qs)
library(pals)
library(CellChat)
library(patchwork)
set.seed(1234)

# define colors (see "color_palettes.Rmd")
majority_celltype_col <- c(
  "Melanoma cells" = "#B00068",
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

# select only primary uveal melanoma samples
Idents(seu) <- "origin"
seu <- subset(seu, idents = "primary")

# select primary tumors with known secondary mutations
Idents(seu) <- "secondary_mutation"
seu <- subset(seu, idents = c("none", "unknown"), invert = TRUE)

# select only celltypes that are present in all conditions
types_use <- c("Melanoma cells", "Endothelial cells", "Pericytes", "Myeloid cells", "T/NK cells", "Glial cells", "B/Plasma cells")
Idents(seu) <- "majority_celltype"
seu <- subset(seu, idents = types_use)

# make function to run CellChat preprocessing
prep_cellchat_obj <- function(seurat_object, celltype_column){
  
  # prepare input
  data.input <- seurat_object[["RNA"]]$data
  Idents(seurat_object) <- celltype_column
  labels <- Idents(seurat_object)
  meta <- data.frame(labels = labels, row.names = names(labels))
  
  # create CellChat object
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
  
  # define data base
  CellChatDB <- CellChatDB.human
  cellchat@DB <- subsetDB(CellChatDB)
  
  # preprocessing
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # compute communication probability 
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  # compute network centrality scores
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "net")
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  return(cellchat)
  
}

# split Seurat object by secondary mutation
seu_list <- SplitObject(seu, split.by = "secondary_mutation")

# run CellChat function for each element in list
chat_list <- c()

for (obj in seu_list){
  
  chat_obj <- prep_cellchat_obj(seurat_object = obj, celltype_column = "majority_celltype")
  
  chat_list[[unique(obj$secondary_mutation)]] <- chat_obj
  
}

# merge into one CellChat object
# this is necessary because CellChat uses a list of objects as input for some functions and the merged object as input for some others
cellchat <- mergeCellChat(chat_list, add.names = names(chat_list))

# order color vector to make sure it's consistent between plots
celltype_col_ordered <- majority_celltype_col[order(factor(names(majority_celltype_col), levels = levels(chat_list[[1]]@idents)))]

# Fig. 4h
netVisual_diffInteraction(
  cellchat, 
  weight.scale = TRUE, 
  measure = "weight", 
  comparison = c("EIF1AX", "BAP1"), 
  color.use = celltype_col_ordered, 
  alpha.edge = 0.8, 
  arrow.size = 0.5, 
  arrow.width = 1.5, 
  title.name = "Differential Interaction Strength - BAP1 vs. EIF1AX", 
  margin = 0.5
  )

netVisual_diffInteraction(
  cellchat, 
  weight.scale = TRUE, 
  measure = "weight", 
  comparison = c("SF3B1", "BAP1"), 
  color.use = celltype_col_ordered, 
  alpha.edge = 0.8, 
  arrow.size = 0.5, 
  arrow.width = 1.5, 
  title.name = "Differential Interaction Strength - BAP1 vs. SF3B1", 
  margin = 0.5
  )

# Fig. 4i
sources <- c("Melanoma cells", "T/NK cells", "Myeloid cells")

netVisual_bubble(
  cellchat, 
  sources.use = sources, 
  targets.use = c("T/NK cells"), 
  comparison = c(1, 2, 3), 
  angle.x = 45, thresh = 0.01, 
  signaling = "MHC-I", 
  remove.isolate = FALSE, 
  line.on = TRUE, 
  color.text = c("#7B1C5DFF", "#EAC541FF", "#1F867BFF")
  ) +
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 30))

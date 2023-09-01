#==============================================================================#
# Author : Angela M. Oill, aoill@tgen.org
# Date: 2023/06/14
# Project: Pediatric CAR-T 
# Description: Figure S06
#==============================================================================#

#==============================================================================#
# Load libraries ----
#==============================================================================#
library(ggplot2)
library(Seurat)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(dittoSeq)
library(scProportionTest)
library(viridis)
library(ComplexHeatmap)


#==============================================================================#
# Read in integrated object ----
#==============================================================================#
product_obj_T <- readRDS("/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/product_T_cell_obj_2023.rds")


#==============================================================================#
# Set up cluster colors ----
#==============================================================================#
t_cell_colors <- c("CD4 Naive" = "#FEEA9A", "CD4 Effector" = "#FECD6A",
                   "CD4 Memory" = "#FDA245", "CD4 Memory (LTB+)" = "#FC6931", 
                   "CD4 Activated" = "#EA2820",
                   "CD4 Proliferating" = "#C20324", "Treg" = "#800026", 
                   
                   "CD8 Naive" = "#D5E3EF", "CD8 Effector" = "#A9C3DE",
                   "CD8 Effector Memory" = "#8C96C6", "CD8 Memory" = "#8A5DAA", 
                   "CD8 Resident Memory-like" = "#831F87",  "CD8 Central Memory (Activated)" = "#4D004B",
                   "CD8 Activated" = "#D3EDCC", "CD8 Proliferating" = "#98D493", 
                   "CD8 Exhausted" = "#4BAF61", 
                   "CD8 T (OX40+)" = "#147E3A", "NKT" = "#00441B"
                   
                   
)


product_tmp_CTs <- as.character(levels(as.factor(product_obj_T@meta.data$ct_T_states)))
product_tmp_color <- setNames(names(t_cell_colors), t_cell_colors) # swap names and values so I can subset
product_tmp_color_sub <- product_tmp_color[product_tmp_color %in% c(product_tmp_CTs)] # subset ct_color to only keep colors in CSF UMAP
product_t_cell_colors_sub <- setNames(names(product_tmp_color_sub), product_tmp_color_sub) # re-swap names and values so that values are the color names 


#==============================================================================#
# Plot ----
#==============================================================================#

#-----------#
## UMAP ----
#-----------#
a <- DimPlot(product_obj_T, group.by = "ct_T_states", 
             reduction = "integrated_sct_umap") +
  scale_color_manual(values = product_t_cell_colors_sub) +
  #ggtitle("") + 
  ggtitle("Product") + 
  NoLegend() #+
#theme(plot.title = element_text(hjust = 0.5, face = "bold",
#                                size = 10))
pa <- LabelClusters(a, id = "ct_T_states", box = TRUE, size = 3.5, label.size = 0.1, 
                    box.padding = 0.5, seed = 309,
                    color = c(rep("black", 5), rep("white", 4)),
                    alpha = 0.8
)

#pa

#----------------------#
## Stacked bar plot ----
#----------------------#
# Add a col to metadata for responding patients and non responding patients
# responding patients (514, 515, 626) and non-responding patients (574, 625)
product_obj_T@meta.data$patient_category <- NA # or any other initialization value
product_obj_T@meta.data$patient_category[product_obj_T@meta.data$UPN == "514"] <- "Responding" 
product_obj_T@meta.data$patient_category[product_obj_T@meta.data$UPN == "515"] <- "Responding" 
product_obj_T@meta.data$patient_category[product_obj_T@meta.data$UPN == "626"] <- "Responding" 
product_obj_T@meta.data$patient_category[product_obj_T@meta.data$UPN == "574"] <- "Non-Responding" 
product_obj_T@meta.data$patient_category[product_obj_T@meta.data$UPN == "625"] <- "Non-Responding" 


pb <- dittoBarPlot(
  object = product_obj_T,
  var = "ct_T_states",
  group.by = "UPN", 
  color.panel = product_t_cell_colors_sub,
  #main = "Product",
  main = "",
  x.reorder = c(3,4, 1,2,5)
) +
  xlab("UPN") +
  NoLegend() +
  theme(#legend.position="bottom",
    plot.title = element_text(hjust = 0.5, face = "bold",
                              size = 10)) +
  geom_vline(xintercept = 2.5, linetype="dashed", size=1)


#--------------------#
## Arrange plots ----
#--------------------#
pdf("/home/aoill/plots/Figure_S07.pdf",
    width = 7.5,
    height = 4)
ggarrange(pa, pb, 
          ncol = 2,
          nrow = 1,
          labels = c("A", "B"),
          widths = c(1, .7))
dev.off()

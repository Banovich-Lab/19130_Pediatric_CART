#==============================================================================#
# Author : Angela M. Oill, aoill@tgen.org
# Date: 2023/06/14
# Project: Pediatric CAR-T 
# Description: Figure S08
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
csf_obj_T <- readRDS("/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/CSF_T_cell_obj_2023.rds")


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


csf_tmp_CTs <- as.character(levels(as.factor(csf_obj_T@meta.data$ct_T_states)))
csf_tmp_color <- setNames(names(t_cell_colors), t_cell_colors) # swap names and values so I can subset
csf_tmp_color_sub <- csf_tmp_color[csf_tmp_color %in% c(csf_tmp_CTs)] # subset ct_color to only keep colors in CSF UMAP
csf_t_cell_colors_sub <- setNames(names(csf_tmp_color_sub), csf_tmp_color_sub) # re-swap names and values so that values are the color names 


#==============================================================================#
# Plot ----
#==============================================================================#
pa <- dittoBarPlot(
  object = csf_obj_T,
  var = "ct_T_states",
  group.by = "CART", 
  color.panel = csf_t_cell_colors_sub,
  main = "CSF"#,
  #x.reorder = c(3,4, 1,2,5)
) +
  xlab("CAR Positivity") +
  NoLegend() +
  theme(legend.position="right",
        plot.title = element_text(hjust = 0.5, face = "bold",
                                  size = 12)) 


prop_test <- sc_utils(csf_obj_T)
prop_test <- permutation_test(
  prop_test, 
  cluster_identity = "ct_T_states",
  sample_1 = "Negative", sample_2 = "Positive",
  sample_identity = "CART"
)

prop_test@results$permutation

pb <- permutation_plot(prop_test) + 
  ggtitle("CSF\nCAR-T Negative vs CAR-T Positive") +
  ylab("") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5, face = "bold",
                                  size = 12)) 

ggarrange(pa, pb)


## Arrange plots ----
pdf("/home/aoill/plots/Figure_S08.pdf",
    width = 9.5,
    height = 4)
ggarrange(pa, pb,
          ncol = 2,
          nrow = 1,
          labels = c("A", "B"))
dev.off()
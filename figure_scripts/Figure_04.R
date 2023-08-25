#==============================================================================#
# Author : Angela M. Oill, aoill@tgen.org
# Date: 2023/06/27
# Project: Pediatric CAR-T 
# Description: Figure 4
#==============================================================================#


#==============================================================================#
# Load libraries ----
#==============================================================================#
library(Seurat)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(ggpubr)
library(dittoSeq)
library(scProportionTest)
library(viridis)
library(scRepertoire)
library(VennDiagram)


#==============================================================================#
# Read in integrated objects ----
#==============================================================================#
csf_obj_T <- readRDS("/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/CSF_T_cell_obj_2023.rds")
pbmc_obj_T <- readRDS("/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/PBMC_T_cell_obj_2023.rds")


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

# I don't want all of the labels to come up on the legend so use the code
# below to subset the ct_color character vector
csf_tmp_CTs <- as.character(levels(as.factor(csf_obj_T@meta.data$ct_T_states)))
csf_tmp_color <- setNames(names(t_cell_colors), t_cell_colors) # swap names and values so I can subset
csf_tmp_color_sub <- csf_tmp_color[csf_tmp_color %in% c(csf_tmp_CTs)] # subset ct_color to only keep colors in CSF UMAP
csf_t_cell_colors_sub <- setNames(names(csf_tmp_color_sub), csf_tmp_color_sub) # re-swap names and values so that values are the color names 

pbmc_tmp_CTs <- as.character(levels(as.factor(pbmc_obj_T@meta.data$ct_T_states)))
pbmc_tmp_color <- setNames(names(t_cell_colors), t_cell_colors) # swap names and values so I can subset
pbmc_tmp_color_sub <- pbmc_tmp_color[pbmc_tmp_color %in% c(pbmc_tmp_CTs)] # subset ct_color to only keep colors in CSF UMAP
pbmc_t_cell_colors_sub <- setNames(names(pbmc_tmp_color_sub), pbmc_tmp_color_sub) # re-swap names and values so that values are the color names 


#==============================================================================#
# Figure 4 ----
#==============================================================================#
# CSF and PBMC UMAPs, point range plots

#------------#
## UMAPS ----
#------------#

### CSF ----
a <- DimPlot(csf_obj_T, group.by = "ct_T_states", 
             reduction = "integrated_sct_umap") +
  scale_color_manual(values = csf_t_cell_colors_sub) +
  ggtitle("CSF - T cells") + NoLegend() #+
#theme(plot.title = element_text(hjust = 0.5, face = "bold",
#                                size = 10))
pa <- LabelClusters(a, id = "ct_T_states", box = TRUE, size = 3.5, label.size = 0.1, 
                    box.padding = 0.5, seed = 309,
                    color = c(rep("black", 3), rep("white", 2), rep("black", 5)),
                    alpha = 0.8
)


### PBMC ----
b <- DimPlot(pbmc_obj_T, group.by = "ct_T_states", 
             reduction = "integrated_sct_umap") +
  scale_color_manual(values = pbmc_t_cell_colors_sub) +
  ggtitle("PBMC - T cells") + NoLegend() #+
#theme(plot.title = element_text(hjust = 0.5, face = "bold",
#                                size = 10))
pb <- LabelClusters(b, id = "ct_T_states", box = TRUE, size = 3.5, label.size = 0.1, 
                    box.padding = 0.5, seed = 309,
                    alpha = 0.8
)


#------------------------#
## Point range plots ----
#------------------------#
### CSF non-res vs res ----
# Proportion test
prop_test <- sc_utils(csf_obj_T)
prop_test <- permutation_test(
  prop_test, 
  cluster_identity = "ct_T_states",
  sample_1 = "Non-Response", sample_2 = "Response",
  sample_identity = "ResponseNon"
)

prop_test@results$permutation # SAVE THIS FOR MANUSCRIPT (PVALUES)
results <- prop_test@results$permutation
#write.csv(results, "/home/aoill/tables/sc_proportion_test_CSF_T_cells_res_non.csv",
#          row.names = F, quote = F)
pc <- permutation_plot(prop_test) + 
  ggtitle("CSF\nNon-Response vs Response") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5, vjust = 2, face = "bold",
                                  size = 10)) 

### PBMCs non-lymph vs lymph ----
prop_test <- sc_utils(pbmc_obj_T)
prop_test <- permutation_test(
  prop_test, 
  cluster_identity = "ct_T_states",
  sample_1 = "Non-Lymphodepleted", sample_2 = "Lymphodepleted",
  sample_identity = "Lymphodepletion"
)

prop_test@results$permutation # SAVE THIS FOR MANUSCRIPT (PVALUES)
results <- prop_test@results$permutation
#write.csv(results, "/home/aoill/tables/sc_proportion_test_PBMC_T_cells_lymph_non.csv",
#          row.names = F, quote = F)
#pe <- permutation_plot(prop_test) + 
pd <- permutation_plot(prop_test) + 
  ggtitle("PBMC\nNon-Lymphodepleted vs Lymphodepleted") +
  ylab("") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5, face = "bold",
                                  size = 10)) 

#---------------------#
## Arrange figure ----
#---------------------#

pdf("/home/aoill/plots/Figure_04.pdf", width = 9.5, height = 9
)
ggarrange(ggarrange(pa, pb, ncol = 2, nrow = 1, align = "hv",
                    labels = c("A", "B")),
          ggarrange(pc,
                    pd,
                    ncol = 2, 
                    common.legend = T,
                    legend = "bottom",
                    align = "hv",
                    labels = c("C", "D")),
          ncol = 1, 
          nrow = 2,
          heights = c(1,.8))
dev.off()


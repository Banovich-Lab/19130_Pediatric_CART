#==============================================================================#
# Author : Angela M. Oill, aoill@tgen.org
# Date: 2023/06/09
# Project: Pediatric CAR-T 
# Description: Figure 3
#==============================================================================#

#==============================================================================#
# Load libraries ----
#==============================================================================#
library(Seurat)
library(ggplot2)
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

## CSF ----
# Read in integrated object
csf_obj <- readRDS("/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/CSF_all_cells_obj_2023.rds")

## PBMCs ----
# Read in integrated object
pbmc_obj <- readRDS("/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/PBMC_all_cells_obj_2023.rds")


#==============================================================================#
# Set up cluster colors ----
#==============================================================================#
ct_colors <- c("CD4+ T" = "#e09cb7", "CD8+ T" = "#d747ad", "NK" = "#a35776", 
               "Treg" = "#d54170", 
               "Monocytes" = "#bcd095", "CD14+ Monocytes" = "#6dce52", "FCGR3A+ Monocytes" = "#697f3f",
               "Macrophages" = "#c8d048",
               "DC" = "#de9e35","cDC" = "#de9e35", "cDC1" = "#927754", "cDC2" = "#e3bb87", 
               "mDC" = "#a0712b",
               "pDC" = "#55daa2",  
               "B cells" = "#4a8773", "B" = "#4a8773", 
               "Plasma" = "#89d7c8", 
               "Proliferating" = "#c7624f"
               
)


# I don't want all of the labels to come up on the legend so use the code
# below to subset the ct_color character vector
csf_tmp_CTs <- as.character(levels(as.factor(csf_obj@meta.data$CT)))
csf_tmp_color <- setNames(names(ct_colors), ct_colors) # swap names and values so I can subset
csf_tmp_color_sub <- csf_tmp_color[csf_tmp_color %in% c(csf_tmp_CTs)] # subset ct_color to only keep colors in CSF UMAP
csf_ct_colors_sub <- setNames(names(csf_tmp_color_sub), csf_tmp_color_sub) # re-swap names and values so that values are the color names 

pbmc_tmp_CTs <- as.character(levels(as.factor(pbmc_obj@meta.data$CT)))
pbmc_tmp_color <- setNames(names(ct_colors), ct_colors) # swap names and values so I can subset
pbmc_tmp_color_sub <- pbmc_tmp_color[pbmc_tmp_color %in% c(pbmc_tmp_CTs)] # subset ct_color to only keep colors in CSF UMAP
pbmc_ct_colors_sub <- setNames(names(pbmc_tmp_color_sub), pbmc_tmp_color_sub) # re-swap names and values so that values are the color names 


#==============================================================================#
# UMAPs ----
#==============================================================================#
#------------#
## CSF ----
#------------#
a <- DimPlot(csf_obj, group.by = "CT", 
             reduction = "integrated_sct_umap") +
  scale_color_manual(values = csf_ct_colors_sub) +
  ggtitle("CSF") + NoLegend()
pa <- LabelClusters(a, id = "CT", box = TRUE, size = 3.5, label.size = 0.1, 
                    box.padding = 0.5, seed = 309,
                    color = c(rep("black", 1), rep("white", 1), rep("white", 2), rep("black", 2), 
                              rep("white", 2), rep("black", 2), rep("white", 1), 
                              rep("black", 1)),
                    alpha = 0.8
)


#------------#
## PBMCs ----
#------------#
b <- DimPlot(pbmc_obj, group.by = "CT", 
             reduction = "integrated_sct_umap") +
  scale_color_manual(values = pbmc_ct_colors_sub) +
  ggtitle("PBMC") + NoLegend()
pb <- LabelClusters(b, id = "CT", box = TRUE, size = 3.5, label.size = 0.1,
                    box.padding = 0.5, seed = 309,
                    color = c(rep("white", 4), rep("black", 3),
                              rep("white", 1), rep("black", 2)),
                    alpha = 0.8
)


#==============================================================================#
# Stacked Bar Plots ----
#==============================================================================#
#------------#
## CSF ----
#------------#
# Re-order split umap by cycle number 
csf_obj@meta.data$UPN_Cycle <- factor(x = csf_obj@meta.data$UPN_Cycle, 
                                      levels = c("514_5", "514_8", "514_11",
                                                 "515_4", "515_8", 
                                                 "574_3", "574_8"))
pc <- dittoBarPlot(
  object = csf_obj,
  var = "CT",
  group.by = "UPN_Cycle", 
  color.panel = csf_ct_colors_sub, 
  main = "CSF",
  x.reorder = c(2,3,1,4,5,6,7)
) + 
  geom_vline(xintercept = 3.5, linetype="dotted", size= 1) + #0.75
  geom_vline(xintercept = 5.5, linetype="dotted", size= 1) +
  NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


#------------#
## PBMCs ----
#------------#
# Re-order split umap by cycle number 
pbmc_obj@meta.data$UPN_Cycle <- factor(x = pbmc_obj@meta.data$UPN_Cycle, 
                                       levels = c("514_5", "514_8", "514_11",
                                                  "515_4", "515_8", 
                                                  "574_3", "574_8",
                                                  "625_1", "625_3", "625_4",
                                                  "626_1", "626_4", "626_7" ))

pd <- dittoBarPlot(
  object = pbmc_obj,
  var = "CT",
  group.by = "UPN_Cycle", 
  color.panel = pbmc_ct_colors_sub,
  main = "PBMC",
  x.reorder = c(2,3,1,4,5,6,7,8,9,10,11,12,13)
) + 
  geom_vline(xintercept = 3.5, linetype="dotted", size=1) + 
  geom_vline(xintercept = 5.5, linetype="dotted", size=1) + 
  geom_vline(xintercept = 7.5, linetype="dotted", size=1) +
  geom_vline(xintercept = 10.5, linetype="dotted", size=1) +
  NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


#==============================================================================#
# Proportion tests/point range plots ----
#==============================================================================#
#------------#
## CSF ----
#------------#
# Proportion test
prop_test_CSF <- sc_utils(csf_obj)
prop_test <- sc_utils(csf_obj)
prop_test <- permutation_test(
  prop_test, 
  cluster_identity = "CT",
  sample_1 = "Non-Response", sample_2 = "Response",
  sample_identity = "ResponseNon"
)

prop_test@results$permutation # SAVE THIS FOR MANUSCRIPT (PVALUES)
results <- prop_test@results$permutation
#write.csv(results, "/home/aoill/tables/sc_proportion_test_CSF_all_res_non.csv",
#          row.names = F, quote = F)
pe <- permutation_plot(prop_test) + 
  ggtitle("CSF\nNon-Response vs Response") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5, vjust = 2, face = "bold",
                                  size = 10)) 

#------------#
## PBMCs ----
#------------#
# Proportion test - Lymphodepletion
prop_test <- sc_utils(pbmc_obj)
prop_test <- permutation_test(
  prop_test, 
  cluster_identity = "CT",
  sample_1 = "Non-Lymphodepleted", sample_2 = "Lymphodepleted",
  sample_identity = "Lymphodepletion"
)

prop_test@results$permutation # SAVE THIS FOR MANUSCRIPT (PVALUES)
results <- prop_test@results$permutation
#write.csv(results, "/home/aoill/tables/sc_proportion_test_PBMC_all_lymph_non.csv",
#          row.names = F, quote = F)
pf <- permutation_plot(prop_test) + 
  ggtitle("PBMC\nNon-Lymphodepleted vs Lymphodepleted") +
  ylab("") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5, face = "bold",
                                  size = 10)) 


#==============================================================================#
# Output figure ----
#==============================================================================#
pdf("/home/aoill/plots/Figure_03.pdf", width = 8.5, height = 11)
ggarrange(
  ggarrange(pa, pb, 
            pc, pd,
            ncol = 2,
            nrow = 2,
            heights = c(1, .7),
            labels = c("A", "B", "C", "D")),
  ggarrange(pe, pf, #pg,
            ncol = 2, nrow = 1, align = "hv", 
            common.legend = T, legend = "bottom", 
            labels = c("E", "F")),
  nrow = 2,
  ncol = 1,
  heights = c(1, .6)
)
dev.off()


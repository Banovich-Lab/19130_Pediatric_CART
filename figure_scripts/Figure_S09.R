#==============================================================================#
# Author : Angela M. Oill, aoill@tgen.org
# Date: 2023/06/14
# Project: Pediatric CAR-T 
# Description: Figure S09
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
# Plot ----
#==============================================================================#
pdf("/home/aoill/plots/Figure_S09.pdf",
    width = 5,
    height = 5)
dittoBarPlot(
  object = csf_obj_T,
  var = "CART",
  group.by = "Frequency_norm_cat_2", 
  color.panel = c("grey", "red"),
  main = "CSF",
  x.reorder = c(2,1)
) +
  xlab("CAR Positivity") +
  NoLegend() +
  theme(#legend.position="bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")) 

dev.off()

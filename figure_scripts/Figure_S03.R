#==============================================================================#
# Author : Angela M. Oill, aoill@tgen.org
# Date: 2023/06/29
# Project: Pediatric CAR-T 
# Description: Figure S03
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
# Read in object ----
#==============================================================================#
obj <- readRDS("/scratch/aoill/projects/CAR-T/00_final/rds_files/All_CSF_PBMC_integrated.rds")
obj <- readRDS("/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/CSF_PBMC_all_cells_obj_2023.rds")


#==============================================================================#
# Unify cell type annotations ----
#==============================================================================#
unique(obj@meta.data$CT)

obj@meta.data$CT_new <- "NA" # or any other initialization value
obj@meta.data$CT_new[obj@meta.data$CT == "T cells"] <- "T cells" 
obj@meta.data$CT_new[obj@meta.data$CT == "CD8+ T"] <- "T cells" 
obj@meta.data$CT_new[obj@meta.data$CT == "CD4+ T"] <- "T cells" 
obj@meta.data$CT_new[obj@meta.data$CT == "Treg"] <- "T cells" 

obj@meta.data$CT_new[obj@meta.data$CT == "Immune"] <- "Meyloid" 
obj@meta.data$CT_new[obj@meta.data$CT == "Monocytes"] <- "Monocytes" 
obj@meta.data$CT_new[obj@meta.data$CT == "CD14+ Monocytes"] <- "Monocytes" 
obj@meta.data$CT_new[obj@meta.data$CT == "FCGR3A+ Monocytes"] <- "Monocytes" 

obj@meta.data$CT_new[obj@meta.data$CT == "Macrophages"] <- "Macrophages" 

obj@meta.data$CT_new[obj@meta.data$CT == "DC"] <- "cDCs" 
obj@meta.data$CT_new[obj@meta.data$CT == "cDC1"] <- "cDCs" 
obj@meta.data$CT_new[obj@meta.data$CT == "cDC2"] <- "cDCs" 
obj@meta.data$CT_new[obj@meta.data$CT == "cDCs"] <- "cDCs" 


obj@meta.data$CT_new[obj@meta.data$CT == "B cells"] <- "B cells" 
obj@meta.data$CT_new[obj@meta.data$CT == "B"] <- "B cells" 

obj@meta.data$CT_new[obj@meta.data$CT == "mDC"] <- "mDCs" 
obj@meta.data$CT_new[obj@meta.data$CT == "pDC"] <- "pDCs" 

obj@meta.data$CT_new[obj@meta.data$CT == "Plasma"] <- "Plasma" 

obj@meta.data$CT_new[obj@meta.data$CT == "Proliferating"] <- "Proliferating" 

obj@meta.data$CT_new[obj@meta.data$CT == "NK"] <- "NK" 


#==============================================================================#
# Set colors for plotting ----
#==============================================================================#
nclust <- length(levels(as.factor(obj@meta.data$CT_new)))
cluster_col <- colorRampPalette(brewer.pal(10, "Paired"))(nclust)
names(cluster_col) <- as.character(levels(as.factor(obj@meta.data$CT_new)))


#==============================================================================#
# Plot and output figure ----
#==============================================================================#
obj_csf <- subset(obj, subset = Sample_Type == "CSF")
a <- DimPlot(obj_csf, group.by = "CT_new", 
             reduction = "integrated_sct_umap",
             label = F) +
  scale_color_manual(values = cluster_col)  +
  ggtitle("CSF") + 
  NoLegend() 
pa <- LabelClusters(a, id = "CT_new", box = TRUE, size = 3.5, label.size = 0.1, 
                    box.padding = 0.5, seed = 309,
                    color = c(rep("white", 1), rep("black", 8)),
                    alpha = 0.8
)


obj_pbmc <- subset(obj, subset = Sample_Type == "PBMC")
b <- DimPlot(obj_pbmc, group.by = "CT_new", 
             reduction = "integrated_sct_umap",
             label = F) +
  scale_color_manual(values = cluster_col)  +
  ggtitle("PBMC") + 
  NoLegend() 
pb <- LabelClusters(b, id = "CT_new", box = TRUE, size = 3.5, label.size = 0.1, 
                    box.padding = 0.5, seed = 309,
                    color = c(rep("white", 1), rep("black", 7)),
                    alpha = 0.8
)


# Stacked barplot
p2 <- dittoBarPlot(
  object = obj,
  var = "CT_new",
  group.by = "Sample_Type", 
  color.panel = cluster_col, # color based on cell types
  main = "",
  #x.reorder = c(2,1)
) #+ coord_flip()


# Proportion test
prop_test <- sc_utils(obj)
prop_test <- permutation_test(
  prop_test, 
  cluster_identity = "CT_new",
  sample_1 = "CSF", sample_2 = "PBMC",
  sample_identity = "Sample_Type"
)

prop_test@results$permutation # SAVE THIS FOR MANUSCRIPT (PVALUES)
results <- prop_test@results$permutation

p3 <- permutation_plot(prop_test) + 
  ggtitle("CSF vs PBMC") +
  ylab("") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5, face = "bold",
                                  size = 14)) 

#write.csv(results, "/home/aoill/tables/sc_proportion_test_CSF_vs_PBMC_all.csv",
#          row.names = F, quote = F)


pdf("/home/aoill/plots/Figure_S03.pdf",
    height = 7, width = 9)
ggarrange(ggarrange(pa, pb), 
          ggarrange(p2, p3, 
                    ncol = 2, nrow = 1, 
                    labels = c("B", "C")),
          ncol = 1, nrow = 2,
          labels = c("A"),
          heights = c(1,.8))
dev.off()

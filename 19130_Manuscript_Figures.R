#==============================================================================#
# Author : Angela M. Oill, aoill@tgen.org
# Date: 2026/01/30
# Project: Pediatric CAR-T 
# Description: Manuscript figures relating to scRNA- and scTCR-seq data
#==============================================================================#


#==============================================================================#
# Load libraries ----
#==============================================================================#
library(Matrix)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(scRepertoire)
library(RColorBrewer)
library(tidyr)
library(ComplexHeatmap)
library(viridis)
library(presto)
library(scProportionTest)
library(dittoSeq)

source("/home/aoill/projects/SingleCellBestPractices/scripts/preprocessing_qc_module.R")
source("/home/aoill/projects/SingleCellBestPractices/scripts/helper_functions_module.R")


#==============================================================================#
# Set variables ----
#==============================================================================#
set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)

## Plot colors ----
ct_colors <- c(
  "T (Memory)" = "#88CCEE", "T (GAPDH+ Glycolysis)" = "#332288", 
  "T (Effector-memory)" = "#44AA99", "T (Undifferentiated)" = "#117733",
  "T (Cytotoxic)" = "#AA4499", 
  "T (Heat-shock)" = "#CC6677", "T (Proliferating)" = "#DDCC77", 
  "Treg" = "#882255",
  "T (Effector)" = "#75871B", "T (Naïve)" = "#DAB2FF", #"#E88C23"
  "T (Central memory)" = "#FF637E", #"#A25933",
  "T (MAIT)" = "#3A4330",
  "T (Activated)" = "#FF2056",
  "NK" = "#0000ff", 
  "NK (CD56-bright, CD16+)" = "#0000ff",  
  "NK (CD56-dim, CD16+)" = "#6767D6", "B" = "#994F00", "pDC" = "#d55e00",
  "Mast cells" = "#FDC745",
  
  "Classical Monocyte" = "#bcd095", "Mo/Mac" = "#6dce52", "Non-classical Monocyte" = "#697f3f",
  "S100 high, HLA class II low Monocyte" = "#64e9b3",
  "Intermediate Monocyte" = "#81e858",
  "Macrophage" = "#64A373", 
  "cDC1" = "#927754", "cDC2" = "#e3bb87", 
  "migDC" = "#de9e35"
)


ct_colors_product <- c("Undifferentiated" = "#117733",
                       "Proliferating" = "#DDCC77", 
                       "Effector" = "#75871B",
                       "Activated" = "#FF2056", 
                       "CAR high" = "#A2F4FD")

## Curated marker genes ----
proliferation <- c("MKI67", "TOP2A")
t_cells <- c("CD3E", "CD3D", "CD8A", "CD4", 
             "FOXP3", "IL2RA", "CTLA4", # Treg
             "HOPX", "NCAM1", "CST7", "NKG7" # NKT
)
nk_cells <- c("GNLY", "NKG7", "KLRD1")
pdc <- c("LILRA4", "CLEC4C", "JCHAIN" # pDCs
)
b_cells <- c("MS4A1", "CD79A", "CD19")


t_cell_states <- unique(c(
  "MKI67", "CD3E", "CD3D", "CD8A", "CD4", 
  "FOXP3", "IL2RA", "CTLA4", 
  "HOPX", "NCAM1", "CST7", "NKG7", 
  "LEF1", "CCR7", "SATB1", "KLF2", "SELL", "IL7R",
  "RORA", "CXCR1", "GZMB", "GZMH", "PRF1", "KLRG1", "GNLY", "NKG7",
  "S100A4", "IL7R", "EOMES", "TCF7", "SELL", "CCR7", "LEF1", "PDCD1",
  "CD69", "IFNG",
  "LEF1", "CCR7", "SATB1", "KLF2", "SELL", "IL7R",
  "CXCR1", "GZMB", "GZMH", "PRF1", "KLRG1", "GNLY", "NKG7", 
  "GZMK", "KLRD1",
  "CCL5", "S100A4", "IL7R", "EOMES", "TCF7", "SELL", "CCR7", "LEF1", "PDCD1",
  "CCR7",
  "ITGAE", "ZNF683",
  "HAVCR2", "TOX", "ENTPD1", "LAG3", "PDCD1", "CXCL13", "LAYN"))

macrophages <- c("LYZ", "MARCO", "FCGR1A", "C1QA", "SPP1")
CD14_monocytes <- c("CD14", "LYZ")
FCGR3A_monocytes <- c("FCGR3A", "MS4A7")
dcs <- c("FCER1A", "CST3", # DCs in general
         "CD8A", "CLEC9A", "ITGAE", "ITGAX", "THBD", 
         #"CD141", 
         "XCR1", # cDC1
         "CD1C", "CD207", "ITGAM", "NOTCH2", "SIRPA", # cDC2
         "LAMP3", "CD83", "CCR7" # mDCs
)


#==============================================================================#
# Load objects ----
#==============================================================================#
csf_mye <- readRDS("/scratch/aoill/projects/CAR-T/00_new_2025/geo_2026/csf_myeloid_seurat_obj.rds")
csf_lym <- readRDS("/scratch/aoill/projects/CAR-T/00_new_2025/geo_2026/csf_lymphoid_seurat_obj.rds")
pbmc_mye <- readRDS("/scratch/aoill/projects/CAR-T/00_new_2025/geo_2026/pbmc_myeloid_seurat_obj.rds")
pbmc_lym <- readRDS("/scratch/aoill/projects/CAR-T/00_new_2025/geo_2026/pbmc_lymphoid_seurat_obj.rds")
product_obj <- readRDS("/scratch/aoill/projects/CAR-T/00_new_2025/geo_2026/product_seurat_obj.rds")
tumor_imm <- readRDS("/scratch/aoill/projects/CAR-T/00_new_2025/geo_2026/tumor_immune_seurat_obj.rds")
csf_pbmc <- readRDS("/scratch/aoill/projects/CAR-T/00_new_2025/geo_2026/csf_and_pbmc_integrated_together_seurat_obj.rds")


#==============================================================================#
# Figure 2A ----
#==============================================================================#
# UMAPs of CSF lymphoid and myeloid lineages. See xxx.R for how 
# csf_myeloid_seurat_obj.rds and csf_lymphoid_seurat_obj.rds were generated
a1 <- DimPlot(csf_lym, group.by = "ct_final_all", 
             reduction = "umap.integrated.rpca") +
  ggtitle("CSF - Lymphoid") + NoLegend() + coord_fixed() +
  scale_color_manual(values = ct_colors)
pa1 <- LabelClusters(a1, id = "ct_final_all", box = TRUE, size = 3.5, label.size = 0.5, 
                    box.padding = 0.1, seed = 309,
                    alpha = 0.8,
                    color = c("black", "white", "black", "white", 
                              "white", "white", "white", "black", "black",
                              "black", "black"),
)

pdf("/home/aoill/plots/00_19130_CART/Figure_02_A_left.pdf",
    width = 4, height = 4)
pa1
dev.off()



a2 <- DimPlot(csf_mye, group.by = "ct_final_all", 
             reduction = "umap.integrated.rpca") +
  ggtitle("CSF - Myeloid") + NoLegend() + coord_fixed() +
  scale_color_manual(values = ct_colors)
pa2 <- LabelClusters(a2, id = "ct_final_all", box = TRUE, size = 3.5, label.size = 0.5, 
                    box.padding = 0.1, seed = 309,
                    alpha = 0.8
)


pdf("/home/aoill/plots/00_19130_CART/Figure_02_A_right.pdf",
    width = 4, height = 4)
pa2
dev.off()


#==============================================================================#
# Figure 2B ----
#==============================================================================#
b1 <- DimPlot(pbmc_lym, group.by = "ct_final_all", 
             reduction = "umap.integrated.rpca") +
  ggtitle("PBMC - Lymphoid") + NoLegend() + coord_fixed() +
  scale_color_manual(values = ct_colors)
pb1 <- LabelClusters(b1, id = "ct_final_all", box = TRUE, size = 3.5, label.size = 0.5, 
                    box.padding = 0.1, seed = 309,
                    alpha = 0.8,
                    color = c("black", "black", "black", "black", 
                              "white", "black", "black", "white", "white",
                              "black", "black", "white", "black"),
)

pdf("/home/aoill/plots/00_19130_CART/Figure_02_B_left.pdf",
    width = 4, 
    height = 6)
pb1
dev.off()


b2 <- DimPlot(pbmc_mye, group.by = "ct_final_all", 
             reduction = "umap.integrated.rpca") +
  ggtitle("CSF - Myeloid") + NoLegend() + coord_fixed() +
  scale_color_manual(values = ct_colors)
pb2 <- LabelClusters(b2, id = "ct_final_all", box = TRUE, size = 3.5, label.size = 0.5, 
                    box.padding = 0.1, seed = 309,
                    alpha = 0.8
)


pdf("/home/aoill/plots/00_19130_CART/Figure_02_B_right.pdf",
    width = 5, height = 5)
pb2
dev.off()


#==============================================================================#
# Figure 2C ----
#==============================================================================#
csf_mye <- JoinLayers(csf_mye)
csf_lym <- JoinLayers(csf_lym)
all_obj_list <- list(csf_lym, csf_mye)
names(all_obj_list) <- c("lymphoid", "myeloid")

all_obj_merge <- merge(x = all_obj_list[[1]], y = all_obj_list[2:length(all_obj_list)])

all_obj_merge@meta.data$UPN_Cycle <- factor(x = all_obj_merge@meta.data$UPN_Cycle, 
                                            levels = c("514_2", "514_5", "514_8", "514_11", "514_13",
                                                       "515_4", "515_8",
                                                       "574_1", "574_3", "574_8",  
                                                       "689_1", "689_4", "689_8", "689_12", "689_17", 
                                                       "692_2", "692_4",
                                                       "705_5", "705_8",  
                                                       "716_2", "716_4", "716_8"))


pdf("/home/aoill/plots/00_19130_CART/Figure_02_C.pdf",
    width = 12,
    height = 4)
dittoBarPlot(
  object = all_obj_merge,
  var = "ct_final_all",
  group.by = "UPN_Cycle", 
  color.panel = ct_colors, 
  main = "CSF",
  #x.reorder = c(2,3,1,4,5,6,7)
  x.reorder = c(3, 4, 5, 1, 2, 
                6, 7, 
                8, 9, 10, 
                11, 14, 15, 12, 13, 
                16, 17,
                18, 19, 
                20, 21, 22
  ),
  var.labels.reorder = c(2, 3, 4, 5, 6, 7, 9, 10, 1, 8, 
                         11, 12, 13, 14, 15, 16, 17, 18)
) + 
  geom_vline(xintercept = 5.5, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 7.5, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 10.5, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 15.5, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 17.5, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 19.5, linetype = "dashed", size = 1) +
  #NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  guides(fill = guide_legend(ncol = 2))
dev.off()


#==============================================================================#
# Figure 2D ----
#==============================================================================#
pbmc_mye <- JoinLayers(pbmc_mye)
pbmc_lym <- JoinLayers(pbmc_lym)

pbmc_all_obj_list <- list(pbmc_lym, pbmc_mye)
names(pbmc_all_obj_list) <- c("lymphoid", "myeloid")

pbmc_all_obj_merge <- merge(x = pbmc_all_obj_list[[1]], y = pbmc_all_obj_list[2:length(pbmc_all_obj_list)])

unique(pbmc_all_obj_merge@meta.data$UPN_Cycle)
unique(pbmc_all_obj_merge@meta.data$ct_final_all)


pdf("/home/aoill/plots/00_19130_CART/Figure_02_D.pdf",
    width = 14,
    height = 4)
dittoBarPlot(
  object = pbmc_all_obj_merge,
  var = "ct_final_all",
  group.by = "UPN_Cycle", 
  color.panel = ct_colors, 
  main = "PBMC",
  x.reorder = c(3, 4, 5, 1, 2, 
                6, 7, 
                8, 9, 10, 
                11, 12, 13, 
                14, 15, 16, 
                17, 20, 18, 19,
                21, 22, 
                23, 24, 25, 
                26, 27, 28
  ),
  var.labels.reorder = c(2, 3, 4, 7, 9, 
                         8, 1, 5, 6, 10,
                         11, 12, 13, 14, 15, 16, 17, 18)
) + 
  geom_vline(xintercept = 5.5, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 7.5, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 10.5, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 13.5, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 16.5, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 20.5, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 22.5, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 25.5, linetype = "dashed", size = 1) +
  #NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  guides(fill = guide_legend(ncol = 2))

dev.off()


#==============================================================================#
# Figure 2E ----
#==============================================================================#
prop_test <- sc_utils(all_obj_merge)
prop_test <- permutation_test(
  prop_test, 
  cluster_identity = "ct_final_all",
  sample_1 = "Non-Lymphodepleted", sample_2 = "Lymphodepleted",
  sample_identity = "Lymphodepletion"
)

prop_test@results$permutation # SAVE THIS FOR MANUSCRIPT (PVALUES)
results <- prop_test@results$permutation
write.csv(results, "/scratch/aoill/tables/00_19130_CART/Table_S08.csv",
          row.names = F, quote = F)

pdf("/home/aoill/plots/00_19130_CART/Figure_02_E.pdf",
    width = 6.25,
    height = 4.5)
permutation_plot(prop_test) + 
  ggtitle("CSF\nNon-Lymphodepleted vs Lymphodepleted") +
  ylab("") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5, face = "bold",
                                  size = 10)) 

dev.off()


#==============================================================================#
# Figure 2F ----
#==============================================================================#
prop_test <- sc_utils(pbmc_all_obj_merge)
prop_test <- permutation_test(
  prop_test, 
  cluster_identity = "ct_final_all",
  sample_1 = "Non-Lymphodepleted", sample_2 = "Lymphodepleted",
  sample_identity = "Lymphodepletion"
)

prop_test@results$permutation # SAVE THIS FOR MANUSCRIPT (PVALUES)
results <- prop_test@results$permutation
write.csv(results, "/scratch/aoill/tables/00_19130_CART/Table_S09.csv",
          row.names = F, quote = F)

pdf("/home/aoill/plots/00_19130_CART/Figure_02_F.pdf",
    width = 6.25,
    height = 4.5)
permutation_plot(prop_test) + 
  ggtitle("PBMC\nNon-Lymphodepleted vs Lymphodepleted") +
  ylab("") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5, face = "bold",
                                  size = 10)) 
dev.off()


#==============================================================================#
# Figure 3A ----
#==============================================================================#
pdf("/home/aoill/plots/00_19130_CART/Figure_03_A.pdf",
    width = 4,
    height = 5
)
DimPlot(csf_lym, group.by = "CART", reduction = "umap.integrated.rpca") +
  scale_color_manual(values = c("grey", "red")) + 
  ggtitle("CAR positivity") +
  guides(color=guide_legend(ncol=1, override.aes = list(size=2))) +
  theme(
    legend.position="right",
    #legend.position="bottom",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12)) + coord_fixed()
dev.off()


#==============================================================================#
# Figure 3B ----
#==============================================================================#
# Not generated in R


#==============================================================================#
# Figure 3C ----
#==============================================================================#
csf_lym_T <- subset(csf_lym, 
                    subset = ct_final_T %in% c("Cytotoxic", "Effector-memory",
                                               "GAPDH+ Glycolysis", "Heat-shock",
                                               "Memory", "Proliferating", "Treg",
                                               "Undifferentiated"
                                               )
                    )
prop_test <- sc_utils(csf_lym_T)
prop_test <- permutation_test(
  prop_test, 
  cluster_identity = "ct_final_T",
  sample_1 = "Negative", sample_2 = "Positive",
  sample_identity = "CART"
)

prop_test@results$permutation # SAVE THIS FOR MANUSCRIPT (PVALUES)
results <- prop_test@results$permutation
write.csv(results, "/scratch/aoill/tables/00_19130_CART/Table_S19.csv",
          row.names = F, quote = F)
pc <- permutation_plot(prop_test) + 
  ggtitle("CSF\nCAR- vs CAR+") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5, vjust = 2, face = "bold",
                                  size = 10)) 
pc

pdf("/home/aoill/plots/00_19130_CART/Figure_03_C.pdf",
    width = 5,
    height = 3.5)
print(pc)
dev.off()


#==============================================================================#
# Figure 3D ----
#==============================================================================#
pdf("/home/aoill/plots/00_19130_CART/Figure_03_D.pdf",
    width = 13 , 
    height = 10)
DimPlot(csf_lym, group.by = "Frequency_norm_cat", reduction = "umap.integrated.rpca",
        split.by = "UPN_Cycle",
        ncol = 7
) +
  scale_color_manual(values = c("Most Expanded (Top 1%)" = "#F0F921FF", 
                                "More Expanded (Top 5% - 1%)" = "#FA9E3BFF",  
                                "Expanded (Top 20% - 5%)" = "#D8576BFF", 
                                "Not Expanded" = alpha("#0D0887FF",0.25),
                                "NA" = alpha("grey", 0.25)
  ), 
  na.value=alpha("grey", 0.25)) + 
  #ggtitle("CSF") +
  ggtitle("Normalized TCR Frequency") +
  guides(color=guide_legend(ncol=1, override.aes = list(size=3))) +
  theme(
    legend.position="right",
    #legend.position="bottom",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12)) + coord_fixed()
dev.off()



#==============================================================================#
# Figure 3E ----
#==============================================================================#
# Plot all
plot_data <- csf_lym_T@meta.data

plot_data$UPN_Cycle <- factor(
  plot_data$UPN_Cycle, 
  levels = c("514_2", "514_5", "514_8", "514_11", "514_13", "515_4", "515_8",
             "574_1", "574_3", "574_8", "689_1", "689_4", "689_8", "689_12", "689_17",
             "692_2", "692_4", "705_5", "705_8", "716_2", "716_4", "716_8")
)


custom_colors <- c(
  "Most Expanded (Top 1%)" = "#F0F921FF", 
  "More Expanded (Top 5% - 1%)" = "#FA9E3BFF", 
  "Expanded (Top 20% - 5%)" = "#D8576BFF",
  "Not Expanded" = "grey"
)


levels(plot_data$Frequency_norm_cat)

plot_data_filtered <- plot_data %>% 
  filter(Frequency_norm_cat %in% c("Most Expanded (Top 1%)", 
                                   "More Expanded (Top 5% - 1%)",
                                   "Expanded (Top 20% - 5%)",
                                   "Not Expanded"))


pdf("/home/aoill/plots/00_19130_CART/Figure_03_E.pdf",
    width = 10,
    height = 3.5)
# Plotting
ggplot(plot_data_filtered, aes(x = UPN_Cycle, fill = Frequency_norm_cat)) +
  # Use geom_bar to automatically count the cells and position="fill" for 100% stacked bar
  geom_bar(position = "fill") +
  
  # Apply the custom colors
  scale_fill_manual(values = custom_colors) +
  
  # Add the vertical lines (based on your original dittoBarPlot lines)
  geom_vline(xintercept = 5.5, linetype = "dotted", size = 1) +
  geom_vline(xintercept = 7.5, linetype = "dotted", size = 1) +
  geom_vline(xintercept = 10.5, linetype = "dotted", size = 1) +
  geom_vline(xintercept = 15.5, linetype = "dotted", size = 1) +
  geom_vline(xintercept = 17.5, linetype = "dotted", size = 1) +
  geom_vline(xintercept = 19.5, linetype = "dotted", size = 1) +
  
  # Customize labels and theme
  labs(
    x = "UPN Cycle",
    y = "Proportion of Cells",
    fill = "Clonal Expansion Category" # Legend Title
  ) +
  theme_classic() + # A clean base theme
  theme(
    # Rotate X-axis labels for readability if UPN_Cycle names are long
    axis.text.x = element_text(angle = 45, hjust = 1),
    # Adjust y-axis to show percentages if desired
    axis.title.y = element_text(vjust = 2),
    # Remove plot title since you set main = ""
    plot.title = element_blank()
  ) + coord_cartesian(ylim = c(0.5, 1))

dev.off()


#==============================================================================#
# Figure 3F ----
#==============================================================================#
prop_test <- sc_utils(csf_lym_T)
prop_test <- permutation_test(
  prop_test, 
  cluster_identity = "ct_final_T",
  sample_1 = "Not Expanded", sample_2 = "Expanded",
  sample_identity = "Frequency_norm_cat_2"
)

prop_test@results$permutation # SAVE THIS FOR MANUSCRIPT (PVALUES)
results <- prop_test@results$permutation
write.csv(results, "/scratch/aoill/tables/00_19130_CART/Table_S20.csv",
          row.names = F, quote = F)
pc <- permutation_plot(prop_test) + 
  ggtitle("CSF\nNot Expanded vs Expanded") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5, vjust = 2, face = "bold",
                                  size = 10)) 
pc

pdf("/home/aoill/plots/00_19130_CART/Figure_03_F.pdf",
    width = 5,
    height = 3.5)
print(pc)
dev.off()


#==============================================================================#
# Figure 3G ----
#==============================================================================#
# Pull out metadata from csf object and only keep the necessary TCR info
csf_lym_T_TCR_meta <- csf_lym_T@meta.data
csf_lym_T_TCR_meta_sub <- csf_lym_T_TCR_meta %>%
  dplyr::select(UPN, Sample_Type, Cycle, CTstrict, 
                clonalFrequency, Frequency_norm,
                Frequency_norm_cat, Frequency_norm_cat_2)
rownames(csf_lym_T_TCR_meta_sub) <- NULL
# remove NAs
csf_lym_T_TCR_meta_subNA <- csf_lym_T_TCR_meta_sub[is.na(csf_lym_T_TCR_meta_sub$clonalFrequency) == F, ]

unique(csf_lym_T_TCR_meta_subNA$UPN)


pbmc_lym_T <- subset(pbmc_lym,
                     subset = ct_final_T %in% c("Effector", "Naïve", "MAIT", 
                                                "Effector-memory", "Memory", 
                                                "GAPDH+ Glycolysis", 
                                                "Central memory", "Treg", 
                                                "Proliferating")
)

pbmc_obj_T_meta <- pbmc_lym_T@meta.data
unique(pbmc_obj_T_meta$UPN)
# [1] 514 574 515 626 625 705 689 716 692

# Pull out TCR name only (CTstrict)
pbmc_obj_T_meta_NA_rm <- pbmc_obj_T_meta[is.na(pbmc_obj_T_meta$clonalFrequency) == F, ]
pbmcs_TCRs_514 <- pbmc_obj_T_meta_NA_rm %>%
  filter(UPN == 514) %>%
  pull("CTstrict") %>% unique()

pbmcs_TCRs_515 <- pbmc_obj_T_meta_NA_rm %>%
  filter(UPN == 515) %>%
  pull("CTstrict") %>% unique()

pbmcs_TCRs_574 <- pbmc_obj_T_meta_NA_rm %>%
  filter(UPN == 574) %>%
  pull("CTstrict") %>% unique()

pbmcs_TCRs_626 <- pbmc_obj_T_meta_NA_rm %>%
  filter(UPN == 626) %>%
  pull("CTstrict") %>% unique()

pbmcs_TCRs_625 <- pbmc_obj_T_meta_NA_rm %>%
  filter(UPN == 625) %>%
  pull("CTstrict") %>% unique()

pbmcs_TCRs_705 <- pbmc_obj_T_meta_NA_rm %>%
  filter(UPN == 705) %>%
  pull("CTstrict") %>% unique()

pbmcs_TCRs_689 <- pbmc_obj_T_meta_NA_rm %>%
  filter(UPN == 689) %>%
  pull("CTstrict") %>% unique()

pbmcs_TCRs_716 <- pbmc_obj_T_meta_NA_rm %>%
  filter(UPN == 716) %>%
  pull("CTstrict") %>% unique()

pbmcs_TCRs_692 <- pbmc_obj_T_meta_NA_rm %>%
  filter(UPN == 692) %>%
  pull("CTstrict") %>% unique()



product_obj_meta <- product_obj@meta.data
product_obj_meta_NA_rm <- product_obj_meta[is.na(product_obj_meta$clonalFrequency) == F, ]
product_TCRs_514 <- product_obj_meta_NA_rm %>%
  filter(UPN == 514) %>%
  pull("CTstrict") %>% unique()

product_TCRs_515 <- product_obj_meta_NA_rm %>%
  filter(UPN == 515) %>%
  pull("CTstrict") %>% unique()

product_TCRs_574 <- product_obj_meta_NA_rm %>%
  filter(UPN == 574) %>%
  pull("CTstrict") %>% unique()

product_TCRs_626 <- product_obj_meta_NA_rm %>%
  filter(UPN == 626) %>%
  pull("CTstrict") %>% unique()

product_TCRs_625 <- product_obj_meta_NA_rm %>%
  filter(UPN == 625) %>%
  pull("CTstrict") %>% unique()

# no product for the other samples


unique(csf_lym_T@meta.data$UPN)
#[1] 514 515 574 689 705 716 692

csf_TCRs_exp_514 <- csf_lym_T_TCR_meta_subNA %>% 
  filter(
    Frequency_norm_cat_2 == "Expanded",
    UPN == 514) %>% 
  pull(CTstrict) %>% unique()

csf_TCRs_exp_515 <- csf_lym_T_TCR_meta_subNA %>% 
  filter(
    Frequency_norm_cat_2 == "Expanded",
    UPN == 515) %>% 
  pull(CTstrict) %>% unique()

csf_TCRs_exp_574 <- csf_lym_T_TCR_meta_subNA %>% 
  filter(
    Frequency_norm_cat_2 == "Expanded",
    UPN == 574) %>% 
  pull(CTstrict) %>% unique()

csf_TCRs_exp_705 <- csf_lym_T_TCR_meta_subNA %>% 
  filter(
    Frequency_norm_cat_2 == "Expanded",
    UPN == 705) %>% 
  pull(CTstrict) %>% unique()

csf_TCRs_exp_689 <- csf_lym_T_TCR_meta_subNA %>% 
  filter(
    Frequency_norm_cat_2 == "Expanded",
    UPN == 689) %>% 
  pull(CTstrict) %>% unique()

csf_TCRs_exp_716 <- csf_lym_T_TCR_meta_subNA %>% 
  filter(
    Frequency_norm_cat_2 == "Expanded",
    UPN == 716) %>% 
  pull(CTstrict) %>% unique()

csf_TCRs_exp_692 <- csf_lym_T_TCR_meta_subNA %>% 
  filter(
    Frequency_norm_cat_2 == "Expanded",
    UPN == 692) %>% 
  pull(CTstrict) %>% unique()


# Make a list of the combinations to be plotted in upset plot
tcr_list_514 <- list(Expanded_CSF = csf_TCRs_exp_514, Product = product_TCRs_514, PBMC = pbmcs_TCRs_514)
tcr_list_515 <- list(Expanded_CSF = csf_TCRs_exp_515, Product = product_TCRs_515, PBMC = pbmcs_TCRs_515)
tcr_list_574 <- list(Expanded_CSF = csf_TCRs_exp_574, Product = product_TCRs_574, PBMC = pbmcs_TCRs_574)

# no product for these samples
tcr_list_689 <- list(Expanded_CSF = csf_TCRs_exp_689, PBMC = pbmcs_TCRs_689)
tcr_list_705 <- list(Expanded_CSF = csf_TCRs_exp_705, PBMC = pbmcs_TCRs_705)
tcr_list_716 <- list(Expanded_CSF = csf_TCRs_exp_716, PBMC = pbmcs_TCRs_716)
tcr_list_692 <- list(Expanded_CSF = csf_TCRs_exp_692, PBMC = pbmcs_TCRs_692)


# Reformat data for plotting
m_514 <- make_comb_mat(tcr_list_514)
mat_514 <-list_to_matrix(tcr_list_514)
m_514 <- make_comb_mat(mat_514)

m_515 <- make_comb_mat(tcr_list_515)
mat_515 <-list_to_matrix(tcr_list_515)
m_515 <- make_comb_mat(mat_515)

m_574 <- make_comb_mat(tcr_list_574)
mat_574 <-list_to_matrix(tcr_list_574)
m_574 <- make_comb_mat(mat_574)

m_689 <- make_comb_mat(tcr_list_689)
mat_689 <-list_to_matrix(tcr_list_689)
m_689 <- make_comb_mat(mat_689)

m_705 <- make_comb_mat(tcr_list_705)
mat_705 <-list_to_matrix(tcr_list_705)
m_705 <- make_comb_mat(mat_705)

m_716 <- make_comb_mat(tcr_list_716)
mat_716 <-list_to_matrix(tcr_list_716)
m_716 <- make_comb_mat(mat_716)

m_692 <- make_comb_mat(tcr_list_692)
mat_692 <-list_to_matrix(tcr_list_692)
m_692 <- make_comb_mat(mat_692)


# Make plots
ptest1 <- UpSet(m_514, column_title = expression(bold("UPN 514")), 
                comb_order = order(comb_degree(m_514), -comb_size(m_514)),
                left_annotation = upset_left_annotation(m_514,
                                                        width = unit(1, "cm")),
                top_annotation = upset_top_annotation(m_514, add_numbers = TRUE,
                                                      annotation_name_rot = 90),
                show_row_names = F,
                set_order = c("Product", "PBMC", "Expanded_CSF")
)
ptest1
p1_ann <- grid.grab()

ptest2 <- UpSet(m_515, column_title = expression(bold("UPN 515")), 
                comb_order = order(comb_degree(m_515), -comb_size(m_515)),
                left_annotation = upset_left_annotation(m_515,
                                                        width = unit(1, "cm")),
                top_annotation = upset_top_annotation(m_515, add_numbers = TRUE,
                                                      annotation_name_rot = 90),
                show_row_names = F,
                set_order = c("Product", "PBMC", "Expanded_CSF")
)
ptest2
p2_ann <- grid.grab()

ptest3 <- UpSet(m_574, column_title = expression(bold("UPN 574")), 
                #comb_order = order(comb_size(m_574)),
                comb_order = order(comb_degree(m_574), -comb_size(m_574)),
                left_annotation = upset_left_annotation(m_574,
                                                        width = unit(1, "cm")),
                top_annotation = upset_top_annotation(m_574, add_numbers = TRUE,
                                                      annotation_name_rot = 90),
                set_order = c("Product", "PBMC", "Expanded_CSF")
                
                #show_row_names = FALSE
)
ptest3
p3_ann <- grid.grab()


# No product in these sample so plot needs to reflect this (remove from upset plot)
ptest4 <- UpSet(m_689, column_title = expression(bold("UPN 689")), 
                comb_order = order(comb_degree(m_689), -comb_size(m_689)),
                left_annotation = upset_left_annotation(m_689,
                                                        width = unit(1, "cm")),
                top_annotation = upset_top_annotation(m_689, add_numbers = TRUE,
                                                      annotation_name_rot = 90),
                set_order = c("PBMC", "Expanded_CSF"),
                show_row_names = FALSE
)
ptest4
p4_ann <- grid.grab()


ptest5 <- UpSet(m_705, column_title = expression(bold("UPN 705")), 
                comb_order = order(comb_degree(m_705), -comb_size(m_705)),
                left_annotation = upset_left_annotation(m_705,
                                                        width = unit(1, "cm")),
                top_annotation = upset_top_annotation(m_705, add_numbers = TRUE,
                                                      annotation_name_rot = 90),
                set_order = c("PBMC", "Expanded_CSF"),
                show_row_names = T
)
ptest5
p5_ann <- grid.grab()


#[1] 514 515 574 689 705 716 692
ptest6 <- UpSet(m_716, column_title = expression(bold("UPN 716")), 
                comb_order = order(comb_degree(m_716), -comb_size(m_716)),
                left_annotation = upset_left_annotation(m_716,
                                                        width = unit(1, "cm")),
                top_annotation = upset_top_annotation(m_716, add_numbers = TRUE,
                                                      annotation_name_rot = 90),
                set_order = c("PBMC", "Expanded_CSF"),
                show_row_names = FALSE
)
ptest6
p6_ann <- grid.grab()


ptest7 <- UpSet(m_692, column_title = expression(bold("UPN 692")), 
                comb_order = order(comb_degree(m_692), -comb_size(m_692)),
                left_annotation = upset_left_annotation(m_692,
                                                        width = unit(1, "cm")),
                top_annotation = upset_top_annotation(m_692, add_numbers = TRUE,
                                                      annotation_name_rot = 90),
                set_order = c("PBMC", "Expanded_CSF"),
                show_row_names = T
)
ptest7
p7_ann <- grid.grab()



pdf("/home/aoill/plots/00_19130_CART/Figure_03_G.pdf",
    width = 10,
    height = 3)
ggarrange(
  p1_ann, p2_ann, p3_ann, 
  nrow = 1,
  ncol = 3,
  align = "hv",
  labels = c(""),
  #widths = c(.8,.8,1)), 
  widths = c(.77,.77,1))
dev.off()


#==============================================================================#
# Extended Data Figure 1A ----
#==============================================================================#
all_obj_merge_Ependymoma <- subset(all_obj_merge, subset = tumor_type == "Ependymoma")

prop_test <- sc_utils(all_obj_merge_Ependymoma)
prop_test <- permutation_test(
  prop_test, 
  cluster_identity = "ct_final_all",
  sample_1 = "Non-Lymphodepleted", sample_2 = "Lymphodepleted",
  sample_identity = "Lymphodepletion"
)

prop_test@results$permutation # SAVE THIS FOR MANUSCRIPT (PVALUES)
results <- prop_test@results$permutation
write.csv(results, "/scratch/aoill/tables/00_19130_CART/Table_S10.csv",
          row.names = F, quote = F)

pdf("/home/aoill/plots/00_19130_CART/Extended_Data_Figure_01_A.pdf",
    width = 6.5, height = 4.5)
permutation_plot(prop_test) + 
  ggtitle("CSF - Ependymoma\nNon-Lymphodepleted vs Lymphodepleted") +
  ylab("") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5, face = "bold",
                                  size = 10)) 
dev.off()


#==============================================================================#
# Extended Data Figure 1B ----
#==============================================================================#
all_obj_merge_DMG <- subset(all_obj_merge, subset = tumor_type == "DMG")

prop_test <- sc_utils(all_obj_merge_DMG)
prop_test <- permutation_test(
  prop_test, 
  cluster_identity = "ct_final_all",
  sample_1 = "Non-Lymphodepleted", sample_2 = "Lymphodepleted",
  sample_identity = "Lymphodepletion"
)

prop_test@results$permutation # SAVE THIS FOR MANUSCRIPT (PVALUES)
results <- prop_test@results$permutation
write.csv(results, "/scratch/aoill/tables/00_19130_CART/Table_S11.csv",
          row.names = F, quote = F)

pdf("/home/aoill/plots/00_19130_CART/Extended_Data_Figure_01_B.pdf",
    width = 6.5, height = 4.5)
permutation_plot(prop_test) + 
  ggtitle("CSF - DMG\nNon-Lymphodepleted vs Lymphodepleted") +
  ylab("") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5, face = "bold",
                                  size = 10)) 
dev.off()


#==============================================================================#
# Extended Data Figure 1C ----
#==============================================================================#
all_obj_merge_res_non <- subset(all_obj_merge, subset = Response_cat %in% c("Response", "Non-response"))
table(all_obj_merge_res_non$Response_cat)

all_obj_merge_res_non_LD <- subset(all_obj_merge_res_non, subset = Lymphodepletion == "Lymphodepleted")
prop_test <- sc_utils(all_obj_merge_res_non_LD)
prop_test <- permutation_test(
  prop_test, 
  cluster_identity = "ct_final_all",
  sample_1 = "Non-response", sample_2 = "Response",
  sample_identity = "Response_cat"
)

prop_test@results$permutation # SAVE THIS FOR MANUSCRIPT (PVALUES)
results <- prop_test@results$permutation
write.csv(results, "/scratch/aoill/tables/00_19130_CART/Table_S12.csv",
          row.names = F, quote = F)

pdf("/home/aoill/plots/00_19130_CART/Extended_Data_Figure_01_C_left.pdf",
    width = 6.5, height = 4.5)
permutation_plot(prop_test) + 
  ggtitle("CSF Lymphodepleted\nNon-Response vs Response") +
  ylab("") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5, face = "bold",
                                  size = 10)) 

dev.off()


all_obj_merge_res_non_nonLD <- subset(all_obj_merge_res_non, subset = Lymphodepletion == "Non-Lymphodepleted")
prop_test <- sc_utils(all_obj_merge_res_non_nonLD)
prop_test <- permutation_test(
  prop_test, 
  cluster_identity = "ct_final_all",
  sample_1 = "Non-response", sample_2 = "Response",
  sample_identity = "Response_cat"
)

prop_test@results$permutation # SAVE THIS FOR MANUSCRIPT (PVALUES)
results <- prop_test@results$permutation
write.csv(results, "/scratch/aoill/tables/00_19130_CART/Table_S13.csv",
          row.names = F, quote = F)

pdf("/home/aoill/plots/00_19130_CART/Extended_Data_Figure_01_C_right.pdf",
    width = 6.5, height = 4.5)
permutation_plot(prop_test) + 
  ggtitle("CSF Non-lymphodepleted\nNon-Response vs Response") +
  ylab("") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5, face = "bold",
                                  size = 10)) 
dev.off()



#==============================================================================#
# Extended Data Figure 1D ----
#==============================================================================#
pbmc_all_obj_merge_res_non <- subset(pbmc_all_obj_merge, subset = Response_cat %in% c("Response", "Non-response"))
table(pbmc_all_obj_merge_res_non$Response_cat)

pbmc_all_obj_merge_res_non_LD <- subset(pbmc_all_obj_merge_res_non, subset = Lymphodepletion == "Lymphodepleted")
prop_test <- sc_utils(pbmc_all_obj_merge_res_non_LD)
prop_test <- permutation_test(
  prop_test, 
  cluster_identity = "ct_final_all",
  sample_1 = "Non-response", sample_2 = "Response",
  sample_identity = "Response_cat"
)

prop_test@results$permutation # SAVE THIS FOR MANUSCRIPT (PVALUES)
results <- prop_test@results$permutation
write.csv(results, "/scratch/aoill/tables/00_19130_CART/Table_S14.csv",
          row.names = F, quote = F)

pdf("/home/aoill/plots/00_19130_CART/Extended_Data_Figure_01_D_left.pdf",
    width = 6.5, height = 4.5)
permutation_plot(prop_test) + 
  ggtitle("PBMC Lymphodepleted\nNon-Response vs Response") +
  ylab("") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5, face = "bold",
                                  size = 10)) 

dev.off()


pbmc_all_obj_merge_res_non_nonLD <- subset(pbmc_all_obj_merge_res_non, subset = Lymphodepletion == "Non-Lymphodepleted")
unique(pbmc_all_obj_merge_res_non_nonLD$Lymphodepletion)
unique(pbmc_all_obj_merge_res_non_nonLD$Response_cat)

pbmc_all_obj_merge_res_non_nonLD <- JoinLayers(pbmc_all_obj_merge_res_non_nonLD)

prop_test <- NULL
prop_test <- sc_utils(pbmc_all_obj_merge_res_non_nonLD)
prop_test <- permutation_test(
  prop_test, 
  cluster_identity = "ct_final_all",
  sample_1 = "Non-response", sample_2 = "Response",
  sample_identity = "Response_cat"
)

prop_test@results$permutation # SAVE THIS FOR MANUSCRIPT (PVALUES)
results <- prop_test@results$permutation
write.csv(results, "/scratch/aoill/tables/00_19130_CART/Table_S15.csv",
          row.names = F, quote = F)

pdf("/home/aoill/plots/00_19130_CART/Extended_Data_Figure_01_D_right.pdf",
    width = 6.5, height = 4.5)
permutation_plot(prop_test) + 
  ggtitle("PBMC Non-lymphodepleted\nNon-Response vs Response") +
  ylab("") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5, face = "bold",
                                  size = 10)) 
dev.off()


#==============================================================================#
# Figure S3 ----
#==============================================================================#
DefaultAssay(csf_lym) <- "SoupX_RNA"
csf_lym_join <- JoinLayers(csf_lym)
markers <- presto::wilcoxauc(csf_lym_join, group_by = "ct_final_all", 
                             assay = "data", seurat_assay = "SoupX_RNA")
tp_markers <- top_markers(markers, n = 5, auc_min = 0.6, pct_in_min = 50, 
                          pct_out_max = 100)

tp_markers

all_tp_markers <- tp_markers %>%
  select(-rank) %>% 
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>%
  .[!is.na(.)]

length(all_tp_markers)

all_tp_markers_all <- unique(c(all_tp_markers, t_cell_states, proliferation, 
                               t_cells, nk_cells, pdc, b_cells))


length(all_tp_markers_all)


csf_obj_cp <- csf_lym_join
DefaultAssay(csf_obj_cp) <- "SoupX_RNA" 
csf_obj_cp <- NormalizeData(csf_obj_cp)
VariableFeatures(csf_obj_cp) <- rownames(csf_obj_cp)
csf_obj_cp <- ScaleData(csf_obj_cp)

p <- DotPlot(csf_obj_cp, features = all_tp_markers_all, 
             group.by = "ct_final_all", scale = TRUE) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
p 


# Reproducing the dotplot using ComplexHeatmap
df <- p$data
head(df)
dim(df)

#rownames(df)

# Adding log-transformed values
#df$log.avg.exp <- log10(df$avg.exp)

# (Scaled) expression levels
exp_mat <- df %>% 
  dplyr::select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame()

row.names(exp_mat) <- exp_mat$features.plot
exp_mat <- exp_mat[,-1] %>% as.matrix()

# The percentage of cells expressing a feature
percent_mat <- df %>% 
  dplyr::select(-avg.exp, -avg.exp.scaled) %>%  
  pivot_wider(names_from = id, values_from = pct.exp) %>% 
  as.data.frame() 

row.names(percent_mat) <- percent_mat$features.plot  
percent_mat <- percent_mat[,-1] %>% as.matrix()

dim(exp_mat)
dim(percent_mat)

# Get an idea of the ranges of the matrix
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99), na.rm=T)

# relace NA with 0 in matrix
exp_mat[!is.finite(exp_mat)] <- 0
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99))

# Any value that is greater than 2 will be mapped to yellow
col_fun = circlize::colorRamp2(c(-1, 0, 2), viridis(20)[c(1,10, 20)])

# Adding annotation colors
#library(scales) 
#hue_pal()(12) # hex <- hue_pal()(12) # there are 12 clusters
#ct_colors_sub_heatmap <- c(list(T = ct_colors_sub))


ct_colors_sub_heatmap <- c(list(Lymphoid = ct_colors))



#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
cluster_anno <- colnames(exp_mat)
column_ha <- HeatmapAnnotation( # HeatmapAnnotation / rowAnnotation
  Lymphoid = cluster_anno,
  col = ct_colors_sub_heatmap, 
  #col = col, 
  na_col = "grey"
)

layer_fun = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.points(x = x, y = y, 
              gp = gpar(col = col_fun(pindex((exp_mat), i, j))), # t(exp_mat)
              size = pindex((percent_mat), i, j)/100 * unit(4, "mm"), # t(percent_mat)
              pch = 16
  )
}

lgd_list = list(
  Legend( labels = c(0,0.25,0.5,0.75,1), title = "pt", type = "points", 
          pch = 16, size = c(0,0.25,0.5,0.75,1) * unit(4,"mm"),
          legend_gp = gpar(col = "black")))

hp <- ComplexHeatmap::Heatmap((exp_mat), # t(exp_mat)
                              heatmap_legend_param=list(title="Scaled expression"),
                              column_title = "Lymphoid", 
                              col = col_fun,
                              rect_gp = gpar(type = "none"),
                              layer_fun = layer_fun,
                              row_names_gp = gpar(fontsize = 8),
                              column_names_gp = gpar(fontsize = 8),
                              #row_km = 4,
                              column_km = 4,
                              top_annotation = column_ha, # top_annotation / left_annotation
                              border = "black")


draw(hp, annotation_legend_list = lgd_list) 


pdf("/home/aoill/plots/00_19130_CART/Figure_S3.pdf",
    width = 7.5, height = 15)
draw(hp, annotation_legend_list = lgd_list) 

dev.off()


#==============================================================================#
# Figure S4 ----
#==============================================================================#
DefaultAssay(csf_mye) <- "SoupX_RNA"
csf_mye_join <- JoinLayers(csf_mye)
markers <- presto::wilcoxauc(csf_mye_join, group_by = "ct_final_all", 
                             assay = "data", seurat_assay = "SoupX_RNA")
tp_markers <- top_markers(markers, n = 10, auc_min = 0.6, pct_in_min = 50, 
                          pct_out_max = 100)

tp_markers

all_tp_markers <- tp_markers %>%
  select(-rank) %>% 
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>%
  .[!is.na(.)]

length(all_tp_markers)

all_tp_markers_all <- unique(c(all_tp_markers, macrophages, 
                               CD14_monocytes, FCGR3A_monocytes, dcs))


length(all_tp_markers_all)


csf_obj_cp <- csf_mye_join
DefaultAssay(csf_obj_cp) <- "SoupX_RNA" 
csf_obj_cp <- NormalizeData(csf_obj_cp)
VariableFeatures(csf_obj_cp) <- rownames(csf_obj_cp)
csf_obj_cp <- ScaleData(csf_obj_cp)

p <- DotPlot(csf_obj_cp, features = all_tp_markers_all, 
             group.by = "ct_final_all", scale = TRUE) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
p 


# Reproducing the dotplot using ComplexHeatmap
df <- p$data
head(df)
dim(df)

#rownames(df)

# Adding log-transformed values
#df$log.avg.exp <- log10(df$avg.exp)

# (Scaled) expression levels
exp_mat <- df %>% 
  dplyr::select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame()

row.names(exp_mat) <- exp_mat$features.plot
exp_mat <- exp_mat[,-1] %>% as.matrix()

# The percentage of cells expressing a feature
percent_mat <- df %>% 
  dplyr::select(-avg.exp, -avg.exp.scaled) %>%  
  pivot_wider(names_from = id, values_from = pct.exp) %>% 
  as.data.frame() 

row.names(percent_mat) <- percent_mat$features.plot  
percent_mat <- percent_mat[,-1] %>% as.matrix()

dim(exp_mat)
dim(percent_mat)

# Get an idea of the ranges of the matrix
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99), na.rm=T)

# relace NA with 0 in matrix
exp_mat[!is.finite(exp_mat)] <- 0
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99))

# Any value that is greater than 2 will be mapped to yellow
col_fun = circlize::colorRamp2(c(-1, 0, 2), viridis(20)[c(1,10, 20)])

# Adding annotation colors
#library(scales) 
#hue_pal()(12) # hex <- hue_pal()(12) # there are 12 clusters
#ct_colors_sub_heatmap <- c(list(Myeloid = ct_colors_sub))

ct_colors_sub_heatmap <- c(list(Myeloid = ct_colors))



#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
cluster_anno <- colnames(exp_mat)
column_ha <- HeatmapAnnotation( # HeatmapAnnotation / rowAnnotation
  Myeloid = cluster_anno,
  col = ct_colors_sub_heatmap, 
  #col = col, 
  na_col = "grey"
)

layer_fun = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.points(x = x, y = y, 
              gp = gpar(col = col_fun(pindex((exp_mat), i, j))), # t(exp_mat)
              size = pindex((percent_mat), i, j)/100 * unit(4, "mm"), # t(percent_mat)
              pch = 16
  )
}

lgd_list = list(
  Legend( labels = c(0,0.25,0.5,0.75,1), title = "pt", type = "points", 
          pch = 16, size = c(0,0.25,0.5,0.75,1) * unit(4,"mm"),
          legend_gp = gpar(col = "black")))

hp <- ComplexHeatmap::Heatmap((exp_mat), # t(exp_mat)
                              heatmap_legend_param=list(title="Scaled expression"),
                              column_title = "Myeloid", 
                              col = col_fun,
                              rect_gp = gpar(type = "none"),
                              layer_fun = layer_fun,
                              row_names_gp = gpar(fontsize = 8),
                              column_names_gp = gpar(fontsize = 8),
                              #row_km = 4,
                              column_km = 5,
                              top_annotation = column_ha, # top_annotation / left_annotation
                              border = "black")


draw(hp, annotation_legend_list = lgd_list) 


pdf("/home/aoill/plots/00_19130_CART/Figure_S4.pdf",
    width = 7, height = 11)
draw(hp, annotation_legend_list = lgd_list) 


dev.off()


#==============================================================================#
# Figure S5 ----
#==============================================================================#
DefaultAssay(pbmc_lym) <- "SoupX_RNA"
pbmc_lym_join <- JoinLayers(pbmc_lym)
markers <- presto::wilcoxauc(pbmc_lym_join, group_by = "ct_final_all", 
                             assay = "data", seurat_assay = "SoupX_RNA")
tp_markers <- top_markers(markers, n = 10, auc_min = 0.6, pct_in_min = 50, 
                          pct_out_max = 100)

tp_markers

as.data.frame(top_markers(markers, n = 25, auc_min = 0.6, pct_in_min = 50, 
                          pct_out_max = 100))


all_tp_markers <- tp_markers %>%
  select(-rank) %>% 
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>%
  .[!is.na(.)]

length(all_tp_markers)


all_tp_markers_all <- unique(c(all_tp_markers, t_cell_states, proliferation, 
                               t_cells, nk_cells, pdc, b_cells))


length(all_tp_markers_all)


pbmc_obj_cp <- pbmc_lym_join
DefaultAssay(pbmc_obj_cp) <- "SoupX_RNA" 
pbmc_obj_cp <- NormalizeData(pbmc_obj_cp)
VariableFeatures(pbmc_obj_cp) <- rownames(pbmc_obj_cp)
pbmc_obj_cp <- ScaleData(pbmc_obj_cp)

p <- DotPlot(pbmc_obj_cp, features = all_tp_markers_all, 
             group.by = "ct_final_all", scale = TRUE) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
p 


# Reproducing the dotplot using ComplexHeatmap
df <- p$data
head(df)
dim(df)

#rownames(df)

# Adding log-transformed values
#df$log.avg.exp <- log10(df$avg.exp)

# (Scaled) expression levels
exp_mat <- df %>% 
  dplyr::select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame()

row.names(exp_mat) <- exp_mat$features.plot
exp_mat <- exp_mat[,-1] %>% as.matrix()

# The percentage of cells expressing a feature
percent_mat <- df %>% 
  dplyr::select(-avg.exp, -avg.exp.scaled) %>%  
  pivot_wider(names_from = id, values_from = pct.exp) %>% 
  as.data.frame() 

row.names(percent_mat) <- percent_mat$features.plot  
percent_mat <- percent_mat[,-1] %>% as.matrix()

dim(exp_mat)
dim(percent_mat)

# Get an idea of the ranges of the matrix
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99), na.rm=T)

# relace NA with 0 in matrix
exp_mat[!is.finite(exp_mat)] <- 0
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99))

# Any value that is greater than 2 will be mapped to yellow
col_fun = circlize::colorRamp2(c(-1, 0, 2), viridis(20)[c(1,10, 20)])

# Adding annotation colors
#library(scales) 
#hue_pal()(12) # hex <- hue_pal()(12) # there are 12 clusters
#ct_colors_sub_heatmap <- c(list(T = ct_colors_sub))


ct_colors_sub_heatmap <- c(list(Lymphoid = ct_colors))



#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
cluster_anno <- colnames(exp_mat)
column_ha <- HeatmapAnnotation( # HeatmapAnnotation / rowAnnotation
  Lymphoid = cluster_anno,
  col = ct_colors_sub_heatmap, 
  #col = col, 
  na_col = "grey"
)

layer_fun = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.points(x = x, y = y, 
              gp = gpar(col = col_fun(pindex((exp_mat), i, j))), # t(exp_mat)
              size = pindex((percent_mat), i, j)/100 * unit(4, "mm"), # t(percent_mat)
              pch = 16
  )
}

lgd_list = list(
  Legend( labels = c(0,0.25,0.5,0.75,1), title = "pt", type = "points", 
          pch = 16, size = c(0,0.25,0.5,0.75,1) * unit(4,"mm"),
          legend_gp = gpar(col = "black")))

hp <- ComplexHeatmap::Heatmap((exp_mat), # t(exp_mat)
                              heatmap_legend_param=list(title="Scaled expression"),
                              column_title = "Lymphoid", 
                              col = col_fun,
                              rect_gp = gpar(type = "none"),
                              layer_fun = layer_fun,
                              row_names_gp = gpar(fontsize = 8),
                              column_names_gp = gpar(fontsize = 8),
                              #row_km = 4,
                              column_km = 4,
                              top_annotation = column_ha, # top_annotation / left_annotation
                              border = "black")


draw(hp, annotation_legend_list = lgd_list) 


pdf("/home/aoill/plots/00_19130_CART/Figure_S5.pdf",
    width = 7.5, height = 16)
draw(hp, annotation_legend_list = lgd_list) 

dev.off()


#==============================================================================#
# Figure S6 ----
#==============================================================================#
DefaultAssay(pbmc_mye) <- "SoupX_RNA"
pbmc_mye_join <- JoinLayers(pbmc_mye)
markers <- presto::wilcoxauc(pbmc_mye_join, group_by = "ct_final_all", 
                             assay = "data", seurat_assay = "SoupX_RNA")
tp_markers <- top_markers(markers, n = 10, auc_min = 0.6, pct_in_min = 50, 
                          pct_out_max = 100)

tp_markers

all_tp_markers <- tp_markers %>%
  select(-rank) %>% 
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>%
  .[!is.na(.)]

length(all_tp_markers)

# not plotting mDCs because they are taking up space
dcs <- c("FCER1A", "CST3", # DCs in general
         #"CD8A", "CLEC9A", "ITGAE", "ITGAX", "THBD", 
         #"CD141", 
         #"XCR1", # cDC1
         "CD1C", "CD207", "ITGAM", "NOTCH2", "SIRPA"#, # cDC2
         #"LAMP3", "CD83", "CCR7" # mDCs
)


all_tp_markers_all <- unique(c(all_tp_markers, #macrophages, 
                               CD14_monocytes, FCGR3A_monocytes, dcs))


length(all_tp_markers_all)


pbmc_obj_cp <- pbmc_mye_join
DefaultAssay(pbmc_obj_cp) <- "SoupX_RNA" 
pbmc_obj_cp <- NormalizeData(pbmc_obj_cp)
VariableFeatures(pbmc_obj_cp) <- rownames(pbmc_obj_cp)
pbmc_obj_cp <- ScaleData(pbmc_obj_cp)

p <- DotPlot(pbmc_obj_cp, features = all_tp_markers_all, 
             group.by = "ct_final_all", scale = TRUE) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
p 


# Reproducing the dotplot using ComplexHeatmap
df <- p$data
head(df)
dim(df)

#rownames(df)

# Adding log-transformed values
#df$log.avg.exp <- log10(df$avg.exp)

# (Scaled) expression levels
exp_mat <- df %>% 
  dplyr::select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame()

row.names(exp_mat) <- exp_mat$features.plot
exp_mat <- exp_mat[,-1] %>% as.matrix()

# The percentage of cells expressing a feature
percent_mat <- df %>% 
  dplyr::select(-avg.exp, -avg.exp.scaled) %>%  
  pivot_wider(names_from = id, values_from = pct.exp) %>% 
  as.data.frame() 

row.names(percent_mat) <- percent_mat$features.plot  
percent_mat <- percent_mat[,-1] %>% as.matrix()

dim(exp_mat)
dim(percent_mat)

# Get an idea of the ranges of the matrix
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99), na.rm=T)

# relace NA with 0 in matrix
exp_mat[!is.finite(exp_mat)] <- 0
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99))

# Any value that is greater than 2 will be mapped to yellow
col_fun = circlize::colorRamp2(c(-1, 0, 2), viridis(20)[c(1,10, 20)])

# Adding annotation colors
#library(scales) 
#hue_pal()(12) # hex <- hue_pal()(12) # there are 12 clusters
#ct_colors_sub_heatmap <- c(list(Myeloid = ct_colors_sub))

ct_colors_sub_heatmap <- c(list(Myeloid = ct_colors))



#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
cluster_anno <- colnames(exp_mat)
column_ha <- HeatmapAnnotation( # HeatmapAnnotation / rowAnnotation
  Myeloid = cluster_anno,
  col = ct_colors_sub_heatmap, 
  #col = col, 
  na_col = "grey"
)

layer_fun = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.points(x = x, y = y, 
              gp = gpar(col = col_fun(pindex((exp_mat), i, j))), # t(exp_mat)
              size = pindex((percent_mat), i, j)/100 * unit(4, "mm"), # t(percent_mat)
              pch = 16
  )
}

lgd_list = list(
  Legend( labels = c(0,0.25,0.5,0.75,1), title = "pt", type = "points", 
          pch = 16, size = c(0,0.25,0.5,0.75,1) * unit(4,"mm"),
          legend_gp = gpar(col = "black")))

hp <- ComplexHeatmap::Heatmap((exp_mat), # t(exp_mat)
                              heatmap_legend_param=list(title="Scaled expression"),
                              column_title = "Myeloid", 
                              col = col_fun,
                              rect_gp = gpar(type = "none"),
                              layer_fun = layer_fun,
                              row_names_gp = gpar(fontsize = 8),
                              column_names_gp = gpar(fontsize = 8),
                              #row_km = 4,
                              column_km = 4,
                              top_annotation = column_ha, # top_annotation / left_annotation
                              border = "black")


draw(hp, annotation_legend_list = lgd_list) 


pdf("/home/aoill/plots/00_19130_CART/Figure_S6.pdf",
    width = 7, height = 11)
draw(hp, annotation_legend_list = lgd_list) 


dev.off()


dcs <- c("FCER1A", "CST3", # DCs in general
         #"CD8A", "CLEC9A", "ITGAE", "ITGAX", "THBD", 
         #"CD141", 
         #"XCR1", # cDC1
         "CD1C", "CD207", "ITGAM", "NOTCH2", "SIRPA", # cDC2
         "LAMP3", "CD83", "CCR7" # mDCs
)


#==============================================================================#
# Figure S7A ----
#==============================================================================#
a <- DimPlot(product_obj, group.by = "ct_final_all", 
             reduction = "umap.integrated.rpca") +
  ggtitle("Product") + NoLegend() + coord_fixed() +
  scale_color_manual(values = ct_colors_product)
pa <- LabelClusters(a, id = "ct_final_all", box = TRUE, size = 3.5, label.size = 0.5, 
                    box.padding = 0.1, seed = 309,
                    alpha = 0.8,
                    color = c("black", "white", "black", "black", 
                              "black"),
)
pa


pdf("/home/aoill/plots/00_19130_CART/Figure_S7_A.pdf",
    width = 5, height = 5)
pa
dev.off()


#==============================================================================#
# Figure S7B ----
#==============================================================================#
DefaultAssay(product_obj) <- "SoupX_RNA"
product_obj <- JoinLayers(product_obj)
markers <- presto::wilcoxauc(product_obj, group_by = "ct_final_all", 
                             assay = "data", seurat_assay = "SoupX_RNA")
tp_markers <- top_markers(markers, n = 10, auc_min = 0.6, pct_in_min = 50, 
                          pct_out_max = 100)

tp_markers


all_tp_markers <- tp_markers %>%
  select(-rank) %>% 
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>%
  .[!is.na(.)]

length(all_tp_markers)


all_tp_markers_all <- unique(c(all_tp_markers, t_cell_states))


length(all_tp_markers_all)


product_obj_cp <- product_obj
DefaultAssay(product_obj_cp) <- "SoupX_RNA" 
product_obj_cp <- NormalizeData(product_obj_cp)
VariableFeatures(product_obj_cp) <- rownames(product_obj_cp)
product_obj_cp <- ScaleData(product_obj_cp)

p <- DotPlot(product_obj_cp, features = all_tp_markers_all, 
             group.by = "ct_final_all", scale = TRUE) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
p 


# Reproducing the dotplot using ComplexHeatmap
df <- p$data
head(df)
dim(df)

#rownames(df)

# Adding log-transformed values
#df$log.avg.exp <- log10(df$avg.exp)

# (Scaled) expression levels
exp_mat <- df %>% 
  dplyr::select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame()

row.names(exp_mat) <- exp_mat$features.plot
exp_mat <- exp_mat[,-1] %>% as.matrix()

# The percentage of cells expressing a feature
percent_mat <- df %>% 
  dplyr::select(-avg.exp, -avg.exp.scaled) %>%  
  pivot_wider(names_from = id, values_from = pct.exp) %>% 
  as.data.frame() 

row.names(percent_mat) <- percent_mat$features.plot  
percent_mat <- percent_mat[,-1] %>% as.matrix()

dim(exp_mat)
dim(percent_mat)

# Get an idea of the ranges of the matrix
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99), na.rm=T)

# relace NA with 0 in matrix
exp_mat[!is.finite(exp_mat)] <- 0
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99))

# Any value that is greater than 2 will be mapped to yellow
col_fun = circlize::colorRamp2(c(-1, 0, 2), viridis(20)[c(1,10, 20)])

# Adding annotation colors
#library(scales) 
#hue_pal()(12) # hex <- hue_pal()(12) # there are 12 clusters
#ct_colors_product_sub_heatmap <- c(list(T = ct_colors_product_sub))


ct_colors_product_sub_heatmap <- c(list(Product = ct_colors_product))



#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
cluster_anno <- colnames(exp_mat)
column_ha <- HeatmapAnnotation( # HeatmapAnnotation / rowAnnotation
  Product = cluster_anno,
  col = ct_colors_product_sub_heatmap, 
  #col = col, 
  na_col = "grey"
)

layer_fun = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.points(x = x, y = y, 
              gp = gpar(col = col_fun(pindex((exp_mat), i, j))), # t(exp_mat)
              size = pindex((percent_mat), i, j)/100 * unit(4, "mm"), # t(percent_mat)
              pch = 16
  )
}

lgd_list = list(
  Legend( labels = c(0,0.25,0.5,0.75,1), title = "pt", type = "points", 
          pch = 16, size = c(0,0.25,0.5,0.75,1) * unit(4,"mm"),
          legend_gp = gpar(col = "black")))

hp <- ComplexHeatmap::Heatmap((exp_mat), # t(exp_mat)
                              heatmap_legend_param=list(title="Scaled expression"),
                              column_title = "T cells", 
                              col = col_fun,
                              rect_gp = gpar(type = "none"),
                              layer_fun = layer_fun,
                              row_names_gp = gpar(fontsize = 8),
                              column_names_gp = gpar(fontsize = 8),
                              #row_km = 4,
                              column_km = 4,
                              top_annotation = column_ha, # top_annotation / left_annotation
                              border = "black")


draw(hp, annotation_legend_list = lgd_list) 


pdf("/home/aoill/plots/Figure_S7_B.pdf",
    width = 7.5, height = 13)
draw(hp, annotation_legend_list = lgd_list) 
dev.off()


#==============================================================================#
# Figure S8A ----
#==============================================================================#
cluster_col <- c("B" = "#A6CEE3",
                 "cDC1" = "#1F78B4",   "cDC2" = "#B2DF8A",  
                 "Macrophage" = "#33A02C", "Monocyte" = "#FB9A99",
                 "NK" = "#E31A1C",  "migDC" = "#FDBF6F",  "Mo/Mac" = "#FF7F00",  
                 "pDC" = "#CAB2D6",
                 "T cells" = "#6A3D9A" 
)


pdf("/home/aoill/plots/00_19130_CART/Figure_S7_B.pdf",
    width = 10,
    height = 5)
p1 <- DimPlot(csf_pbmc, group.by = "ct_final_all_simplified", 
              reduction = "umap.integrated.rpca",
              label = F, raster=FALSE,
              split.by = "Sample_Type") +
  scale_color_manual(values = cluster_col) #+ NoLegend()
pa <- LabelClusters(p1, id = "ct_final_all_simplified", box = TRUE, size = 3.5, label.size = 0.5, 
                    box.padding = 0.1, seed = 309,
                    alpha = 0.8
)
pa
dev.off()


#==============================================================================#
# Figure S8B ----
#==============================================================================#
pdf("/home/aoill/plots/00_19130_CART/Figure_S7_B.pdf",
    width = 4,
    height = 3)
dittoBarPlot(
  object = csf_pbmc,
  var = "ct_final_all_simplified",
  group.by = "Sample_Type", 
  color.panel = cluster_col, # color based on cell types
  main = "",
  #x.reorder = c(2,1)
) #+ coord_flip()
dev.off()


#==============================================================================#
# Figure S8C ----
#==============================================================================#
# this is table S16
all_CT_results_one_df <- read_csv("/home/aoill/tables/csf_pbmc_correlations_2026_01_13_exp_fix.csv")

plot_data <- all_CT_results_one_df %>%
  filter(FDR <= 0.05) %>%
  group_by(`Cell Type`) %>%
  dplyr::summarise(GeneCount = n()) %>%
  ungroup()


pdf("/home/aoill/plots/00_19130_CART/Figure_S7_C.pdf",
    width = 4,
    height = 4.5)
ggplot(plot_data, aes(x = reorder(`Cell Type`, -GeneCount), y = GeneCount)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = GeneCount), vjust = -0.5) + # This adds the numbers (500, 100) above bars
  labs(
    x = "Cell Type",
    y = "Significant Gene Count",
    title = "Significant Genes by Cell Type"
  ) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


#==============================================================================#
# Figure S10A ----
#==============================================================================#
pdf("/home/aoill/plots/00_19130_CART/Figure_S10_A.pdf",
    width = 13 , 
    height = 10)
DimPlot(csf_lym, group.by = "CART", reduction = "umap.integrated.rpca",
        split.by = "UPN_Cycle",
        ncol = 7) +
  scale_color_manual(values = c("grey", "red")) + 
  ggtitle("CAR positivity") +
  guides(color=guide_legend(ncol=1, override.aes = list(size=2))) +
  theme(
    legend.position="right",
    #legend.position="bottom",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12)) + coord_fixed()
dev.off()


#==============================================================================#
# Figure S10B ----
#==============================================================================#
upn_colors <- c(
  "514" = "#e834eb",
  "515" = "#9d4f9e",
  "574" = "#c2087e",
  
  "689" = "#17b4bd",
  "692" = "#0a708c",
  "705" = "#0810ff",
  "716" = "#7a95f0"
)

proportion_results <- csf_lym_T@meta.data %>%
  # Group the data by the sample column
  dplyr::group_by(UPN_Cycle) %>% 
  # Summarize the data for each sample:
  dplyr::summarise(
    # Count the total number of cells in the sample
    Total_Cells = n(), 
    # Count the number of positive cells
    Positive_Cells = sum(CART == "Positive"), 
    # Calculate the proportion
    Proportion_Positive = Positive_Cells / Total_Cells 
  ) %>%
  # Ungroup the data (good practice after summarizing)
  ungroup()

#colnames(proportion_results)
proportion_results$UPN_Cycle <- factor(x = proportion_results$UPN_Cycle, 
                                       levels = c("514_2", "514_5", "514_8", "514_11", "514_13",
                                                  "515_4", "515_8",
                                                  "574_1", "574_3", "574_8",  
                                                  "689_1", "689_4", "689_8", "689_12", "689_17", 
                                                  "692_2", "692_4",
                                                  "705_5", "705_8",  
                                                  "716_2", "716_4", "716_8"))

proportion_results <- proportion_results %>%
  separate(col = UPN_Cycle, 
           into = c("UPN", "Cycle"), 
           sep = "_")
proportion_results$UPN_Cycle <- paste(proportion_results$UPN, proportion_results$Cycle, sep = "_")

proportion_results$Cycle <- as.numeric(proportion_results$Cycle)

# add in lymphodepletion info
proportion_results$Lymphodepletion <- NA
proportion_results$Lymphodepletion[proportion_results$UPN == "514"] <- "Non-Lymphodepleted"
proportion_results$Lymphodepletion[proportion_results$UPN == "515"] <- "Non-Lymphodepleted"
proportion_results$Lymphodepletion[proportion_results$UPN == "574"] <- "Non-Lymphodepleted"
proportion_results$Lymphodepletion[proportion_results$UPN == "689"] <- "Lymphodepleted"
proportion_results$Lymphodepletion[proportion_results$UPN == "692"] <- "Lymphodepleted"
proportion_results$Lymphodepletion[proportion_results$UPN == "705"] <- "Lymphodepleted"
proportion_results$Lymphodepletion[proportion_results$UPN == "716"] <- "Lymphodepleted"


pdf("/home/aoill/plots/00_19130_CART/Figure_S10_B.pdf",
    width = 8, height = 3)
ggplot(proportion_results, aes(x = Cycle, y = Proportion_Positive, color = UPN)) +
  geom_point() +
  geom_line() +
  #geom_smooth(method = "lm", se = FALSE) + 
  ylab("Proportion CAR+ T cells") +
  scale_x_continuous(breaks = seq(from = 1, to = max(proportion_results$Cycle), by = 1)) + 
  scale_color_manual(values = upn_colors) +
  facet_wrap(~ Lymphodepletion)
dev.off()


#==============================================================================#
# Figure S13A ----
#==============================================================================#
a <- DimPlot(tumor_imm, group.by = "ct_final", 
             reduction = "umap.integrated.rpca") +
  ggtitle("Immune") + NoLegend() + 
  coord_fixed() +
  scale_color_manual(values = ct_colors) 
pa <- LabelClusters(a, id = "ct_final", box = TRUE, size = 3.5, label.size = 0.5, 
                    box.padding = 1, seed = 309,
                    alpha = 0.8
)

pdf("/home/aoill/plots/00_19130_CART/Figure_S13_A.pdf",
    width = 5, height = 5)
pa
dev.off()


#==============================================================================#
# Figure S13B ----
#==============================================================================#
# CART UMAP
tumor_imm_tcells <- subset(tumor_imm, subset = ct_final == "T (Activated)")
pb <- DimPlot(tumor_imm_tcells, reduction = "umap.integrated.rpca", 
              group.by = "CART",
              raster=FALSE,
              label = F) + coord_fixed() +
  scale_color_manual(values = c("grey", "red")) + theme(legend.position = "bottom")


pdf("/home/aoill/plots/00_19130_CART/Figure_S13_B.pdf",
    width = 5, height = 5)
pb
dev.off()


#==============================================================================#
# Figure S13C ----
#==============================================================================#
pc <- DimPlot(tumor_imm_tcells, 
              group.by = "cloneSize", 
              raster=FALSE) + 
  scale_color_manual(values = c(#"#F0F921FF", "#FA9E3BFF", 
    "#D8576BFF", 
    "#287D8EFF", alpha("#0D0887FF",0.25), 
    alpha("lightgrey", 0.1)), na.value=alpha("lightgrey", 0.1)) + 
  coord_fixed() + 
  theme(legend.position = "bottom") + guides(
    color = guide_legend(ncol = 2,
                         override.aes = list(size = 3))
  ) 

pdf("/home/aoill/plots/00_19130_CART/Figure_S13_C.pdf",
    width = 5, height = 5)
pc
dev.off()


pdf("/home/aoill/plots/00_19130_CART/Figure_S13_ABC.pdf",
    width = 11,
    height = 4)
ggarrange(pa, pb, pc, 
          ncol = 3,
          align = "hv")
dev.off()


#==============================================================================#
# Figure S13D ----
#==============================================================================#
# Subset TCRs for 574 from tumor
tmp_meta <- tumor_imm@meta.data
tmp_meta_sub <- tmp_meta %>%
  dplyr::select(UPN, Sample_Type, CTstrict, clonalFrequency)
rownames(tmp_meta) <- NULL
tmp_meta_sub <- na.omit(tmp_meta_sub)
# There could be duplicated rows in this metadata (clonotypes with freq >1 will 
# come up multiple times)
tmp_meta_sub_unique <- unique(tmp_meta_sub)

tcrs_574 <- tmp_meta_sub_unique %>%
  dplyr::filter(UPN == "574") %>%
  dplyr::pull(CTstrict) 


# Pull out metadata from csf object and only keep the necessary TCR info
csf_obj_TCR_meta <- csf_lym@meta.data
csf_obj_TCR_meta_sub <- csf_obj_TCR_meta %>%
  dplyr::select(UPN, Sample_Type, Cycle, Day, CTstrict, 
                clonalFrequency, Frequency_norm,
                Frequency_norm_cat, Frequency_norm_cat_2)
rownames(csf_obj_TCR_meta_sub) <- NULL
# remove NAs
csf_obj_TCR_meta_subNA <- csf_obj_TCR_meta_sub[is.na(csf_obj_TCR_meta_sub$clonalFrequency) == F, ]

# Get expanded and unexpanded TCRs from CSF
csf_TCRs_exp_574 <- csf_obj_TCR_meta_subNA %>% 
  filter(
    Frequency_norm_cat_2 == "Expanded",
    UPN == 574) %>% 
  pull(CTstrict) %>% unique()

csf_TCRs_574_unexp <- csf_obj_TCR_meta_subNA %>%
  filter(UPN == 574) %>%
  filter(!(CTstrict %in% csf_TCRs_exp_574)) %>%
  pull("CTstrict") %>% unique()

# Get TCRs from PBMCs and product
pbmc_obj_meta <- pbmc_lym@meta.data
product_obj_meta <- product_obj@meta.data

pbmc_obj_meta_NA_rm <- pbmc_obj_meta[is.na(pbmc_obj_meta$clonalFrequency) == F, ]
pbmcs_TCRs_574 <- pbmc_obj_meta_NA_rm %>%
  filter(UPN == 574) %>%
  pull("CTstrict") %>% unique()

product_obj_meta_NA_rm <- product_obj_meta[is.na(product_obj_meta$clonalFrequency) == F, ]
product_TCRs_574 <- product_obj_meta_NA_rm %>%
  filter(UPN == 574) %>%
  pull("CTstrict") %>% unique()

# Make list of TCRs to plot
tcr_list_574_all <- list(Tumor = tcrs_574, Expanded_CSF = csf_TCRs_exp_574,
                         Unexpanded_CSF = csf_TCRs_574_unexp,
                         Product = product_TCRs_574, PBMC = pbmcs_TCRs_574)
m_574_all <- make_comb_mat(tcr_list_574_all)
mat_574_all <-list_to_matrix(tcr_list_574_all)
m_574_all <- make_comb_mat(mat_574_all)

ptest3 <- UpSet(m_574_all, 
                #column_title = expression(bold("Expanded TCRs in CSF")), 
                column_title = expression(bold("UPN 574")), 
                set_order = c("Product", "PBMC", "Unexpanded_CSF", "Expanded_CSF", "Tumor"),
                comb_order = order(comb_degree(m_574_all), -comb_size(m_574_all)),
                left_annotation = upset_left_annotation(m_574_all,
                                                        width = unit(1, "cm")),
                top_annotation = upset_top_annotation(m_574_all, add_numbers = TRUE,
                                                      annotation_name_rot = 90),
                #show_row_names = FALSE
)
ptest3

p3_ann <- grid.grab()


pdf("/home/aoill/plots/00_19130_CART/Figure_S13_D.pdf",
    width = 6, height = 3)
ptest3
dev.off() 


#==============================================================================#
# Figure S13E ----
#==============================================================================#
pe <- dittoBarPlot(
  object = tumor_imm,
  var = "ct_final",
  group.by = "UPN_Tumor_info", 
  color.panel = ct_colors, 
  main = "",
  x.reorder = c(3,2,1),
  #var.labels.reorder = c(7,1,2,3,4,5,6,8)
) + 
  xlab("") +
  #geom_vline(xintercept = 2.5, linetype="dashed", size= 1) + 
  NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

pe

pdf("/home/aoill/plots/00_19130_CART/Figure_S13_E.pdf",
    width = 3, height = 5)
pe
dev.off()


#==============================================================================#
# Figure S14 ----
#==============================================================================#
# These are Seurat object generated from the processing scripts
seurat_list_csf <- readRDS("/scratch/aoill/projects/CAR-T/00_new_2025/rds_files/CSF_seurat_list_all_2025_08_25.rds")

seurat_merge_csf <- merge(x = seurat_list_csf[[1]], y = seurat_list_csf[2:length(seurat_list_csf)])
seurat_merge_csf@meta.data$New_Ident <- "SeuratProject"
Idents(seurat_merge_csf) <- "New_Ident"

seurat_list_pbmc <- readRDS("/scratch/aoill/projects/CAR-T/00_new_2025/rds_files/PBMC_seurat_list_all_2025_08_26.rds")

seurat_merge_pbmc <- merge(x = seurat_list_pbmc[[1]], y = seurat_list_pbmc[2:length(seurat_list_pbmc)])
seurat_merge_pbmc@meta.data$New_Ident <- "SeuratProject"
Idents(seurat_merge_pbmc) <- "New_Ident"

seurat_list_product <- readRDS("/scratch/aoill/projects/CAR-T/00_new_2025/rds_files/product_seurat_list_unfiltered.rds")

seurat_merge_product <- merge(x = seurat_list_product[[1]], y = seurat_list_product[2:length(seurat_list_product)])
seurat_merge_product@meta.data$New_Ident <- "SeuratProject"
Idents(seurat_merge_product) <- "New_Ident"

seurat_list_tumor <- readRDS("/scratch/aoill/projects/CAR-T/00_new_2025/rds_files/tumor_seurat_merge_unfiltered.rds")
seurat_merge_tumor <- seurat_list_tumor
seurat_merge_tumor@meta.data$New_Ident <- "SeuratProject"
Idents(seurat_merge_tumor) <- "New_Ident"


pdf("/home/aoill/plots/00_19130_CART/Figure_S14.pdf", width = 12, height = 6)

par(mfrow=c(2,4))

# CSF %mt log(nFeature_RNA)
smoothScatter(seurat_merge_csf@meta.data$percent.mt, log(seurat_merge_csf@meta.data$nFeature_RNA),
              xlab = "% MT", ylab = "log(nFeature_RNA)",
              main = "CSF")
abline(h = log(650), v = 10)
text(15,log(750), "nFeature_RNA = 650,\npercent.mt = 10", adj = c(0, -.1))

# PBMC %mt log(nFeature_RNA)
smoothScatter(seurat_merge_pbmc@meta.data$percent.mt, log(seurat_merge_pbmc@meta.data$nFeature_RNA),
              xlab = "% MT", ylab = "log(nFeature_RNA)",
              main = "PBMCs")
abline(h = log(650), v = 10)
text(15,log(750), "nFeature_RNA = 650,\npercent.mt = 10", adj = c(0, -.1))

# Product %mt log(nFeature_RNA)
smoothScatter(seurat_merge_product@meta.data$percent.mt, log(seurat_merge_product@meta.data$nFeature_RNA),
              xlab = "% MT", ylab = "log(nFeature_RNA)",
              main = "Product")
abline(h = log(1300), v = 10)
text(15,log(1400), "nFeature_RNA = 1300,\npercent.mt = 10", adj = c(0, -.1))

# Tumor %mt log(nFeature_RNA)
smoothScatter(seurat_merge_tumor@meta.data$percent.mt, log(seurat_merge_tumor@meta.data$nFeature_RNA),
              xlab = "% MT", ylab = "log(nFeature_RNA)",
              main = "Tumor")
abline(h = log(1500), v = 10)
text(15,log(1600), "nFeature_RNA = 1500,\npercent.mt = 10", adj = c(0, -.1))



# CSF %mt log(nCount_RNA)
smoothScatter(seurat_merge_csf@meta.data$percent.mt, log(seurat_merge_csf@meta.data$nCount_RNA),
              xlab = "% MT", ylab = "log(nCount_RNA)",
              main = "CSF")
abline(h = log(1200), v = 10)
text(15,log(1300), "nCount_RNA = 1200,\npercent.mt = 10", adj = c(0, -.1))

# PBMC %mt log(nCount_RNA)
smoothScatter(seurat_merge_pbmc@meta.data$percent.mt, log(seurat_merge_pbmc@meta.data$nCount_RNA),
              xlab = "% MT", ylab = "log(nCount_RNA)",
              main = "PBMCs")
abline(h = log(1200), v = 10)
text(15,log(1300), "nCount_RNA = 1200,\npercent.mt = 10", adj = c(0, -.1))

# Product %mt log(nCount_RNA)
smoothScatter(seurat_merge_product@meta.data$percent.mt, log(seurat_merge_product@meta.data$nCount_RNA),
              xlab = "% MT", ylab = "log(nCount_RNA)",
              main = "Product")
abline(h = log(2500), v = 10)
text(15,log(2600), "nCount_RNA = 2500,\npercent.mt = 10", adj = c(0, -.1))

# Tumor %mt log(nCount_RNA)
smoothScatter(seurat_merge_tumor@meta.data$percent.mt, log(seurat_merge_tumor@meta.data$nCount_RNA),
              xlab = "% MT", ylab = "log(nCount_RNA)",
              main = "Tumor")
abline(h = log(2300), v = 10)
text(15,log(2400), "nCount_RNA = 2300,\npercent.mt = 10", adj = c(0, -.1))

dev.off()



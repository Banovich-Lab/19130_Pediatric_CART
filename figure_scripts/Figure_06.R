#==============================================================================#
# Author : Angela M. Oill, aoill@tgen.org
# Date: 2023/06/29
# Project: Pediatric CAR-T 
# Description: Figure 6
#==============================================================================#


#==============================================================================#
# Load Libraries, helper modules and functions ----
#==============================================================================#
library(Seurat)
library(scRepertoire)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(dittoSeq)
library(scProportionTest)
library(UpSetR)


#==============================================================================#
# Read in object ----
#==============================================================================#
obj <- readRDS("/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/Tumor_immune_cells_obj_2023.rds")
obj_Tcells <- readRDS("/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/Tumor_T_cells_obj_2023.rds")

csf_obj_T <- readRDS("/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/CSF_T_cell_obj_2023.rds")
pbmc_obj_T <- readRDS("/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/PBMC_T_cell_obj_2023.rds")
product_obj_T <- readRDS("/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/product_T_cell_obj_2023.rds")


#==============================================================================#
# Set up colors ----
#==============================================================================#
ct_colors <- c("CD4+ T" = "#e09cb7", "CD8+ T" = "#d747ad", "NK" = "#a35776", 
               "Treg" = "#d54170", 
               "T cells" = "#de7aac",
               "Myeloid" = "#5fa54b",
               "Monocytes" = "#bcd095", "CD14+ Monocytes" = "#6dce52", "FCGR3A+ Monocytes" = "#697f3f",
               "Macrophages" = "#c8d048",
               "DC" = "#de9e35","cDC" = "#de9e35", "cDCs" = "#de9e35", "cDC1" = "#927754", "cDC2" = "#e3bb87", 
               "mDC" = "#a0712b",
               "pDC" = "#55daa2",  
               "B cells" = "#4a8773", "B" = "#4a8773", 
               "Plasma" = "#89d7c8", 
               "Proliferating" = "#c7624f"
               
)

tumor_tmp_CTs <- as.character(levels(as.factor(obj@meta.data$CT)))
tumor_tmp_color <- setNames(names(ct_colors), ct_colors) # swap names and values so I can subset
tumor_tmp_color_sub <- tumor_tmp_color[tumor_tmp_color %in% c(tumor_tmp_CTs)] # subset ct_color to only keep colors in CSF UMAP
tumor_ct_colors_sub <- setNames(names(tumor_tmp_color_sub), tumor_tmp_color_sub) # re-swap names and values so that values are the color names 

# TCR colors
colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
                                            "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                            "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))


#==============================================================================#
# Plot ----
#==============================================================================#

#-----------#
## UMAP ----
#-----------#
a <- DimPlot(obj, group.by = "CT", 
             reduction = "integrated_sct_umap") +
  scale_color_manual(values = tumor_ct_colors_sub) +
  ggtitle("Tumor - Immune Cells") + NoLegend()
pa <- LabelClusters(a, id = "CT", box = TRUE, size = 3.5, label.size = 0.1, 
                    box.padding = 0.5, seed = 309,
                    color = c("black", "white", rep("black", 4)),
                    alpha = 0.8
)


#--------------------------------------------------------#
## Stacked bar plot by patient and pre/post treatment ----
#--------------------------------------------------------#
pe <- dittoBarPlot(
  object = obj,
  var = "CT",
  group.by = "UPN_Tumor_info", 
  color.panel = tumor_ct_colors_sub, 
  main = "",
  x.reorder = c(3,2,1)
) + 
  xlab("") +
  #geom_vline(xintercept = 2.5, linetype="dashed", size= 1) + 
  NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


#-------------------#
## T cell UMAPS ----
#-------------------#
pb <- DimPlot(obj_Tcells, group.by = "T_cell_subset", 
              reduction = "integrated_sct_umap",
              label = F) +
  ggtitle("T cells") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        legend.text=element_text(size=10))


pc <- DimPlot(obj_Tcells, group.by = "CART") + 
  scale_color_manual(values=c("grey", "red")) +
  ggtitle("CAR Positivity") +
  theme(legend.position="right",
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        legend.text=element_text(size=10))


pd <- DimPlot(obj_Tcells, group.by = "cloneType", reduction = "integrated_sct_umap") +
  scale_color_manual(values = c("Medium (5 < X <= 20)" = "#FA9E3BFF", 
                                "Small (1 < X <= 5)" = "#D8576BFF", 
                                "Single (0 < X <= 1)" = "#0D0887FF", 
                                "NA" = "darkgrey")) + 
  ggtitle("TCR Frequency") +
  theme(
    legend.position="right",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.text=element_text(size=10)) +
  guides(color=guide_legend(ncol=1, override.aes = list(size=2)))


#-----------------------------#
## Upset plot for UPN 574 ----
#-----------------------------#
# Subset TCRs for 574 from tumor
tmp_meta <- obj_Tcells@meta.data
tmp_meta_sub <- tmp_meta %>%
  dplyr::select(UPN, Sample_Type, CTstrict, Frequency)
rownames(tmp_meta) <- NULL
tmp_meta_sub <- na.omit(tmp_meta_sub)
# There could be duplicated rows in this metadata (clonotypes with freq >1 will 
# come up multiple times)
tmp_meta_sub_unique <- unique(tmp_meta_sub)

tcrs_574 <- tmp_meta_sub_unique %>%
  dplyr::filter(UPN == "574") %>%
  dplyr::pull(CTstrict) 


# Pull out metadata from csf object and only keep the necessary TCR info
csf_obj_T_TCR_meta <- csf_obj_T@meta.data
csf_obj_T_TCR_meta_sub <- csf_obj_T_TCR_meta %>%
  dplyr::select(UPN, Sample_Type, Cycle, Cycle_Day, CTstrict, 
                Frequency, Frequency_norm,
                Frequency_norm_cat, Frequency_norm_cat_2)
rownames(csf_obj_T_TCR_meta_sub) <- NULL
# remove NAs
csf_obj_T_TCR_meta_subNA <- csf_obj_T_TCR_meta_sub[is.na(csf_obj_T_TCR_meta_sub$Frequency) == F, ]

# Get expanded and unexpanded TCRs from CSF
csf_TCRs_exp_574 <- csf_obj_T_TCR_meta_subNA %>% 
  filter(
    Frequency_norm_cat_2 == "Expanded",
    UPN == 574) %>% 
  pull(CTstrict) %>% unique()

csf_TCRs_574_unexp <- csf_obj_T_TCR_meta_subNA %>%
  filter(UPN == 574) %>%
  filter(!(CTstrict %in% csf_TCRs_exp_574)) %>%
  pull("CTstrict") %>% unique()

# Get TCRs from PBMCs and product
pbmc_obj_T_meta <- pbmc_obj_T@meta.data
product_obj_T_meta <- product_obj_T@meta.data

pbmc_obj_T_meta_NA_rm <- pbmc_obj_T_meta[is.na(pbmc_obj_T_meta$Frequency) == F, ]
pbmcs_TCRs_574 <- pbmc_obj_T_meta_NA_rm %>%
  filter(UPN == 574) %>%
  pull("CTstrict") %>% unique()

product_obj_T_meta_NA_rm <- product_obj_T_meta[is.na(product_obj_T_meta$Frequency) == F, ]
product_TCRs_574 <- product_obj_T_meta_NA_rm %>%
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

plot_upsets <- ggarrange(p3_ann,
                         ncol = 1,
                         nrow = 1,
                         labels = c("F")
)


#==============================================================================#
# Arrange plots and output pdf ----
#==============================================================================#
pdummy <- ggplot() + theme_void()

pdf("/home/aoill/plots/Figure_06", width = 7.5,
    height = 11)
ggarrange(ggarrange(ggarrange(pa, pdummy, ncol = 1, nrow = 2, heights = c(1,.14)), 
                    pe, labels = c("A", "E"),
                    widths = c(1,.7)), 
          ggarrange(pb, pc, pd, ncol = 3, nrow = 1, 
                    align = "hv", legend = "bottom",
                    labels = c("B", "C", "D")),
          plot_upsets,
          ncol = 1, nrow = 3,
          heights = c(1, .75, .7))
dev.off()


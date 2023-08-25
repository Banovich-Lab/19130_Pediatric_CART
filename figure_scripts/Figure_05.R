#==============================================================================#
# Author : Angela M. Oill, aoill@tgen.org
# Date: 2023/06/27
# Project: Pediatric CAR-T 
# Description: Figure 5
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
library(UpSetR)


#==============================================================================#
# Read in integrated objects ----
#==============================================================================#
csf_obj_T <- readRDS("/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/CSF_T_cell_obj_2023.rds")
pbmc_obj_T <- readRDS("/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/PBMC_T_cell_obj_2023.rds")
product_obj_T <- readRDS("/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/product_T_cell_obj_2023.rds")


#==============================================================================#
# Make plots ----
#==============================================================================#
# This will be TCR data for CSF
csf_obj_T@meta.data$UPN_Cycle <- factor(x = csf_obj_T@meta.data$UPN_Cycle, 
                                        levels = c("514_5", "514_8", "514_11",
                                                   "515_4", "515_8", 
                                                   "574_3", "574_8"))

#----------------------#
## Split UMAP CART ----
#----------------------#
pa <- DimPlot(csf_obj_T, group.by = "CART", reduction = "integrated_sct_umap",
              split.by = "UPN_Cycle",
              ncol = 7) +
  scale_color_manual(values = c("grey", "red")) + 
  ggtitle("CAR positivity") +
  guides(color=guide_legend(ncol=1, override.aes = list(size=2))) +
  theme(
    legend.position="right",
    #legend.position="bottom",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12))


#-------------------------#
## Split UMAP TCR freq ----
#-------------------------#
pb <- DimPlot(csf_obj_T, group.by = "Frequency_norm_cat", reduction = "integrated_sct_umap",
              split.by = "UPN_Cycle",
              ncol = 7) +
  scale_color_manual(values = c("#F0F921FF", "#FA9E3BFF", "#D8576BFF", alpha("#0D0887FF",0.25)), na.value=alpha("grey", 0.25)) + 
  #ggtitle("CSF") +
  ggtitle("Normalized TCR Frequency") +
  guides(color=guide_legend(ncol=1, override.aes = list(size=3))) +
  theme(
    legend.position="right",
    #legend.position="bottom",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12)) 



#---------------------#
## Proportion test ----
#---------------------#
# Non-exp vs exp
prop_test <- sc_utils(csf_obj_T)
prop_test <- permutation_test(
  prop_test, 
  cluster_identity = "ct_T_states",
  sample_1 = "Not Expanded", sample_2 = "Expanded",
  sample_identity = "Frequency_norm_cat_2"
)

prop_test@results$permutation

pc <- permutation_plot(prop_test) + 
  ggtitle("CSF\nNot Expanded vs Expanded") +
  ylab("") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5, face = "bold",
                                  size = 10)) 


#------------------#
## Upset plots ----
#------------------#
# Pull out metadata from csf object and only keep the necessary TCR info
csf_obj_T_TCR_meta <- csf_obj_T@meta.data
csf_obj_T_TCR_meta_sub <- csf_obj_T_TCR_meta %>%
  dplyr::select(UPN, Sample_Type, Cycle, Cycle_Day, CTstrict, 
                Frequency, Frequency_norm,
                Frequency_norm_cat, Frequency_norm_cat_2)
rownames(csf_obj_T_TCR_meta_sub) <- NULL
# remove NAs
csf_obj_T_TCR_meta_subNA <- csf_obj_T_TCR_meta_sub[is.na(csf_obj_T_TCR_meta_sub$Frequency) == F, ]


# Pull out metadata from PBMCs and Product (to get TCR clonotype info later)
pbmc_obj_T_meta <- pbmc_obj_T@meta.data
product_obj_T_meta <- product_obj_T@meta.data

# Pull out TCR name only (CTstrict)
pbmc_obj_T_meta_NA_rm <- pbmc_obj_T_meta[is.na(pbmc_obj_T_meta$Frequency) == F, ]
pbmcs_TCRs_514 <- pbmc_obj_T_meta_NA_rm %>%
  filter(UPN == 514) %>%
  pull("CTstrict") %>% unique()

pbmcs_TCRs_515 <- pbmc_obj_T_meta_NA_rm %>%
  filter(UPN == 515) %>%
  pull("CTstrict") %>% unique()

pbmcs_TCRs_574 <- pbmc_obj_T_meta_NA_rm %>%
  filter(UPN == 574) %>%
  pull("CTstrict") %>% unique()


product_obj_T_meta_NA_rm <- product_obj_T_meta[is.na(product_obj_T_meta$Frequency) == F, ]
product_TCRs_514 <- product_obj_T_meta_NA_rm %>%
  filter(UPN == 514) %>%
  pull("CTstrict") %>% unique()

product_TCRs_515 <- product_obj_T_meta_NA_rm %>%
  filter(UPN == 515) %>%
  pull("CTstrict") %>% unique()

product_TCRs_574 <- product_obj_T_meta_NA_rm %>%
  filter(UPN == 574) %>%
  pull("CTstrict") %>% unique()


# Get expanded TCRs from CSF
csf_TCRs_exp_514 <- csf_obj_T_TCR_meta_subNA %>% 
  filter(
    Frequency_norm_cat_2 == "Expanded",
    UPN == 514) %>% 
  pull(CTstrict) %>% unique()

csf_TCRs_exp_515 <- csf_obj_T_TCR_meta_subNA %>% 
  filter(
    Frequency_norm_cat_2 == "Expanded",
    UPN == 515) %>% 
  pull(CTstrict) %>% unique()

csf_TCRs_exp_574 <- csf_obj_T_TCR_meta_subNA %>% 
  filter(
    Frequency_norm_cat_2 == "Expanded",
    UPN == 574) %>% 
  pull(CTstrict) %>% unique()


# Make a list of the combinations to be plotted in upset plot
tcr_list_514 <- list(Expanded_CSF = csf_TCRs_exp_514, Product = product_TCRs_514, PBMC = pbmcs_TCRs_514)
tcr_list_515 <- list(Expanded_CSF = csf_TCRs_exp_515, Product = product_TCRs_515, PBMC = pbmcs_TCRs_515)
tcr_list_574 <- list(Expanded_CSF = csf_TCRs_exp_574, Product = product_TCRs_574, PBMC = pbmcs_TCRs_574)


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


# Make plots
ptest1 <- UpSet(m_514, column_title = expression(bold("UPN 514")), 
                comb_order = order(comb_degree(m_514), -comb_size(m_514)),
                left_annotation = upset_left_annotation(m_514,
                                                        width = unit(1, "cm")),
                top_annotation = upset_top_annotation(m_514, add_numbers = TRUE,
                                                      annotation_name_rot = 90),
                show_row_names = FALSE
)
ptest1
p1_ann <- grid.grab()

ptest2 <- UpSet(m_515, column_title = expression(bold("UPN 515")), 
                comb_order = order(comb_degree(m_514), -comb_size(m_514)),
                left_annotation = upset_left_annotation(m_515,
                                                        width = unit(1, "cm")),
                top_annotation = upset_top_annotation(m_515, add_numbers = TRUE,
                                                      annotation_name_rot = 90),
                show_row_names = FALSE
)
ptest2
p2_ann <- grid.grab()

ptest3 <- UpSet(m_574, column_title = expression(bold("UPN 574")), 
                #comb_order = order(comb_size(m_574)),
                comb_order = order(comb_degree(m_574), -comb_size(m_574)),
                left_annotation = upset_left_annotation(m_574,
                                                        width = unit(1, "cm")),
                top_annotation = upset_top_annotation(m_574, add_numbers = TRUE,
                                                      annotation_name_rot = 90)
                #show_row_names = FALSE
)
ptest3
p3_ann <- grid.grab()


#------------------------------------#
## Arrange plots and output pdf ----
#------------------------------------#
pdf("/home/aoill/plots/Figure_05.pdf", width = 11,
    height = 8)
ggarrange(ggarrange(pa, pb,
                    ncol = 1,
                    nrow = 2,
                    align = "hv",
                    labels = c("A", "B")),
          ggarrange(pc, ggarrange(
            p1_ann, p2_ann, p3_ann, 
            nrow = 1,
            ncol = 3,
            align = "hv",
            labels = c("D"),
            #widths = c(.8,.8,1)), 
            widths = c(.65,.65,1)), 
            ncol = 2, widths = c(.5, 1), #align = "v",
            labels = c("C")),
          ncol = 1,
          nrow = 2,
          heights = c(1, .6)
)
dev.off()


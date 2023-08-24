#==============================================================================#
# Author : Angela M. Oill, aoill@tgen.org
# Date: 2023/06/09
# Project: Pediatric CAR-T 
# Description: Integration of tumor, CSF, and PBMC samples
#==============================================================================#


#==============================================================================#
# Load Libraries, helper modules and functions ----
#==============================================================================#
#library(Matrix)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(dittoSeq)
library(scProportionTest)
library(scales)

source("/home/aoill/projects/SingleCellBestPractices/scripts/preprocessing_qc_module.R")
source("/home/aoill/projects/SingleCellBestPractices/scripts/helper_functions_module.R")


#==============================================================================#
# Read in objects ----
#==============================================================================#
# This data has already been integrated but want to merge all data together 
# and integrate
# Tumor immune cells only
tumor_obj <- readRDS("/scratch/aoill/projects/CAR-T/00_final/rds_files/Tumor_integrated_immune_cells_20230515.rds")
DimPlot(tumor_obj, group.by = "CT", label = T)

csf_obj <- readRDS("/scratch/aoill/projects/CAR-T/00_final/rds_files/CSF_integrated_CTs_CART_TCR_20230445.rds")
DimPlot(csf_obj, group.by = "CT", label = T)

pbmc_obj <- readRDS("/scratch/aoill/projects/CAR-T/00_final/rds_files/pbmc_integrated_final_20230505_CT_TCR_other_meta.rds")
DimPlot(pbmcs_obj, group.by = "CT", label = T)


#==============================================================================#
# Prep objects into a list for integration ----
#==============================================================================#
#-------------------------------------#
## Subset 574 and 625 from objects ----
#-------------------------------------#
# Tumor only has these samples
unique(tumor_obj@meta.data$UPN)

# Pull out 574 from CSF object
unique(csf_obj@meta.data$UPN)
csf_obj_sample_sub <- subset(csf_obj, subset = UPN == "574")


# Pull out 574 and 625 from PBMC object
unique(pbmc_obj@meta.data$UPN)
pbmc_obj_sample_sub <- subset(pbmc_obj, subset = UPN %in% c("574", "625"))


#---------------------------------------#
## Split each object by Cell_Prefix ----
#---------------------------------------#
tumor_obj_batches_list <- SplitObject(tumor_obj, split.by = "Cell_Prefix")
csf_obj_sample_sub_batches_list <- SplitObject(csf_obj_sample_sub, split.by = "Cell_Prefix")
pbmc_obj_sample_sub_batches_list <- SplitObject(pbmc_obj_sample_sub, split.by = "Cell_Prefix")


#---------------------------------------#
## Combine each list ----
#---------------------------------------#
batches_all <- c(tumor_obj_batches_list, csf_obj_sample_sub_batches_list, pbmc_obj_sample_sub_batches_list)


#==============================================================================#
# Integration ----
#==============================================================================#
# Re-run SCTransform
batches_all <- lapply(batches_all, function(xx){
  # Rerunning SCTransform
  DefaultAssay(xx) <- "SoupX_RNA"
  xx <- SCTransform(xx, 
                    method = "glmGamPoi",
                    vars.to.regress = c("S.Score", "G2M.Score"), 
                    vst.flavor = "v2",
                    verbose = T)
  
})


# Integrate
all_integrated_obj <- basic_sct_rpca_integration(batches_all, npcs=50, k_weight=100, numfeatures = 2500)

## SAVE ##
saveRDS(all_integrated_obj, "/scratch/aoill/projects/CAR-T/00_final/rds_files/Tumors_CSF_PBMC_574_625_integrated.rds")
#all_integrated_obj <- readRDS("/scratch/aoill/projects/CAR-T/00_final/rds_files/Tumors_CSF_PBMC_574_625_integrated.rds")


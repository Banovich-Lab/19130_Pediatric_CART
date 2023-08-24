#==============================================================================#
# Author : Angela M. Oill, aoill@tgen.org
# Date: 2023/06/09
# Project: Pediatric CAR-T 
# Description: Integration of all CSF and PBMC samples
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
# Set variables ----
#==============================================================================#
set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)


#==============================================================================#
# Read in objects ----
#==============================================================================#
csf_obj <- readRDS("/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/CSF_all_cells_obj_2023.rds")

pbmc_obj <- readRDS("/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/PBMC_all_cells_obj_2023.rds")

#==============================================================================#
# Prep objects into a list for integration ----
#==============================================================================#
#---------------------------------------#
## Split each object by Cell_Prefix ----
#---------------------------------------#
csf_obj_batches_list <- SplitObject(csf_obj, split.by = "Cell_Prefix")
pbmc_obj_batches_list <- SplitObject(pbmc_obj, split.by = "Cell_Prefix")


#---------------------------------------#
## Combine each list ----
#---------------------------------------#
batches_all <- c(csf_obj_batches_list, pbmc_obj_batches_list)


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
saveRDS(all_integrated_obj, "/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/CSF_PBMC_all_cells_obj_2023.rds")


#==============================================================================#
# Author: Angela M. Oill, aoill@tgen.org
# Date: 2023/06/09
# Project: Pediatric CAR-T 
# Description: Figure S10
#==============================================================================#


#==============================================================================#
# Load Libraries, helper modules and functions ----
#==============================================================================#
source("/home/aoill/projects/SingleCellBestPractices/scripts/preprocessing_qc_module.R")


#==============================================================================#
# Set variables ----
#==============================================================================#
set.seed(1234)

# Path to google sheet with metadata for samples
metapath <- "https://docs.google.com/spreadsheets/d/1gRL53qgRBApRHxznwTK_17I1cRlAL0GGf8nMikOJle0/edit#gid=0"
sheet_name <- "Sheet1"


#==============================================================================#
# Prep metadata ----
#==============================================================================#

#--------------------#
## Read metadata ----
#--------------------#
sample_metadata <- load_metadata(metadata_path = metapath, data_type = "googlesheet", sheet_name = sheet_name)

#---------------------------#
## Separate CSF samples ----
#---------------------------#
sample_metadata_not_multiplexed_csf <- sample_metadata %>% 
  filter(Multiplexed == "No") %>%
  filter(Sample_Type == "CSF")

#----------------------------------------------------------------#
## Separate multiplexed samples (PBMCs and product samples) ----
#----------------------------------------------------------------#
sample_metadata_multiplexed <- sample_metadata %>% filter(Multiplexed == "Yes")

sample_metadata_not_multiplexed <- sample_metadata %>% filter(Multiplexed == "No")

#-----------------------------#
## Separate tumor samples ----
#-----------------------------#
sample_metadata_not_multiplexed_tumor <- sample_metadata %>% 
  filter(Multiplexed == "No") %>%
  filter(Sample_Type == "Tumor")


#==============================================================================#
# Read in Seurat data as a list ----
#==============================================================================#

#-----------#
## CSF ----
#-----------#
seurat_list_csf <- prep_seurat_list(
  sample_metadata_not_multiplexed_csf, 
  batch_ID = "Cell_Prefix", 
  cellRanger_path = "CellRanger_path", 
  cell_ID_prefix = "Cell_Prefix", 
  run_soupX = FALSE) # Run SoupX after filtering

#-----------------------------------------------#
## Multiplexed samples (PBMCs and Products) ----
#-----------------------------------------------#
# This function also runs HTODemux for demultiplexing
seurat_list_demultiplexed <- prep_seurat_list_multiplexed(
  sample_metadata_multiplexed, 
  batch_ID = "Cell_Prefix", 
  cellRanger_path = "CellRanger_path", 
  cell_ID_prefix = "Cell_Prefix", 
  CellHashing_Ab = "CellHashing_Ab")

# Keep only singlets
seurat_list_demultiplexed_singlets <- filter_doublets_cellhashing(seurat_list_demultiplexed)

# Since we are performing sample specific filtering, separate out PBMCs and 
# product prior to QC visualization

#------------#
## PBMCs ----
#------------#
# Keep only PBMC samples
pbmc_seurat_list_demultiplexed_singlets <- list()
for (i in 1:length(seurat_list_demultiplexed_singlets)) {
  print(i)
  sub_i <- subset(seurat_list_demultiplexed_singlets[[i]], subset = Sample_Type == "PBMC")
  pbmc_seurat_list_demultiplexed_singlets <- c(pbmc_seurat_list_demultiplexed_singlets, sub_i)
}
names(pbmc_seurat_list_demultiplexed_singlets) <- c("Batch37_filtered", "Batch38_filtered", "Batch39_filtered", "Batch40_filtered")


seurat_list_pbmc <- pbmc_seurat_list_demultiplexed_singlets

#--------------#
## Product ----
#--------------#
# Keep only product samples
product_seurat_list_demultiplexed_singlets <- list()
for (i in 1:length(seurat_list_demultiplexed_singlets)) {
  print(i)
  sub_i <- subset(seurat_list_demultiplexed_singlets[[i]], subset = Sample_Type == "Product")
  product_seurat_list_demultiplexed_singlets <- c(product_seurat_list_demultiplexed_singlets, sub_i)
}
names(product_seurat_list_demultiplexed_singlets) <- c("Batch37_filtered", "Batch38_filtered", "Batch39_filtered")

# There was a product sample that wasn't multiplexed. Read in this sample and 
# add it to the product sample list
sample_metadata_not_multiplexed_626 <- sample_metadata_not_multiplexed %>% filter(FID_GEXFB == "F05565")
seurat_list_626 <- prep_seurat_list(sample_metadata_not_multiplexed_626, batch_ID = "FID_GEXFB", cellRanger_path = "CellRanger_path", cell_ID_prefix = "Cell_Prefix", run_soupX = FALSE)

seurat_list_product <- c(product_seurat_list_demultiplexed_singlets, seurat_list_626)


#--------------#
## Tumors ----
#--------------#
seurat_list_tumor <- prep_seurat_list(
  sample_metadata_not_multiplexed_tumor, 
  batch_ID = "Cell_Prefix", 
  cellRanger_path = "CellRanger_path", 
  cell_ID_prefix = "Cell_Prefix", 
  run_soupX = FALSE) # Run SoupX after filtering


#==============================================================================#
# Merge Suerat list ----
#==============================================================================#

#-----------#
## CSF ----
#-----------#
seurat_merge_csf <- merge(x = seurat_list_csf[[1]], y = seurat_list_csf[2:length(seurat_list_csf)])
seurat_merge_csf@meta.data$New_Ident <- "SeuratProject"
Idents(seurat_merge_csf) <- "New_Ident"

#------------#
## PBMCs ----
#------------#
seurat_merge_pbmc <- merge(x = seurat_list_pbmc[[1]], y = seurat_list_pbmc[2:length(seurat_list_pbmc)])
seurat_merge_pbmc@meta.data$New_Ident <- "SeuratProject"
Idents(seurat_merge_pbmc) <- "New_Ident"

#--------------#
## Product ----
#--------------#
seurat_merge_product <- merge(x = seurat_list_product[[1]], y = seurat_list_product[2:length(seurat_list_product)])
seurat_merge_product@meta.data$New_Ident <- "SeuratProject"
Idents(seurat_merge_product) <- "New_Ident"

#--------------#
## Tumors ----
#--------------#
seurat_merge_tumor <- merge(x = seurat_list_tumor[[1]], y = seurat_list_tumor[2:length(seurat_list_tumor)])
seurat_merge_tumor@meta.data$New_Ident <- "SeuratProject"
Idents(seurat_merge_tumor) <- "New_Ident"


#==============================================================================#
# View QC metrics ----
#==============================================================================#
pdf("/home/aoill/plots/Figure_S10.pdf", width = 12, height = 6)

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

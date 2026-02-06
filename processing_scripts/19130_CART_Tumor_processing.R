#==============================================================================#
# Authors: 
#     Angela M. Oill, aoill@tgen.org
#     Rishika Mudunuri, rmudunuri@tgen.org 
# Date: 2025
# Project: Pediatric CAR-T 
# Description: Preprocessing and integration of tumor samples
#==============================================================================#


#==============================================================================#
# Load Libraries, helper modules and functions ----
#==============================================================================#
library(Seurat)
library(scCustomize)
library(Matrix)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(DropletUtils)
library(scRepertoire)
library(RColorBrewer)
library(scater)
library(randomcoloR)
library(presto)

source("/scratch/rmudunuri/cart_project/preprocessing_qc_module.R")
source("/scratch/rmudunuri/cart_project/helper_functions_module.R")

#==============================================================================#
# Set variables ----
#==============================================================================#
set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)

#==============================================================================#
# Prep metadata ----
#==============================================================================#

#--------------------#
## Read metadata ----
#--------------------#
sample_metadata <- read.csv("/scratch/rmudunuri/cart_project/19130_all_batches_metadata_new.csv")

#-----------------------------------#
## Separate multiplexed samples ----
#-----------------------------------#
sample_metadata_not_multiplexed_tumor <- sample_metadata %>% 
  filter(Multiplexed == "No") %>%
  filter(Sample_Type == "Tumor")


#==============================================================================#
# Read in Seurat data as a list ----
#==============================================================================#
seurat_list <- prep_seurat_list(
  sample_metadata_not_multiplexed_tumor, 
  batch_ID = "Cell_Prefix", 
  cellRanger_path = "CellRanger_path", 
  cell_ID_prefix = "Cell_Prefix", 
  run_soupX = FALSE) # Run SoupX after filtering


#==============================================================================#
# Filter ----
#==============================================================================#
seurat_merge <- merge(x = seurat_list[[1]], y = seurat_list[2:length(seurat_list)])

#saveRDS(seurat_merge, "/scratch/rmudunuri/cart_project/Tumor/Second_version/seurat_merge_unfiltered.rds")

# Tumor %mt log(nFeature_RNA)
smoothScatter(seurat_merge@meta.data$percent.mt, log(seurat_merge@meta.data$nFeature_RNA),
              xlab = "% MT", ylab = "log(nFeature_RNA)",
              main = "CSF")
abline(h = log(1500), v = 10)
text(15,log(1500), "nFeature_RNA = 1500,\npercent.mt = 10", adj = c(0, -.1))


# Tumor %mt log(nCount_RNA)
smoothScatter(seurat_merge@meta.data$percent.mt, log(seurat_merge@meta.data$nCount_RNA),
              xlab = "% MT", ylab = "log(nCount_RNA)",
              main = "CSF")
abline(h = log(2300), v = 10)
text(15,log(2500), "nCount_RNA = 2300,\npercent.mt = 10", adj = c(0, -.1))

seurat_list_filtered <- filter_manual(seurat_list, 
                                      pt_mt = 10, 
                                      nFeature = 1500,
                                      nCount = 2300) 

seurat_list_separated <- list()
for (i in 1:length(seurat_list_filtered)) {
  print(i)
  
  pedi_batch <- seurat_list_filtered[[i]]
  pedi_batch@meta.data$Sample_ID <- paste(pedi_batch@meta.data$UPN, 
                                          pedi_batch@meta.data$Sample_Type_More,
                                          #pedi_batch@meta.data$Cycle, 
                                          #pedi_batch@meta.data$Day, 
                                          pedi_batch@meta.data$Batch,
                                          sep = "_")
  pedi_batch <- Seurat::SplitObject(pedi_batch, split.by = "Sample_ID")
  seurat_list_separated <- c(seurat_list_separated, pedi_batch)
}


for (i in 1:length(seurat_list_separated)) {
  print(paste0("Converting sample ", names(seurat_list_separated[i])))
  obj.sub <- seurat_list_separated[[i]]
  DropletUtils::write10xCounts(path = paste0("/scratch/rmudunuri/cart_project/Tumor/Second_version/soupX/demultiplexed_",
                                             names(seurat_list_separated[i])), 
                               x = obj.sub[["RNA"]]@layers$counts, 
                               barcodes = colnames(obj.sub[["RNA"]]), # cell names
                               gene.id = rownames(obj.sub[["RNA"]]), # Gene identifiers, one per row of X
                               type = "sparse", version="3", overwrite = T)
}


#-----------------------------------------------------#
## Get number of cells before and after filtering ----
#-----------------------------------------------------#
seurat_merge <- merge(x = seurat_list[[1]], y = seurat_list[2:length(seurat_list)])
seurat_merge@meta.data$UPN_Cycle_Day_Sample_Type_Batch <- paste0(seurat_merge@meta.data$UPN, "_",
                                                                 seurat_merge@meta.data$Sample_Type_More, "_",
                                                                 seurat_merge@meta.data$Batch)
table(seurat_merge@meta.data$UPN_Cycle_Day_Sample_Type_Batch)
cell_numbers_unfiltered <- as.data.frame(table(seurat_merge@meta.data$UPN_Cycle_Day_Sample_Type_Batch))
colnames(cell_numbers_unfiltered) <- c("Sample", "Unfiltered")

seurat_merge <- merge(x = seurat_list_filtered[[1]], y = seurat_list_filtered[2:length(seurat_list_filtered)])
seurat_merge@meta.data$UPN_Cycle_Day_Sample_Type_Batch <- paste0(seurat_merge@meta.data$UPN, "_",
                                                                 seurat_merge@meta.data$Sample_Type_More, "_",
                                                                 seurat_merge@meta.data$Batch)
table(seurat_merge@meta.data$UPN_Cycle_Day_Sample_Type_Batch)

cell_numbers_filtered <- as.data.frame(table(seurat_merge@meta.data$UPN_Cycle_Day_Sample_Type_Batch))
colnames(cell_numbers_filtered) <- c("Sample", "Filtered")

cell_numbers_all <- join(cell_numbers_unfiltered, cell_numbers_filtered)
cell_numbers_all$prop_keep <- cell_numbers_all$Filtered/cell_numbers_all$Unfiltered


#cell_numbers_all <- cell_numbers_all %>% separate_wider_delim(Sample, delim = "_", names = c("UPN", "Sample_Type_More", "Batch"))
sum(cell_numbers_all$Unfiltered)
# 47120
sum(cell_numbers_all$Filtered)
# 16687

#print(cell_numbers_all, n = nrow(cell_numbers_all))

#==============================================================================#
# Run SoupX ----
#==============================================================================#
# Extract product samples from these above objects
sample_list <- names(seurat_list_separated)

seurat_list_separated_SoupX <- sapply(sample_list,  function(i){
  print(i)
  # Read in count and droplet data
  d10x_toc <- Read10X(paste0("/scratch/rmudunuri/cart_project/Tumor/Second_version/soupX/demultiplexed_", i))
  
  # Subset metadata based on upn, cycle, day, batch and then grab cellRanger_path
  upn_id <- str_split(i, "_")[[1]][1]
  sample_type_id_1 <- str_split(i, "_")[[1]][2]
  sample_type_id_2 <- str_split(i, "_")[[1]][3]
  sample_type_id_3 <- str_split(i, "_")[[1]][4]
  sample_type_id_4 <- str_split(i, "_")[[1]][5]
  batch_id <- str_split(i, "_")[[1]][6]
  sample_type_more <- paste(sample_type_id_1, sample_type_id_2, sample_type_id_3, sample_type_id_4, sep = "_")
  sample_metadata_not_multiplexed_tumor_i_CR_path <- sample_metadata_not_multiplexed_tumor %>%
    filter(UPN == upn_id) %>% filter(Sample_Type_More == sample_type_more) %>% 
    filter(Batch == batch_id) %>% 
    pull(CellRanger_path)
  
  d10x_tod <- Read10X(paste0(sample_metadata_not_multiplexed_tumor_i_CR_path, "/outs/raw_feature_bc_matrix/"))
  
  # Run SoupX
  sc <- SoupChannel(d10x_tod, d10x_toc, calcSoupProfile = FALSE) 
  sc <- estimateSoup(sc)
  toc_seu <- CreateSeuratObject(d10x_toc)
  toc_seu <- SCTransform(toc_seu, vst.flavor = "v2")
  toc_seu <- RunPCA(toc_seu)
  toc_seu <- RunUMAP(toc_seu, dims = 1:15) 
  toc_seu <- FindNeighbors(toc_seu, dims = 1:15) 
  toc_seu <- FindClusters(toc_seu, resolution = 0.5) 
  ## Add meta data to soupX object
  sc <- setClusters(sc, setNames(toc_seu$seurat_clusters, rownames(toc_seu@meta.data)))
  ## Estimate contamination (automated method)
  message(paste0("Getting autoEstCont for: ", i))
  print(paste0("Getting autoEstCont for: ", i))
  sc <- autoEstCont(sc, forceAccept=TRUE) 
  out <- adjustCounts(sc)
  
  # Create Seurat object using corrected data
  d10x_seu <- CreateSeuratObject(out, assay = "SoupX_RNA")
  d10x_seu[["RNA"]] <- toc_seu@assays[["RNA"]]
  d10x_seu <- PercentageFeatureSet(d10x_seu, pattern = "^MT-", col.name = "percent.mt_RNA", assay = "RNA")
  d10x_seu <- PercentageFeatureSet(d10x_seu, pattern = "^RP[SL]|^MRP[SL]", col.name = "percent.ribo_RNA", assay = "RNA")
  d10x_seu <- PercentageFeatureSet(d10x_seu, pattern = "^MT-", col.name = "percent.mt_SoupX_RNA", assay = "SoupX_RNA")
  d10x_seu <- PercentageFeatureSet(d10x_seu, pattern = "^RP[SL]|^MRP[SL]", col.name = "percent.ribo_SoupX_RNA", assay = "SoupX_RNA")
  
})

names(seurat_list_separated_SoupX) <- sample_list

#seurat_list_separated_SoupX <- seurat_list_separated_SoupX
# Add metadata
batch_ID = "FID_GEXFB"
cellRanger_path = "CellRanger_path"
for (i in names(seurat_list_separated_SoupX)) {
  upn_id <- str_split(i, "_")[[1]][1] 
  
  sample_type_id_1 <- str_split(i, "_")[[1]][2]
  sample_type_id_2 <- str_split(i, "_")[[1]][3]
  sample_type_id_3 <- str_split(i, "_")[[1]][4]
  sample_type_id_4 <- str_split(i, "_")[[1]][5]
  sample_type_more <- paste(sample_type_id_1, sample_type_id_2, sample_type_id_3, sample_type_id_4, sep = "_")
  
  batch_type_id <- str_split(i, "_")[[1]][6]
  
  meta.data <- sample_metadata %>% filter(UPN == upn_id) %>% 
    filter(Sample_Type_More == sample_type_more) %>% 
    filter(Batch == batch_type_id)
  
  # Add sample metadata
  for (m in colnames(meta.data)){
    if (m == cellRanger_path) {
      next
    } else {
      seurat_list_separated_SoupX[[i]]@meta.data[[m]] = meta.data[[m]]
    }
  }
  seurat_list_separated_SoupX[[i]]@meta.data$Sample_ID <- paste(seurat_list_separated_SoupX[[i]]@meta.data$UPN, 
                                                                seurat_list_separated_SoupX[[i]]@meta.data$Sample_Type_More,
                                                                #seurat_list_separated_SoupX[[i]]@meta.data$Cycle, 
                                                                #seurat_list_separated_SoupX[[i]]@meta.data$Day, 
                                                                seurat_list_separated_SoupX[[i]]@meta.data$Batch,
                                                                sep = "_")
  
}

#==============================================================================#
# Run and add cell cycle score ----
#==============================================================================#
add_cell_cycle_score_2 <- function(sample_seurat_list, rna_assay = "RNA"){
  # Add cell cycle score
  s.genes <- cc.genes.updated.2019$s.genes
  g2m.genes <- cc.genes.updated.2019$g2m.genes
  
  for (i in 1:length(sample_seurat_list)){
    message(names(sample_seurat_list)[i])
    DefaultAssay(sample_seurat_list[[i]]) <- rna_assay
    sample_seurat_list[[i]] <- NormalizeData(sample_seurat_list[[i]], verbose = FALSE, assay = rna_assay)
    sample_seurat_list[[i]] <- CellCycleScoring(sample_seurat_list[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident = F) 
    DefaultAssay(sample_seurat_list[[i]]) <- rna_assay
  }
  
  return(sample_seurat_list)
}

seurat_list_separated_SoupX <- add_cell_cycle_score_2(seurat_list_separated_SoupX, rna_assay = "SoupX_RNA")


#-----------------------------------------#
## Re-normalize with cell cycle scores ----
#-----------------------------------------#
for (i in 1:length(seurat_list_separated_SoupX)) {
  DefaultAssay(seurat_list_separated_SoupX[[i]]) = "SoupX_RNA"
  seurat_list_separated_SoupX[[i]] = SCTransform(seurat_list_separated_SoupX[[i]], 
                                                 method = "glmGamPoi", 
                                                 vars.to.regress = c("S.Score", "G2M.Score"),
                                                 vst.flavor = "v2",
                                                 verbose = T)
}


#==============================================================================#
# Run DoubletFinder by batch ----
#==============================================================================#
# merge and re split the object by batch instead of each unique sample (UPN,cycle,day,sample type, batch)
seurat_list_separated_SoupX_merge <- merge(x = seurat_list_separated_SoupX[[1]], y = seurat_list_separated_SoupX[2:length(seurat_list_separated_SoupX)])
seurat_list_separated_SoupX_batch <- Seurat::SplitObject(seurat_list_separated_SoupX_merge, split.by = "Batch")


dblrate <- 0.15 # doublet rate

for (i in 1:length(seurat_list_separated_SoupX_batch)) {
  print(paste("Starting analysis on ",names(seurat_list_separated_SoupX_batch)[i], sep = ""))
  # Pre-process Seurat object (sctransform)
  DefaultAssay(seurat_list_separated_SoupX_batch[[i]]) <- "SoupX_RNA"
  seurat_list_separated_SoupX_batch[[i]] <- SCTransform(seurat_list_separated_SoupX_batch[[i]], 
                                                        assay = "SoupX_RNA", 
                                                        method = "glmGamPoi", 
                                                        vars.to.regress = c("S.Score", "G2M.Score"),
                                                        vst.flavor = "v2",
                                                        verbose = F)
  seurat_list_separated_SoupX_batch[[i]] <- RunPCA(seurat_list_separated_SoupX_batch[[i]])
  npcs <- min(get_pcs(seurat_list_separated_SoupX_batch[[i]]))
  seurat_list_separated_SoupX_batch[[i]] <- RunUMAP(seurat_list_separated_SoupX_batch[[i]], dims = 1:npcs)
  seurat_list_separated_SoupX_batch[[i]] <- FindNeighbors(seurat_list_separated_SoupX_batch[[i]], dims = 1:npcs, verbose = F)
  seurat_list_separated_SoupX_batch[[i]] <- FindClusters(seurat_list_separated_SoupX_batch[[i]], resolution = 1, verbose = F)
  
  
  # pK Identification (no ground-truth) 
  # For some reason, after running the preprocessing steps RNA and SoupX RNA 
  # gets split, not sure why this is happening but I am joining it all here
  # before running steps for doubletfinder
  DefaultAssay(seurat_list_separated_SoupX_batch[[i]]) <- "SoupX_RNA"
  seurat_list_separated_SoupX_batch[[i]] <- JoinLayers(seurat_list_separated_SoupX_batch[[i]])
  DefaultAssay(seurat_list_separated_SoupX_batch[[i]]) <- "RNA"
  seurat_list_separated_SoupX_batch[[i]] <- JoinLayers(seurat_list_separated_SoupX_batch[[i]])
  DefaultAssay(seurat_list_separated_SoupX_batch[[i]]) <- "SCT"
  sweep.res.list_seu_test <- paramSweep(seurat_list_separated_SoupX_batch[[i]], PCs = 1:npcs, sct = TRUE)
  sweep.stats_seu_test <- summarizeSweep(sweep.res.list_seu_test, GT = FALSE)
  bcmvn_seu_test <- find.pK(sweep.stats_seu_test)
  
  # Select optimal pK 
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn_seu_test[which.max(bcmvn_seu_test$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  
  # Homotypic Doublet Proportion Estimate 
  homotypic.prop <- modelHomotypic(seurat_list_separated_SoupX_batch[[i]]@meta.data$seurat_clusters)
  nExp_poi <- round(dblrate*nrow(seurat_list_separated_SoupX_batch[[i]]@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  
  # Run DoubletFinder
  seurat_list_separated_SoupX_batch[[i]] <- doubletFinder(seu = seurat_list_separated_SoupX_batch[[i]], 
                                                          PCs = 1:npcs,
                                                          pK = optimal.pk, 
                                                          nExp = nExp_poi, 
                                                          sct = TRUE
  )
  # Run DoubletFinder again w/ pANN
  pANN_col_pos <- ncol(seurat_list_separated_SoupX_batch[[i]]@meta.data) - 1
  seurat_list_separated_SoupX_batch[[i]] <- doubletFinder(seu = seurat_list_separated_SoupX_batch[[i]], 
                                                          PCs = 1:npcs, 
                                                          pK = optimal.pk, 
                                                          nExp = nExp_poi.adj, 
                                                          reuse.pANN = colnames(seurat_list_separated_SoupX_batch[[i]]@meta.data)[pANN_col_pos],
                                                          sct = TRUE
  )
}

# Change colnames
for (i in 1:length(seurat_list_separated_SoupX_batch)) {
  pANN_col_pos <- ncol(seurat_list_separated_SoupX_batch[[i]]@meta.data) - 2
  names(seurat_list_separated_SoupX_batch[[i]]@meta.data)[pANN_col_pos] <- "pANN_col"
  names(seurat_list_separated_SoupX_batch[[i]]@meta.data)[pANN_col_pos+1] <- "DF.class_col"
  names(seurat_list_separated_SoupX_batch[[i]]@meta.data)[pANN_col_pos+2] <- "doublet_finder"
  
}



for (i in 1:length(seurat_list_separated_SoupX_batch)) {
  DefaultAssay(seurat_list_separated_SoupX_batch[[i]]) = "SoupX_RNA"
  seurat_list_separated_SoupX_batch[[i]] = SCTransform(seurat_list_separated_SoupX_batch[[i]], 
                                                       assay = "SoupX_RNA", 
                                                       method = "glmGamPoi", 
                                                       vars.to.regress = c("S.Score", "G2M.Score"),
                                                       vst.flavor = "v2",
                                                       verbose = T)
}


# then re-merge and re-split by how I will integrate
seurat_list_separated_SoupX_batch_merge <- merge(x = seurat_list_separated_SoupX_batch[[1]], 
                                                 y = seurat_list_separated_SoupX_batch[2:length(seurat_list_separated_SoupX_batch)]#,
                                                 #add.cell.ids = names(seurat_list_separated_SoupX_batch)
)

seurat_list_separated_SoupX_batch_merge@meta.data$UPN_Cycle_Day_Sample_Type_Batch <- paste(
  seurat_list_separated_SoupX_batch_merge@meta.data$UPN,
  #seurat_list_separated_SoupX_batchh_merge@meta.data$Cycle, 
  #seurat_list_separated_SoupX_batchh_merge@meta.data$Day, 
  "NA_", 
  "NA_",
  seurat_list_separated_SoupX_batch_merge@meta.data$Sample_Type, 
  seurat_list_separated_SoupX_batch_merge@meta.data$Batch, 
  sep = "_")
seurat_list_separated_SoupX_batch_DoubletFinder <- Seurat::SplitObject(seurat_list_separated_SoupX_batch_merge, split.by = "UPN_Cycle_Day_Sample_Type_Batch")


#==============================================================================#
# Remove RB and MT genes ----
#==============================================================================#
# No doublet finder, new cell cycle scoring
seurat_list_separated_SoupX_batch_DoubletFinder_noRBSMT <- list()
message("Removing ribosomal and mitochondrial genes from Seurat object")
for (i in 1:length(seurat_list_separated_SoupX_batch_DoubletFinder)) {
  RBMTgenes <- grep(pattern = "^RP[SL]|^MRP[SL]|^MT-", x = rownames(seurat_list_separated_SoupX_batch_DoubletFinder[[i]]@assays$SoupX_RNA$counts), value = TRUE, invert = TRUE)
  seurat_list_separated_SoupX_batch_DoubletFinder_noRBSMT[[i]] = DietSeurat(seurat_list_separated_SoupX_batch_DoubletFinder[[i]], features = RBMTgenes, counts = TRUE, data = FALSE, scale.data = FALSE)
  
}
names(seurat_list_separated_SoupX_batch_DoubletFinder_noRBSMT) <- names(seurat_list_separated_SoupX_batch_DoubletFinder)


#==============================================================================#
# Integration, all cells ----
#==============================================================================#

#------------------------------------------------------------------------------#
## Integrate Tumor ----
#------------------------------------------------------------------------------#
# https://satijalab.org/seurat/articles/seurat5_integration
tumor_merge <- merge(x = seurat_list_separated_SoupX_batch_DoubletFinder_noRBSMT[[1]], 
                     y = seurat_list_separated_SoupX_batch_DoubletFinder_noRBSMT[2:length(seurat_list_separated_SoupX_batch_DoubletFinder_noRBSMT)])

# REMOVE DOUBLETS HERE
tumor_merge_singlets <- subset(tumor_merge, subset = doublet_finder == "Singlet")
table(tumor_merge@meta.data$doublet_finder)
table(tumor_merge_singlets@meta.data$doublet_finder)
tumor_merge <- tumor_merge_singlets


DefaultAssay(tumor_merge) <- "SoupX_RNA"
tumor_merge <- JoinLayers(tumor_merge)
tumor_merge

# Split layers by RUN
tumor_merge[["SoupX_RNA"]] <- split(tumor_merge[["SoupX_RNA"]], f = tumor_merge$UPN_Cycle_Day_Sample_Type_Batch)
tumor_merge


Sys.time()

# SCTransform is performed per sample, since they are split into layers
tumor_merge <- SCTransform(tumor_merge, assay = "SoupX_RNA", 
                                         vst.flavor = "v2",
                                         vars.to.regress = c("S.Score", "G2M.Score"),
                                         return.only.var.genes = FALSE,
                                         variable.features.n = 1000)
tumor_merge <- RunPCA(tumor_merge)
npcs <- min(get_pcs(tumor_merge))
npcs # 18 (1000) 
Sys.time()

# Integrate the layers using the SCT values
tumor_merge <- IntegrateLayers(object = tumor_merge,
                               method = RPCAIntegration,
                               assay = "SCT", # either specify here or run default assay to SCT
                               orig.reduction = "pca",
                               new.reduction = "integrated.rpca",
                               verbose = FALSE,
                               normalization.method = "SCT",
                               dims = 1:npcs)

# # re-join layers after integration
# ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]]) # do I need to do this

# Find neighbors and clusters, and create UMAP
tumor_merge <- FindNeighbors(tumor_merge, dims = 1:npcs, verbose = FALSE, reduction = "integrated.rpca")
tumor_merge <- RunUMAP(tumor_merge, dims = 1:npcs, verbose = FALSE,
                       reduction = "integrated.rpca", reduction.name = "umap.integrated.rpca")
tumor_merge <- FindClusters(tumor_merge, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), verbose = FALSE)


saveRDS(tumor_merge, "/scratch/aoill/projects/CAR-T/00_new_2025/geo_2026/tumor_seurat_obj.rds")


#==============================================================================#
# Integration for Immune cells ----
#==============================================================================#
Idents(tumor_merge) <- "SCT_snn_res.0.3"
# get number of cells in each cluster
table(Idents(tumor_merge))

## Merging low cell count samples 

# Get cell names from the clusters that I think are immune cells
tumor_merge_immune <- subset(tumor_merge, idents = c(4, 7, 9, 11))
table(tumor_merge_immune$SCT_snn_res.0.3)

table(tumor_merge_immune@meta.data$UPN_Cycle_Day_Sample_Type_Batch)

#Using Batch for integration due to low cell count in some samples
table(tumor_merge_immune@meta.data$Batch)

# Drop levels
tumor_merge_immune@meta.data[["SCT_snn_res.0.3"]] <- as.factor(tumor_merge_immune@meta.data[["SCT_snn_res.0.3"]])
tumor_merge_immune@meta.data[["SCT_snn_res.0.3"]] <- droplevels(tumor_merge_immune@meta.data[["SCT_snn_res.0.3"]])

DefaultAssay(tumor_merge_immune) <- "SoupX_RNA"
tumor_merge_immune <- JoinLayers(tumor_merge_immune)
tumor_merge_immune

# Split layers by RUN
tumor_merge_immune[["SoupX_RNA"]] <- split(tumor_merge_immune[["SoupX_RNA"]], f = tumor_merge_immune$Batch)
tumor_merge_immune

# SCTransform is performed per sample, since they are split into layers
tumor_merge_immune <- SCTransform(tumor_merge_immune, assay = "SoupX_RNA", 
                                  vst.flavor = "v2",
                                  vars.to.regress = c("S.Score", "G2M.Score"),
                                  return.only.var.genes = FALSE,
                                  variable.features.n = 1000)

tumor_merge_immune <- RunPCA(tumor_merge_immune)
npcs <- min(get_pcs(tumor_merge_immune))
npcs # 12 (1000)
Sys.time()

#table(pbmc_merge_lym@meta.data$UPN_Cycle_Day_Sample_Type_Batch)
# Integrate the layers using the SCT values
tumor_merge_immune <- IntegrateLayers(object = tumor_merge_immune, 
                                      method = RPCAIntegration,
                                      assay = "SCT", # either specify here or run default assay to SCT
                                      orig.reduction = "pca", 
                                      new.reduction = "integrated.rpca",
                                      verbose = F,
                                      normalization.method = "SCT", 
                                      dims = 1:npcs
)

tumor_merge_immune <- FindNeighbors(tumor_merge_immune, dims = 1:npcs, verbose = FALSE, reduction = "integrated.rpca")
tumor_merge_immune <- RunUMAP(tumor_merge_immune, dims = 1:npcs, verbose = FALSE, 
                              reduction = "integrated.rpca", reduction.name = "umap.integrated.rpca")
tumor_merge_immune <- FindClusters(tumor_merge_immune, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), verbose = FALSE)

#------------------------------------------------------------------------------#
## Add CART information to Immune object ----
#------------------------------------------------------------------------------#
DefaultAssay(tumor_merge_immune) <- "RNA"
tumor_merge_immune <- JoinLayers(tumor_merge_immune)

tumor_merge_immune@meta.data$IL13OPCounts <- tumor_merge_immune[["RNA"]]$counts["IL13OP",]
tumor_merge_immune@meta.data$CART <- ifelse(
  (tumor_merge_immune@meta.data$IL13OPCounts >= 3), "Positive", "Negative")

table(tumor_merge_immune@meta.data$CART)

DimPlot(tumor_merge_immune, reduction = "umap.integrated.rpca", group.by = "CART", cols = c("gray", "red"))


#==============================================================================#
# Add TCR data to T cell object ----
#==============================================================================#
## Add TCR data determined from scRepertoire ----
DefaultAssay(tumor_merge_immune) <- "SCT"

## Load the contig data ----
## Tumor
# Three samples, one doesn't have TCR data
F05627 <- read.csv("/tgen_labs/banovich/SingleCell/CellRanger/5_0_0/Ensemble_98/PipelineData/Projects/BCTCSF/CellRangerOuts/VDJ/BCTCSF_0099_1_BR_Whole_T7_X5TCR_F05627_HNJG2DSX5/outs/filtered_contig_annotations.csv")
F05630 <- read.csv("/tgen_labs/banovich/SingleCell/CellRanger/5_0_0/Ensemble_98/PipelineData/Projects/BCTCSF/CellRangerOuts/VDJ/BCTCSF_0115_1_BR_Whole_T5_X5TCR_F05630_HNJG2DSX5/outs/filtered_contig_annotations.csv")


## Process the contig data ----
seurat_obj <- tumor_merge_immune

seurat_obj <- subset(seurat_obj, subset = FID_GEXFB != "F05918-GEX_F05922-CAR")

seurat_list <- SplitObject(seurat_obj, split.by = "FID_GEXFB")

# Strip the barcode extra labeling in seurat_list since the imported contig's 
# cell IDs are not formatted with a prefix (like Batch37_BARCODEID)
seurat_list_cells_renamed <- list()
for (i in names(seurat_list)){
  print(i)
  new_cell_IDs_tmp <- t(as.data.frame(str_split(rownames(seurat_list[[i]]@meta.data), "_")))
  last_col_num <- max(ncol(new_cell_IDs_tmp))
  new_cell_IDs <-  as.character(new_cell_IDs_tmp[,last_col_num])
  renamed_assay <- RenameCells(
    seurat_list[[i]],
    new.names = new_cell_IDs
  )
  seurat_list_cells_renamed <- c(seurat_list_cells_renamed, renamed_assay)
}
names(seurat_list_cells_renamed) <- names(seurat_list)


# Merge all contigs and contig lists into one contig list
contig_list_tumor <- list(F05627, F05630)


names(contig_list_tumor)
names(contig_list_tumor) <- c("F05627_574_Tumor", "F05630_625_Tumor")


## Combine the contigs ----
# Will want to use the parameter removeNA set to true to remove any cell barcode 
# with an NA value in at least one of the chains
combined <- combineTCR(contig_list_tumor, 
                       samples =  c("F05627_574_Tumor", "F05630_625_Tumor"),
                       #cells ="T-AB",
                       filterMulti = FALSE, removeNA = TRUE
)


## Add TCR info to the Seurat objects ----

### Add a new barcode column in the TCR list to match the Seurat objects 
# Replace the barcode column in combined which matching barcode to seurat object 
# Path to google sheet with metadata for your samples
#metapath <- "https://docs.google.com/spreadsheets/d/1fmEuOFXm893mfS2T1zpEsMrN7MpFc_Aym24ZODbw2bA/edit#gid=0"
#sheet_name <- "Sheet1"

sample_metadata <- read.csv("/scratch/rmudunuri/cart_project/19130_all_batches_metadata_new.csv")

# I think it will do the frequency correctly if I add the following information 
# as columns to all entries in this list: UPN, Cycle, Day, Sample_Type,
# UPN_Cycle_Day, UPN_Cycle_Day_Sample_Type
combined2 <- lapply(names(combined), function(i){
  j <- combined[[i]]
  k <- unique(j$sample)
  l <- strsplit(k, "_")[[1]][1] # this is the TCRseq_ID in my metadata file
  # add new barcode column
  # also make an original barcode column (without anything in front of barcode)
  j$tcr_barcode <- j$barcode
  j$orig_barcode <- sapply(strsplit(j$barcode, "_"), `[`, 4) # position in the string where the original barcode starts
  cell_id_prefix <- sample_metadata %>% filter(TCRseq_ID == l) %>% pull(Cell_Prefix) %>% unique()
  j$barcode <- paste0(cell_id_prefix, "_", j$orig_barcode)
  j$TCR_ID <- sapply(strsplit(j$sample, "_"), `[`, 1)
  j$UPN <- sapply(strsplit(j$sample, "_"), `[`, 2)
  j$Cycle <- "NA"
  j$Day <- "NA"
  j$Sample_Type <- sapply(strsplit(j$sample, "_"), `[`, 5)
  j$UPN_Cycle_Day_Sample_Type <- paste(j$UPN, j$Cycle, j$Day, j$Sample_Type, sep = "_")
  j$TCR_ID_UPN_Cycle_Day_Sample_Type <- paste(j$TCR_ID, j$UPN, j$Cycle, j$Day, j$Sample_Type, sep = "_")
  j
})
names(combined2) <- c("F05627_574_Tumor", "F05630_625_Tumor")

saveRDS(combined2, "/scratch/aoill/projects/CAR-T/00_new_2025/geo_2026/tumor_scRepertiore_output.rds")


#==============================================================================#
# Add cell typing info ----
#==============================================================================#
# Cluster 5 looks like tumor so remove from final immune object
tumor_ids <- tumor_merge_immune_TCR@meta.data %>% 
  filter(SCT_snn_res.0.3 == "5") %>%
  rownames()

tumor_imm_notumor <- subset(tumor_merge_immune_TCR, cells = tumor_ids, invert = TRUE)
unique(tumor_imm_notumor@meta.data$ct_subclusters_num)


tumor_imm_notumor@meta.data$ct_final <- tumor_imm_notumor@meta.data$ct_subclusters_num # or any other initialization value
tumor_imm_notumor@meta.data$ct_final[tumor_imm_notumor@meta.data$SCT_snn_res.0.3 == "0"] <- "Mo/Mac"
tumor_imm_notumor@meta.data$ct_final[tumor_imm_notumor@meta.data$SCT_snn_res.0.3 == "2"] <- "Mo/Mac"
tumor_imm_notumor@meta.data$ct_final[tumor_imm_notumor@meta.data$SCT_snn_res.0.3 == "1"] <- "Mo/Mac"
tumor_imm_notumor@meta.data$ct_final[tumor_imm_notumor@meta.data$SCT_snn_res.0.3 == "4"] <- "Mo/Mac"
tumor_imm_notumor@meta.data$ct_final[tumor_imm_notumor@meta.data$SCT_snn_res.0.3 == "8"] <- "Mo/Mac"

tumor_imm_notumor@meta.data$ct_final[tumor_imm_notumor@meta.data$SCT_snn_res.0.3 == "3"] <- "T (Activated)"

tumor_imm_notumor@meta.data$ct_final[tumor_imm_notumor@meta.data$SCT_snn_res.0.3 == "9"] <- "pDC"

tumor_imm_notumor@meta.data$ct_final[tumor_imm_notumor@meta.data$SCT_snn_res.0.3 == "7"] <- "Mast cells"
tumor_imm_notumor@meta.data$ct_final[tumor_imm_notumor@meta.data$SCT_snn_res.0.3 == "6"] <- "cDC2"


#==============================================================================#
# Add additional metadata ----
#==============================================================================#
tumor_merge_immune_TCR@meta.data$tumor_type <- NA # or any other initialization value
tumor_merge_immune_TCR@meta.data$tumor_type[tumor_imm@meta.data$UPN == "574"] <- "Ependymoma"
tumor_merge_immune_TCR@meta.data$tumor_type[tumor_imm@meta.data$UPN == "625"] <- "Ependymoma" 


saveRDS(tumor_merge_immune_TCR, "/scratch/aoill/projects/CAR-T/00_new_2025/geo_2026/tumor_immune_seurat_obj.rds")

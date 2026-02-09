#==============================================================================#
# Author : Angela M. Oill, aoill@tgen.org
# Date: 2025/04/21
# Project: Pediatric CAR-T 
# Description: Pre-processing and integration of PBMC samples
#==============================================================================#


#==============================================================================#
# Load Libraries, helper modules and functions ----
#==============================================================================#
library(Matrix)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(DropletUtils)
library(scRepertoire)
library(scGSVA)
library(RColorBrewer)
library(tidyr)


source("/home/aoill/projects/SingleCellBestPractices/scripts/preprocessing_qc_module.R")
source("/home/aoill/projects/SingleCellBestPractices/scripts/helper_functions_module.R")


#==============================================================================#
# Set variables ----
#==============================================================================#
set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)

# Path to google sheet with metadata for your samples
#metapath <- "https://docs.google.com/spreadsheets/d/1gRL53qgRBApRHxznwTK_17I1cRlAL0GGf8nMikOJle0/edit#gid=0"
metapath <- "/home/aoill/projects/CAR-T/00_new_2025/19130_all_batches_metadata_new.csv"
sheet_name <- "Sheet1"


#==============================================================================#
# Prep metadata ----
#==============================================================================#

#--------------------#
## Read metadata ----
#--------------------#
sample_metadata <- load_metadata(metadata_path = metapath, data_type = "csv", sheet_name = sheet_name)

#-----------------------------------#
## Separate multiplexed samples ----
#-----------------------------------#
# original data
sample_metadata_multiplexed_original <- sample_metadata %>% 
  filter(Multiplexed == "Yes") %>% 
  filter(Batch %in% c(37, 38, 39, 40))

# new data
sample_metadata_multiplexed_new <- sample_metadata %>% 
  filter(Multiplexed == "Yes") %>% 
  filter(Batch %in% c(1, 7, 8, 9, 10))

# no pbmc samples that weren't multiplexed
#sample_metadata_not_multiplexed <- sample_metadata %>% filter(Multiplexed == "No")


#==============================================================================#
# Demultiplex ----
#==============================================================================#

#--------------------------#
## Hash demultiplexing ----
#--------------------------#
seurat_list_demultiplexed <- prep_seurat_list_multiplexed(
  sample_metadata_multiplexed_original, 
  batch_ID = "Cell_Prefix", 
  cellRanger_path = "CellRanger_path", 
  cell_ID_prefix = "Cell_Prefix", 
  CellHashing_Ab = "CellHashing_Ab")


# Keep only singlets
seurat_list_demultiplexed_singlets <- filter_doublets_cellhashing(seurat_list_demultiplexed)

# Keep only PBMC samples
pbmc_seurat_list_demultiplexed_singlets <- list()
for (i in 1:length(seurat_list_demultiplexed_singlets)) {
  print(i)
  sub_i <- subset(seurat_list_demultiplexed_singlets[[i]], subset = Sample_Type == "PBMC")
  pbmc_seurat_list_demultiplexed_singlets <- c(pbmc_seurat_list_demultiplexed_singlets, sub_i)
}
names(pbmc_seurat_list_demultiplexed_singlets) <- c("Batch37_filtered", "Batch38_filtered", "Batch39_filtered", "Batch40_filtered")


pbmc_seurat_list <- pbmc_seurat_list_demultiplexed_singlets


#--------------------#
## New data ----
#--------------------#
# This data is multiplexed and was demultiplexed using demuxlet. I have a file 
# for each batch with which cells are singlets/doublets/ambiguous and which 
# sample it belongs to (demuxlet.best). I will need to join this to the metadata
# based on the cell ID then only keep singlets and split the objects by sample.
# After this proceed like 
metadata = sample_metadata_multiplexed_new
batch_ID = "Cell_Prefix"
cellRanger_path = "CellRanger_path"
cell_ID_prefix = "Cell_Prefix"

sample_list <- unique(metadata[[batch_ID]])


sample_seurat_list <- lapply(sample_list, function(i){
  message(i)
  print(paste("Reading in ", 
              length(unique(paste0(metadata[which(metadata[[batch_ID]]==i),][[cellRanger_path]], "/outs/filtered_feature_bc_matrix/"))), " file.",
              sep = ""
  )
  )
  print(unique(paste0(metadata[which(metadata[[batch_ID]]==i),][[cellRanger_path]], "/outs/filtered_feature_bc_matrix/")))
  sample_10x_data <- Read10X(unique(paste0(metadata[which(metadata[[batch_ID]]==i),][[cellRanger_path]], "/outs/filtered_feature_bc_matrix/")))
  print("Read in 10X file")
  seurat_object <- CreateSeuratObject(counts = sample_10x_data)
  seurat_object <- PercentageFeatureSet(seurat_object, pattern = "^MT-", col.name = "percent.mt_RNA", assay = "RNA")
  seurat_object <- PercentageFeatureSet(seurat_object, pattern = "^RP[SL]|^MRP[SL]", col.name = "percent.ribo_RNA", assay = "RNA")
  seurat_object@meta.data$orig_cellID <- rownames(seurat_object@meta.data)
  print("Added in info to objet")
  
  # add souporcell results here
  soc_fn <- paste("/scratch/aoill/projects/CAR-T/00_new_2025/souporcell/", tolower(i), "_results.tsv", sep = "")
  soc_results <- read.table(soc_fn, header = T)
  print("Read in souporcell results")
  
  
  # Make a BARCODE column in metadata so that I can join the demuxlet results
  seurat_object@meta.data$barcode <- seurat_object@meta.data$orig_cellID
  seurat_object@meta.data$barcode <- sub("^[^_]+_", "", seurat_object@meta.data$barcode)
  
  dim(seurat_object@meta.data)
  
  rn_meta <- rownames(seurat_object@meta.data)
  seurat_object@meta.data <- left_join(seurat_object@meta.data, soc_results)
  rownames(seurat_object@meta.data) <- rn_meta
  dim(seurat_object@meta.data)
  
  # Rename cell IDs, adding a prefix specified by the user
  seurat_object <- RenameCells(seurat_object,
                               new.names = paste0(metadata[which(metadata[[batch_ID]]==i),][[cell_ID_prefix]], 
                                                  "_", colnames(seurat_object)))
  
  # Keep only singlets
  #seurat_object <- subset(seurat_object, subset = DROPLET.TYPE == "SNG")
  seurat_object <- subset(seurat_object, subset = status == "singlet")
  
  # Add a column called UPN from BEST.GUESS column
  # this was from demuxlet outs, doesnt need to be done with the soupor
  seurat_object@meta.data$sample_ID <- sub(",.*", "", seurat_object@meta.data$sample_ID)
  
  
  # Add UPN column based on ID in BEST.GUESS
  upn_map <- c(
    "PEDCAR_0002_1_PB_Whole_C1" = 514,
    "PEDCAR_0003_1_PB_Whole_C1" = 574,
    "PEDCAR_0004_1_PB_Whole_C1" = 689,
    "PEDCAR_0005_1_PB_Whole_C1" = 692,
    "PEDCAR_0006_1_PB_Whole_C1" = 705,
    "PEDCAR_0007_1_PB_Whole_C1" = 716
  )
  
  # Add UPN column to Seurat metadata
  seurat_object@meta.data$UPN <- upn_map[seurat_object@meta.data$sample_ID]
  
  # Add other metadata to object
  batch.meta.data <- metadata[which(metadata[[batch_ID]]==i),]
  rn_s_meta <- rownames(seurat_object@meta.data)
  seurat_object@meta.data <- left_join(seurat_object@meta.data, batch.meta.data)
  rownames(seurat_object@meta.data) <- rn_s_meta
  # then I should have a list of seurat object, one for each batch
  table(seurat_object@meta.data$UPN)
  
  return(seurat_object)
})
names(sample_seurat_list) <- sample_list


# Add together and then look at smooth scatter plots
pbmc_seurat_list_all <- c(pbmc_seurat_list, sample_seurat_list)


#==============================================================================#
# Filter ----
#==============================================================================#

#----------------------------------------#
## Merge Seurat and visualize quality ----
#----------------------------------------#
pbmc_seurat_merge <- merge(x = pbmc_seurat_list_all[[1]], y = pbmc_seurat_list_all[2:length(pbmc_seurat_list_all)])
pbmc_seurat_merge@meta.data$New_Ident <- "SeuratProject"
Idents(pbmc_seurat_merge) <- "New_Ident"

# %mt log(nFeature_RNA)
smoothScatter(pbmc_seurat_merge@meta.data$percent.mt, log(pbmc_seurat_merge@meta.data$nFeature_RNA),
              xlab = "% MT", ylab = "log(nFeature_RNA)",
              main = "PBMCs")
abline(h = log(650), v = 10)
text(15,log(750), "nFeature_RNA = 650,\npercent.mt = 10", adj = c(0, -.1))


# %mt log(nCount_RNA)
smoothScatter(pbmc_seurat_merge@meta.data$percent.mt, log(pbmc_seurat_merge@meta.data$nCount_RNA),
              xlab = "% MT", ylab = "log(nCount_RNA)",
              main = "PBMCs")
abline(h = log(1200), v = 10)
text(15,log(1300), "nCount_RNA = 1200,\npercent.mt = 10", adj = c(0, -.1))

# I think the same filter as before works

#-------------#
## Filter ----
#-------------#
pbmc_seurat_list_filtered <- filter_manual(pbmc_seurat_list_all, 
                                           pt_mt = 10, 
                                           nFeature = 650,
                                           nCount = 1200)

pbmc_seurat_list_separated <- list()
for (i in 1:length(pbmc_seurat_list_filtered)) {
  print(i)
  
  pedi_batch <- pbmc_seurat_list_filtered[[i]]
  pedi_batch@meta.data$Sample_ID <- paste(pedi_batch@meta.data$UPN, 
                                          pedi_batch@meta.data$Sample_Type,
                                          pedi_batch@meta.data$Cycle, 
                                          pedi_batch@meta.data$Day, 
                                          pedi_batch@meta.data$Batch,
                                          sep = "_")
  pedi_batch <- Seurat::SplitObject(pedi_batch, split.by = "Sample_ID")
  pbmc_seurat_list_separated <- c(pbmc_seurat_list_separated, pedi_batch)
}


for (i in 1:length(pbmc_seurat_list_separated)) {
  print(paste0("Converting sample ", names(pbmc_seurat_list_separated[i])))
  obj.sub <- pbmc_seurat_list_separated[[i]]
  DropletUtils::write10xCounts(path = paste0("/scratch/aoill/projects/CAR-T/00_new_2025/soupX/demultiplexed_",names(pbmc_seurat_list_separated[i])), 
                               x = obj.sub[["RNA"]]@layers$counts, 
                               barcodes = colnames(obj.sub[["RNA"]]), # cell names
                               gene.id = rownames(obj.sub[["RNA"]]), # Gene identifiers, one per row of X
                               type = "sparse", version="3", overwrite = T)
}


#-----------------------------------------------------#
## Get number of cells before and after filtering ----
#-----------------------------------------------------#
pbmc_seurat_merge <- merge(x = pbmc_seurat_list_all[[1]], y = pbmc_seurat_list_all[2:length(pbmc_seurat_list_all)])
pbmc_seurat_merge@meta.data$UPN_Cycle_Day_Sample_Type_Batch <- paste0(pbmc_seurat_merge@meta.data$UPN, "_",
                                                                      pbmc_seurat_merge@meta.data$Cycle,  "_", 
                                                                      pbmc_seurat_merge@meta.data$Day, "_",
                                                                      pbmc_seurat_merge@meta.data$Sample_Type, "_",
                                                                      pbmc_seurat_merge@meta.data$Batch)
cell_numbers <- as.data.frame(table(pbmc_seurat_merge@meta.data$UPN_Cycle_Day_Sample_Type_Batch))
colnames(cell_numbers) <- c("Sample", "Unfiltered")

pbmc_seurat_merge <- merge(x = pbmc_seurat_list_filtered[[1]], y = pbmc_seurat_list_filtered[2:length(pbmc_seurat_list_filtered)])
pbmc_seurat_merge@meta.data$UPN_Cycle_Day_Sample_Type_Batch <- paste0(pbmc_seurat_merge@meta.data$UPN, "_",
                                                                      pbmc_seurat_merge@meta.data$Cycle,  "_", 
                                                                      pbmc_seurat_merge@meta.data$Day, "_",
                                                                      pbmc_seurat_merge@meta.data$Sample_Type, "_",
                                                                      pbmc_seurat_merge@meta.data$Batch)
cell_numbers_filtered <- as.data.frame(table(pbmc_seurat_merge@meta.data$UPN_Cycle_Day_Sample_Type_Batch))
colnames(cell_numbers_filtered) <- c("Sample", "Filtered")

cell_numbers_all <- join(cell_numbers, cell_numbers_filtered)
cell_numbers_all$prop_keep <- cell_numbers_all$Filtered/cell_numbers_all$Unfiltered


cell_numbers_all <- cell_numbers_all %>% separate_wider_delim(Sample, delim = "_", names = c("UPN", "Cycle", "Day", "Sample_Type", "Batch"))
sum(cell_numbers_all$Unfiltered)
# 204171
sum(cell_numbers_all$Filtered)
# 191556


# Cell numbers of new samples
cell_numbers_all %>% filter(Batch %in% c("1", "7", "8", "9", "10")) %>% 
  dplyr::select(Unfiltered) %>% sum()
# 172141
cell_numbers_all %>% filter(Batch %in% c("1", "7", "8", "9", "10")) %>% 
  dplyr::select(Filtered) %>% sum()
# 166942

print(cell_numbers_all, n = nrow(cell_numbers_all))


#==============================================================================#
# Run SoupX ----
#==============================================================================#
# Filter out 689_PBMC_8_1_8 before running 
name_to_remove <- "689_PBMC_8_1_8"
pbmc_seurat_list_separated_subset <- pbmc_seurat_list_separated[names(pbmc_seurat_list_separated) != name_to_remove]


# Extract product samples from these above objects
sample_list <- names(pbmc_seurat_list_separated_subset)


pbmc_seurat_list_separated_SoupX <- sapply(sample_list,  function(i){
  print(i)
  # Read in count and droplet data
  d10x_toc <- Read10X(paste0("/scratch/aoill/projects/CAR-T/00_new_2025/soupX/demultiplexed_", i))
  
  # Need to read in batch specific empty droplet file
  upn_id <- str_split(i, "_")[[1]][1]
  sample_type_id <- str_split(i, "_")[[1]][2]
  cycle_type_id <- str_split(i, "_")[[1]][3]
  day_type_id <- str_split(i, "_")[[1]][4]
  batch_id <- str_split(i, "_")[[1]][5]
  if (batch_id %in% c("37", "38", "39", "40")) {
    batch_path <- sample_metadata_multiplexed_original %>% 
      filter(UPN == upn_id) %>%
      filter(Sample_Type == sample_type_id) %>%
      filter(Cycle == cycle_type_id) %>%
      filter(Day == day_type_id) %>%
      filter(Batch == batch_id) %>%
      pull(CellRanger_path) %>% unique()
    d10x_tod <- Read10X(paste0(batch_path, "/outs/raw_feature_bc_matrix/"))
    
    # Run SoupX
    sc <- SoupChannel(d10x_tod[[1]], d10x_toc, calcSoupProfile = FALSE) 
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
    
    
  } else { # this is the new demultiplexed data
    batch_path <- sample_metadata_multiplexed_new %>% filter(Batch == batch_id) %>%
      pull(CellRanger_path) %>% unique()
    d10x_tod <- Read10X(paste0(batch_path, "/outs/raw_feature_bc_matrix/"))
    
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
  }
  
  
  
})
names(pbmc_seurat_list_separated_SoupX) <- sample_list


# Add metadata
#batch_ID = "FID_GEXFB"
cellRanger_path = "CellRanger_path"
for (i in names(pbmc_seurat_list_separated_SoupX)) {
  upn_id <- str_split(i, "_")[[1]][1]
  sample_type_id <- str_split(i, "_")[[1]][2]
  cycle_type_id <- str_split(i, "_")[[1]][3]
  day_type_id <- str_split(i, "_")[[1]][4]
  batch_type_id <- str_split(i, "_")[[1]][5]
  #batch_id <- str_split(i, "_")[[1]][5]
  meta.data <- sample_metadata %>% filter(UPN == upn_id) %>% 
    filter(Sample_Type == sample_type_id) %>% 
    filter(Cycle == cycle_type_id) %>% 
    filter(Day == day_type_id) %>% 
    filter(Batch == batch_type_id)
  
  # Add sample metadata
  for (m in colnames(meta.data)){
    if (m == cellRanger_path) {
      next
    } else {
      pbmc_seurat_list_separated_SoupX[[i]]@meta.data[[m]] = meta.data[[m]]
    }
  }
  pbmc_seurat_list_separated_SoupX[[i]]@meta.data$Sample_ID <- paste(pbmc_seurat_list_separated_SoupX[[i]]@meta.data$UPN, 
                                                                     pbmc_seurat_list_separated_SoupX[[i]]@meta.data$Sample_Type,
                                                                     pbmc_seurat_list_separated_SoupX[[i]]@meta.data$Cycle, 
                                                                     pbmc_seurat_list_separated_SoupX[[i]]@meta.data$Day, 
                                                                     pbmc_seurat_list_separated_SoupX[[i]]@meta.data$Batch,
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
# I do not think I need to run this step
for (i in 1:length(seurat_list_separated_SoupX)) {
  DefaultAssay(seurat_list_separated_SoupX[[i]]) = "SoupX_RNA"
  seurat_list_separated_SoupX[[i]] = SCTransform(seurat_list_separated_SoupX[[i]], 
                                                      method = "glmGamPoi", 
                                                      vars.to.regress = c("S.Score", "G2M.Score"),
                                                      vst.flavor = "v2",
                                                      verbose = T, 
                                                      assay = "SoupX_RNA")
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
  seurat_list_separated_SoupX_batch_merge@meta.data$Cycle, 
  seurat_list_separated_SoupX_batch_merge@meta.data$Day, 
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
# Integrate ----
#==============================================================================#
# https://satijalab.org/seurat/articles/seurat5_integration
pbmc_merge_rpca <- merge(x = seurat_list_separated_SoupX_batch_DoubletFinder_noRBSMT[[1]], 
                        y = seurat_list_separated_SoupX_batch_DoubletFinder_noRBSMT[2:length(seurat_list_separated_SoupX_batch_DoubletFinder_noRBSMT)])

# REMOVE DOUBLETS HERE
pbmc_merge_rpca_singlets <- subset(pbmc_merge_rpca, subset = doublet_finder == "Singlet")
table(pbmc_merge_rpca@meta.data$doublet_finder)
table(pbmc_merge_rpca_singlets@meta.data$doublet_finder)
pbmc_merge_rpca <- pbmc_merge_rpca_singlets

DefaultAssay(pbmc_merge_rpca) <- "SoupX_RNA"
pbmc_merge_rpca <- JoinLayers(pbmc_merge_rpca)
#pbmc_merge_rpca

# Split layers by RUN
pbmc_merge_rpca[["SoupX_RNA"]] <- split(pbmc_merge_rpca[["SoupX_RNA"]], f = pbmc_merge_rpca$UPN_Cycle_Day_Sample_Type_Batch)
#pbmc_merge_rpca


Sys.time()

# SCTransform is performed per sample, since they are split into layers
pbmc_merge_rpca <- SCTransform(pbmc_merge_rpca, assay = "SoupX_RNA", 
                              vst.flavor = "v2",
                              vars.to.regress = c("S.Score", "G2M.Score"),
                              return.only.var.genes = FALSE,
                              variable.features.n = 2500)
pbmc_merge_rpca <- RunPCA(pbmc_merge_rpca)
npcs <- min(get_pcs(pbmc_merge_rpca))
npcs # 18 (2500)
Sys.time()

min(table(pbmc_merge_rpca$UPN_Cycle_Day_Sample_Type_Batch))
# Integrate the layers using the SCT values
pbmc_merge_rpca <- IntegrateLayers(object = pbmc_merge_rpca, 
                                  method = RPCAIntegration,
                                  assay = "SCT", # either specify here or run default assay to SCT
                                  orig.reduction = "pca", 
                                  new.reduction = "integrated.rpca",
                                  verbose = FALSE,
                                  normalization.method = "SCT", 
                                  dims = 1:npcs,
                                  k.weight = 50)

# Find neighbors and clusters, and create UMAP
pbmc_merge_rpca <- FindNeighbors(pbmc_merge_rpca, dims = 1:npcs, verbose = FALSE, reduction = "integrated.rpca")
pbmc_merge_rpca <- RunUMAP(pbmc_merge_rpca, dims = 1:npcs, verbose = FALSE, 
                          reduction = "integrated.rpca", reduction.name = "umap.integrated.rpca")
pbmc_merge_rpca <- FindClusters(pbmc_merge_rpca, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), verbose = FALSE)


#------------------------------------------------------------------------------#
## Add metadata ----
#------------------------------------------------------------------------------#
pbmc_merge_rpca@meta.data$UPN_Cycle <- paste(pbmc_merge_rpca@meta.data$UPN, pbmc_merge_rpca@meta.data$Cycle, sep = "_")


pbmc_merge_rpca@meta.data$Response_cat <- "NA" # or any other initialization value
pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "514_2"] <- "Baseline"
pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "514_5"] <- "Non-response"
pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "514_8"] <- "Response"
pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "514_11"] <- "Non-response"
pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "514_13"] <- "Non-response"

pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "515_4"] <- "Response"
pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "515_8"] <- "Non-response"

pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "574_1"] <- "Baseline"
pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "574_3"] <- "Non-response"
pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "574_8"] <- "Non-response"

pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "689_1"] <- "Baseline"
pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "689_4"] <- "Response"
pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "689_8"] <- "Response"
pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "689_12"] <- "Non-response"
pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "689_17"] <- "Non-response"

pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "692_2"] <- "Baseline"
pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "692_4"] <- "Response"

pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "705_1"] <- "Baseline"
pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "705_5"] <- "Non-response"
pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "705_8"] <- "Non-response"

pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "716_2"] <- "Baseline"
pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "716_4"] <- "Response"
pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "716_8"] <- "Response"

pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "625_1"] <- "Baseline" 
pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "625_3"] <- "Non-response" 
pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "625_4"] <- "Non-response" 

pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "626_1"] <- "Baseline" 
pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "626_4"] <- "Response" 
pbmc_merge_rpca@meta.data$Response_cat[pbmc_merge_rpca@meta.data$UPN_Cycle == "626_7"] <- "Non-response" 

table(pbmc_merge_rpca@meta.data$Response_cat)


# add disease info
pbmc_merge_rpca@meta.data$tumor_type <- NA # or any other initialization value
pbmc_merge_rpca@meta.data$tumor_type[pbmc_merge_rpca@meta.data$UPN == "514"] <- "Ependymoma"
pbmc_merge_rpca@meta.data$tumor_type[pbmc_merge_rpca@meta.data$UPN == "574"] <- "Ependymoma"
pbmc_merge_rpca@meta.data$tumor_type[pbmc_merge_rpca@meta.data$UPN == "625"] <- "Ependymoma" 
pbmc_merge_rpca@meta.data$tumor_type[pbmc_merge_rpca@meta.data$UPN == "689"] <- "Ependymoma"

pbmc_merge_rpca@meta.data$tumor_type[pbmc_merge_rpca@meta.data$UPN == "515"] <- "DMG"
pbmc_merge_rpca@meta.data$tumor_type[pbmc_merge_rpca@meta.data$UPN == "626"] <- "DMG" 
pbmc_merge_rpca@meta.data$tumor_type[pbmc_merge_rpca@meta.data$UPN == "692"] <- "DMG"
pbmc_merge_rpca@meta.data$tumor_type[pbmc_merge_rpca@meta.data$UPN == "716"] <- "DMG"

pbmc_merge_rpca@meta.data$tumor_type[pbmc_merge_rpca@meta.data$UPN == "705"] <- "Glioblastoma"


# Add lymphodepletion status
unique(pbmc_merge_rpca@meta.data$UPN)
#[1] 514 574 515 626 625 705 689 716 692

pbmc_merge_rpca@meta.data$Lymphodepletion <- NA # or any other initialization value
pbmc_merge_rpca@meta.data$Lymphodepletion[pbmc_merge_rpca@meta.data$UPN == "514"] <- "Non-Lymphodepleted"
pbmc_merge_rpca@meta.data$Lymphodepletion[pbmc_merge_rpca@meta.data$UPN == "515"] <- "Non-Lymphodepleted"
pbmc_merge_rpca@meta.data$Lymphodepletion[pbmc_merge_rpca@meta.data$UPN == "574"] <- "Non-Lymphodepleted" 
pbmc_merge_rpca@meta.data$Lymphodepletion[pbmc_merge_rpca@meta.data$UPN == "625"] <- "Lymphodepleted" 
pbmc_merge_rpca@meta.data$Lymphodepletion[pbmc_merge_rpca@meta.data$UPN == "626"] <- "Lymphodepleted" 
pbmc_merge_rpca@meta.data$Lymphodepletion[pbmc_merge_rpca@meta.data$UPN == "689"] <- "Lymphodepleted" 
pbmc_merge_rpca@meta.data$Lymphodepletion[pbmc_merge_rpca@meta.data$UPN == "705"] <- "Lymphodepleted" 
pbmc_merge_rpca@meta.data$Lymphodepletion[pbmc_merge_rpca@meta.data$UPN == "716"] <- "Lymphodepleted" 
pbmc_merge_rpca@meta.data$Lymphodepletion[pbmc_merge_rpca@meta.data$UPN == "692"] <- "Lymphodepleted" 



#------------------------------------------------------------------------------#
## Sub-cluster ----
#------------------------------------------------------------------------------#
# Cluster 6
# Some erthrocyte, some lymphoid
Idents(pbmc_merge_rpca) <- "SCT_snn_res.0.2"
pbmc_merge_sub <- FindSubCluster(pbmc_merge_rpca, "6", resolution = 0.2, graph.name = "SCT_snn")


# Lym/Mye integration ----
#Lymphoid: 3, 4, 6_5 , 6_6, 7, 8, 10
#Myeloid: 0, 1, 2, 5, 9, 12



#==============================================================================#
# Integrate Lymphoid, pass 1 ----
#==============================================================================#
# pbmc_merge_sub
lymphoid <- c("3", "4", "6_5", "6_6", "7", "8", "10")
# we are using sub.cluster 
Idents(pbmc_merge_sub) <- "sub.cluster"
# get number of cells in each cluster
table(Idents(pbmc_merge_sub))
#6_1    10     7    12     4     3   6_0     2     9     8     5   6_4     1     0    11   6_2   6_6   6_5   6_3 
#3569  3886 10157  2475 14068 15474  3928 17274  6101  7268 12187   880 22652 38536  2546  1407   215   704  1289 

# Get cell names from the clusters that I think are lymphoid cells
pbmc_merge_sub_lym <- subset(pbmc_merge_sub, idents = lymphoid)
table(pbmc_merge_sub_lym@meta.data$sub.cluster)
#10     3     4   6_5   6_6     7     8 
#3886 15474 14068   704   215 10157  7268 

# Drop levels
pbmc_merge_sub_lym@meta.data[["sub.cluster"]] <- as.factor(pbmc_merge_sub_lym@meta.data[["sub.cluster"]])
pbmc_merge_sub_lym@meta.data[["sub.cluster"]] <- droplevels(pbmc_merge_sub_lym@meta.data[["sub.cluster"]])

DefaultAssay(pbmc_merge_sub_lym) <- "SoupX_RNA"
pbmc_merge_sub_lym <- JoinLayers(pbmc_merge_sub_lym)
pbmc_merge_sub_lym

# Split layers by RUN
pbmc_merge_sub_lym[["SoupX_RNA"]] <- split(pbmc_merge_sub_lym[["SoupX_RNA"]], f = pbmc_merge_sub_lym$UPN_Cycle_Day_Sample_Type_Batch)
pbmc_merge_sub_lym


# SCTransform is performed per sample, since they are split into layers
pbmc_merge_sub_lym <- SCTransform(pbmc_merge_sub_lym, assay = "SoupX_RNA", 
                                  vst.flavor = "v2",
                                  vars.to.regress = c("S.Score", "G2M.Score"),
                                  return.only.var.genes = FALSE,
                                  variable.features.n = 2500)

pbmc_merge_sub_lym <- RunPCA(pbmc_merge_sub_lym)
npcs <- min(get_pcs(pbmc_merge_sub_lym))
npcs # 13 (2500)
Sys.time()


min(table(pbmc_merge_sub_lym$UPN_Cycle_Day_Sample_Type_Batch))
#[1] 68

# Integrate the layers using the SCT values
pbmc_merge_sub_lym <- IntegrateLayers(object = pbmc_merge_sub_lym, 
                                      method = RPCAIntegration,
                                      assay = "SCT", # either specify here or run default assay to SCT
                                      orig.reduction = "pca", 
                                      new.reduction = "integrated.rpca",
                                      verbose = FALSE,
                                      normalization.method = "SCT", 
                                      dims = 1:npcs,
                                      k.weight = 45
)

# Error: k.weight (100) is set larger than the number of cells in the smallest object (60). Please choose a smaller k.weight.


# Find neighbors and clusters, and create UMAP
pbmc_merge_sub_lym <- FindNeighbors(pbmc_merge_sub_lym, dims = 1:npcs, verbose = FALSE, reduction = "integrated.rpca")
pbmc_merge_sub_lym <- RunUMAP(pbmc_merge_sub_lym, dims = 1:npcs, verbose = FALSE, 
                              reduction = "integrated.rpca", reduction.name = "umap.integrated.rpca")
pbmc_merge_sub_lym <- FindClusters(pbmc_merge_sub_lym, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), verbose = FALSE)


#------------------------------------------------------------------------------#
## Sub-cluster ----
#------------------------------------------------------------------------------#
# Cluster 5
Idents(pbmc_merge_sub_lym) <- "SCT_snn_res.0.2"
pbmc_merge_sub_lym_sub <- FindSubCluster(pbmc_merge_sub_lym, "5", resolution = 0.15, graph.name = "SCT_snn")


#Lymphoid: 0, 1, 2, 3, 5_2, 5_3, 6
#Myeloid: 5_4
#Doublets/low quality: 4, 5_0, 5_1, 5_5


#==============================================================================#
# Integrate Myeloid, pass 1 ----
#==============================================================================#
# Myeloid with 5_4 from lymphoid pass 1 object
# pbmc_merge_sub_lym sub.cluster
# Myeliod clusters from full object pbmc_merge_sub
myeloid <- c("0", "1", "2", "5", "9", "12")
pbmc_merge_sub@meta.data$cell_ID <- rownames(pbmc_merge_sub@meta.data)
myeloid_IDs_full_obj <- pbmc_merge_sub@meta.data %>% filter(sub.cluster %in% myeloid) %>% 
  dplyr::pull(cell_ID) 

# Get cell IDs from lymphoid pass 1 object
# pbmc_merge_sub_lym_sub sub.cluster 5_4
pbmc_merge_sub_lym_sub@meta.data$cell_ID <- rownames(pbmc_merge_sub_lym_sub@meta.data)
myeloid_IDs_lym_obj <- pbmc_merge_sub_lym_sub@meta.data %>% filter(sub.cluster %in% c("5_4")) %>% 
  dplyr::pull(cell_ID) 

myeloid_IDs_all <- c(myeloid_IDs_full_obj, myeloid_IDs_lym_obj)
length(myeloid_IDs_all)
# 99514
# Get cell names from the clusters that I think are myeloid cells
pbmc_merge_sub_mye <- subset(pbmc_merge_sub, subset = cell_ID %in% myeloid_IDs_all)
nrow(pbmc_merge_sub_mye@meta.data)
# 99514

DefaultAssay(pbmc_merge_sub_mye) <- "SoupX_RNA"
pbmc_merge_sub_mye <- JoinLayers(pbmc_merge_sub_mye)
#pbmc_merge_sub_mye

# Split layers by RUN
pbmc_merge_sub_mye[["SoupX_RNA"]] <- split(pbmc_merge_sub_mye[["SoupX_RNA"]], f = pbmc_merge_sub_mye$UPN_Cycle_Day_Sample_Type_Batch)
#pbmc_merge_sub_mye


# SCTransform is performed per sample, since they are split into layers
pbmc_merge_sub_mye <- SCTransform(pbmc_merge_sub_mye, assay = "SoupX_RNA", 
                                  vst.flavor = "v2",
                                  vars.to.regress = c("S.Score", "G2M.Score"),
                                  return.only.var.genes = FALSE,
                                  variable.features.n = 2500)

pbmc_merge_sub_mye <- RunPCA(pbmc_merge_sub_mye)
npcs <- min(get_pcs(pbmc_merge_sub_mye))
npcs # 15 (2500)
Sys.time()

min(table(pbmc_merge_sub_mye$UPN_Cycle_Day_Sample_Type_Batch))
#[1] 32

# Integrate the layers using the SCT values
pbmc_merge_sub_mye <- IntegrateLayers(object = pbmc_merge_sub_mye, 
                                      method = RPCAIntegration,
                                      assay = "SCT", # either specify here or run default assay to SCT
                                      orig.reduction = "pca", 
                                      new.reduction = "integrated.rpca",
                                      verbose = FALSE,
                                      normalization.method = "SCT", 
                                      dims = 1:npcs,
                                      k.weight = 25
)

# Error: k.weight (100) is set larger than the number of cells in the smallest object (60). Please choose a smaller k.weight.


# Find neighbors and clusters, and create UMAP
pbmc_merge_sub_mye <- FindNeighbors(pbmc_merge_sub_mye, dims = 1:npcs, verbose = FALSE, reduction = "integrated.rpca")
pbmc_merge_sub_mye <- RunUMAP(pbmc_merge_sub_mye, dims = 1:npcs, verbose = FALSE, 
                              reduction = "integrated.rpca", reduction.name = "umap.integrated.rpca")
pbmc_merge_sub_mye <- FindClusters(pbmc_merge_sub_mye, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), verbose = FALSE)

# clust 6 looks like T cells 


#==============================================================================#
# Integrate Lymphoid, pass 2 ----
#==============================================================================#
# Lymphoid pass 2: Add cluster 6 from pass 1 myeloid object and re-integrate
# Get cell IDs from lymphoid pass 1 object
# pbmc_merge_sub_lym_sub sub.cluster 5_4
pbmc_merge_sub_lym_sub@meta.data$cell_ID <- rownames(pbmc_merge_sub_lym_sub@meta.data)
lymphoid_IDs_lym_obj <- pbmc_merge_sub_lym_sub@meta.data %>% 
  filter(sub.cluster %in% c("0", "1", "2", "3", "5_2", "5_3", "6")) %>% 
  dplyr::pull(cell_ID) 
length(lymphoid_IDs_lym_obj)
# 43390

pbmc_merge_sub_mye@meta.data$cell_ID <- rownames(pbmc_merge_sub_mye@meta.data)
lymphoid_IDs_mye_obj <- pbmc_merge_sub_mye@meta.data %>% 
  filter(SCT_snn_res.0.2 %in% c("6")) %>% 
  dplyr::pull(cell_ID) 
length(lymphoid_IDs_mye_obj)
# 152

# Get cell IDs from myeloid pass 1 object
lymphoid_IDs_all <- c(lymphoid_IDs_lym_obj, lymphoid_IDs_mye_obj)
length(lymphoid_IDs_all)
# 43542
# Get cell names from the clusters that I think are myeloid cells
pbmc_merge_sub_lym_pass2 <- subset(pbmc_merge_sub, subset = cell_ID %in% lymphoid_IDs_all)
nrow(pbmc_merge_sub_lym_pass2@meta.data)
# 43542

DefaultAssay(pbmc_merge_sub_lym_pass2) <- "SoupX_RNA"
pbmc_merge_sub_lym_pass2 <- JoinLayers(pbmc_merge_sub_lym_pass2)
#pbmc_merge_sub_lym_pass2

# Split layers by RUN
pbmc_merge_sub_lym_pass2[["SoupX_RNA"]] <- split(pbmc_merge_sub_lym_pass2[["SoupX_RNA"]], f = pbmc_merge_sub_lym_pass2$UPN_Cycle_Day_Sample_Type_Batch)


# SCTransform is performed per sample, since they are split into layers
pbmc_merge_sub_lym_pass2 <- SCTransform(pbmc_merge_sub_lym_pass2, assay = "SoupX_RNA", 
                                        vst.flavor = "v2",
                                        vars.to.regress = c("S.Score", "G2M.Score"),
                                        return.only.var.genes = FALSE,
                                        variable.features.n = 2500)

pbmc_merge_sub_lym_pass2 <- RunPCA(pbmc_merge_sub_lym_pass2)
npcs <- min(get_pcs(pbmc_merge_sub_lym_pass2))
npcs # 15 (2500)
Sys.time()


min(table(pbmc_merge_sub_lym_pass2$UPN_Cycle_Day_Sample_Type_Batch))
#[1] 44

# Integrate the layers using the SCT values
pbmc_merge_sub_lym_pass2 <- IntegrateLayers(object = pbmc_merge_sub_lym_pass2, 
                                            method = RPCAIntegration,
                                            assay = "SCT", # either specify here or run default assay to SCT
                                            orig.reduction = "pca", 
                                            new.reduction = "integrated.rpca",
                                            verbose = FALSE,
                                            normalization.method = "SCT", 
                                            dims = 1:npcs,
                                            k.weight = 35
)

# Error: k.weight (100) is set larger than the number of cells in the smallest object (60). Please choose a smaller k.weight.

# Find neighbors and clusters, and create UMAP
pbmc_merge_sub_lym_pass2 <- FindNeighbors(pbmc_merge_sub_lym_pass2, dims = 1:npcs, verbose = FALSE, reduction = "integrated.rpca")
pbmc_merge_sub_lym_pass2 <- RunUMAP(pbmc_merge_sub_lym_pass2, dims = 1:npcs, verbose = FALSE, 
                                    reduction = "integrated.rpca", reduction.name = "umap.integrated.rpca")
pbmc_merge_sub_lym_pass2 <- FindClusters(pbmc_merge_sub_lym_pass2, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), verbose = FALSE)



#------------------------------------------------------------------------------#
## Sub-cluster ----
#------------------------------------------------------------------------------#
# remove cluster 6 and re-integrate 
# Cluster 4
Idents(pbmc_merge_sub_lym_pass2) <- "SCT_snn_res.0.2"
pbmc_merge_sub_lym_pass2_sub <- FindSubCluster(pbmc_merge_sub_lym_pass2, "4", resolution = 0.2, graph.name = "SCT_snn")

#4_0: B cells
#4_1: doublet
#4_2: pDC
#4_3: doublet


#==============================================================================#
# Integrate Lymphoid, pass 3 ----
#==============================================================================#
# Lymphoid pass 3, remove doublet clusters 
# pbmc_merge_sub_lym_pass2_sub
#4_0: B cells
#4_1: doublet
#4_2: pDC
#4_3: doublet
lymphoid_pass3 <- c("0", "1", "2", "3", "4_0", "4_2", "5")
# we are using sub.cluster 
Idents(pbmc_merge_sub_lym_pass2_sub) <- "sub.cluster"
# get number of cells in each cluster
table(Idents(pbmc_merge_sub_lym_pass2_sub))
#4_0     0     2   4_1     3     1     6     5   4_2   4_3 
#2634 11807  9116   465  7075 10110   573  1088   415   259 

# Get cell names from the clusters that I think are lymphoid cells
pbmc_merge_lym_pass3 <- subset(pbmc_merge_sub_lym_pass2_sub, idents = lymphoid_pass3)
table(pbmc_merge_lym_pass3@meta.data$sub.cluster)
#0     1     2     3   4_0   4_2     5 
#11807 10110  9116  7075  2634   415  1088 

# Drop levels
pbmc_merge_lym_pass3@meta.data[["sub.cluster"]] <- as.factor(pbmc_merge_lym_pass3@meta.data[["sub.cluster"]])
pbmc_merge_lym_pass3@meta.data[["sub.cluster"]] <- droplevels(pbmc_merge_lym_pass3@meta.data[["sub.cluster"]])

DefaultAssay(pbmc_merge_lym_pass3) <- "SoupX_RNA"
pbmc_merge_lym_pass3 <- JoinLayers(pbmc_merge_lym_pass3)
#pbmc_merge_lym_pass3

# Split layers by RUN
pbmc_merge_lym_pass3[["SoupX_RNA"]] <- split(pbmc_merge_lym_pass3[["SoupX_RNA"]], f = pbmc_merge_lym_pass3$UPN_Cycle_Day_Sample_Type_Batch)
#pbmc_merge_lym_pass3


# SCTransform is performed per sample, since they are split into layers
pbmc_merge_lym_pass3 <- SCTransform(pbmc_merge_lym_pass3, assay = "SoupX_RNA", 
                                    vst.flavor = "v2",
                                    vars.to.regress = c("S.Score", "G2M.Score"),
                                    return.only.var.genes = FALSE,
                                    variable.features.n = 2500)

pbmc_merge_lym_pass3 <- RunPCA(pbmc_merge_lym_pass3)
npcs <- min(get_pcs(pbmc_merge_lym_pass3))
npcs # 14 (2500)
Sys.time()


min(table(pbmc_merge_lym_pass3$UPN_Cycle_Day_Sample_Type_Batch))
#[1] 42

# Integrate the layers using the SCT values
pbmc_merge_lym_pass3 <- IntegrateLayers(object = pbmc_merge_lym_pass3, 
                                        method = RPCAIntegration,
                                        assay = "SCT", # either specify here or run default assay to SCT
                                        orig.reduction = "pca", 
                                        new.reduction = "integrated.rpca",
                                        verbose = FALSE,
                                        normalization.method = "SCT", 
                                        dims = 1:npcs,
                                        k.weight = 35
)



# Find neighbors and clusters, and create UMAP
pbmc_merge_lym_pass3 <- FindNeighbors(pbmc_merge_lym_pass3, dims = 1:npcs, verbose = FALSE, reduction = "integrated.rpca")
pbmc_merge_lym_pass3 <- RunUMAP(pbmc_merge_lym_pass3, dims = 1:npcs, verbose = FALSE, 
                                reduction = "integrated.rpca", reduction.name = "umap.integrated.rpca")
pbmc_merge_lym_pass3 <- FindClusters(pbmc_merge_lym_pass3, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), verbose = FALSE)


#------------------------------------------------------------------------------#
## Sub-cluster ----
#------------------------------------------------------------------------------#
# Cluster 5
# subcluster 5 to separate pDCs from b cells 
Idents(pbmc_merge_lym_pass3_5) <- "SCT_snn_res.0.2"
pbmc_merge_lym_pass3_5_sub <- FindSubCluster(pbmc_merge_lym_pass3_5, "5", resolution = 0.26, graph.name = "SCT_snn")

#5_0: B cells
#5_1: B cells
#5_2: doublets
#5_3: pdc


#==============================================================================#
# Integrate Lymphoid, pass 4 ----
#==============================================================================#
# Lymphoid pass 4 remove doublet cluster
#5_0: B cells
#5_1: B cells
#5_2: doublets
#5_3: pdc
lymphoid_pass4 <- c("0", "1", "2", "3", "4", "5_0", "5_1", "5_3", "6")
# we are using sub.cluster 
Idents(pbmc_merge_lym_pass3_sub) <- "sub.cluster"
# get number of cells in each cluster
table(Idents(pbmc_merge_lym_pass3_sub))
#5_0     4     1     0   5_1     3   5_2     2     6   5_3 
#1297  4781  8845 10482  1040  6389   474  7534  1154   249 

# Get cell names from the clusters that I think are lymphoid cells
pbmc_merge_lym_pass4 <- subset(pbmc_merge_lym_pass3_sub, idents = lymphoid_pass4)
table(pbmc_merge_lym_pass4@meta.data$sub.cluster)
#0     1     2     3     4   5_0   5_1   5_3     6 
#10482  8845  7534  6389  4781  1297  1040   249  1154 

# Drop levels
pbmc_merge_lym_pass4@meta.data[["sub.cluster"]] <- as.factor(pbmc_merge_lym_pass4@meta.data[["sub.cluster"]])
pbmc_merge_lym_pass4@meta.data[["sub.cluster"]] <- droplevels(pbmc_merge_lym_pass4@meta.data[["sub.cluster"]])

DefaultAssay(pbmc_merge_lym_pass4) <- "SoupX_RNA"
pbmc_merge_lym_pass4 <- JoinLayers(pbmc_merge_lym_pass4)
#pbmc_merge_lym_pass4

# Split layers by RUN
pbmc_merge_lym_pass4[["SoupX_RNA"]] <- split(pbmc_merge_lym_pass4[["SoupX_RNA"]], f = pbmc_merge_lym_pass4$UPN_Cycle_Day_Sample_Type_Batch)
#pbmc_merge_lym_pass4


# SCTransform is performed per sample, since they are split into layers
pbmc_merge_lym_pass4 <- SCTransform(pbmc_merge_lym_pass4, assay = "SoupX_RNA", 
                                    vst.flavor = "v2",
                                    vars.to.regress = c("S.Score", "G2M.Score"),
                                    return.only.var.genes = FALSE,
                                    variable.features.n = 2500)

pbmc_merge_lym_pass4 <- RunPCA(pbmc_merge_lym_pass4)
npcs <- min(get_pcs(pbmc_merge_lym_pass4))
npcs # 14 (2500)
Sys.time()

min(table(pbmc_merge_lym_pass4$UPN_Cycle_Day_Sample_Type_Batch))
#[1] 42

# Integrate the layers using the SCT values
pbmc_merge_lym_pass4 <- IntegrateLayers(object = pbmc_merge_lym_pass4, 
                                        method = RPCAIntegration,
                                        assay = "SCT", # either specify here or run default assay to SCT
                                        orig.reduction = "pca", 
                                        new.reduction = "integrated.rpca",
                                        verbose = FALSE,
                                        normalization.method = "SCT", 
                                        dims = 1:npcs,
                                        k.weight = 35
)



# Find neighbors and clusters, and create UMAP
pbmc_merge_lym_pass4 <- FindNeighbors(pbmc_merge_lym_pass4, dims = 1:npcs, verbose = FALSE, reduction = "integrated.rpca")
pbmc_merge_lym_pass4 <- RunUMAP(pbmc_merge_lym_pass4, dims = 1:npcs, verbose = FALSE, 
                                reduction = "integrated.rpca", reduction.name = "umap.integrated.rpca")
pbmc_merge_lym_pass4 <- FindClusters(pbmc_merge_lym_pass4, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), verbose = FALSE)



# Label 
pbmc_merge_lym_pass4@meta.data$ct_clusters <- NA
pbmc_merge_lym_pass4@meta.data$ct_clusters[
  pbmc_merge_lym_pass4@meta.data$SCT_snn_res.0.2 %in% c("0")] <- "T cells"
pbmc_merge_lym_pass4@meta.data$ct_clusters[
  pbmc_merge_lym_pass4@meta.data$SCT_snn_res.0.2 %in% c("1")] <- "NK/NKT"
pbmc_merge_lym_pass4@meta.data$ct_clusters[
  pbmc_merge_lym_pass4@meta.data$SCT_snn_res.0.2 %in% c("2")] <- "T cells"
pbmc_merge_lym_pass4@meta.data$ct_clusters[
  pbmc_merge_lym_pass4@meta.data$SCT_snn_res.0.2 %in% c("3")] <- "T cells"
pbmc_merge_lym_pass4@meta.data$ct_clusters[
  pbmc_merge_lym_pass4@meta.data$SCT_snn_res.0.2 %in% c("4")] <- "B"
pbmc_merge_lym_pass4@meta.data$ct_clusters[
  pbmc_merge_lym_pass4@meta.data$SCT_snn_res.0.2 %in% c("5")] <- "T cells"
pbmc_merge_lym_pass4@meta.data$ct_clusters[
  pbmc_merge_lym_pass4@meta.data$SCT_snn_res.0.2 %in% c("6")] <- "pDC"


#==============================================================================#
# Integrate T cells, pass 1 ----
#==============================================================================#
t_cells <- c("0", "2", "3", "5")
# we are using SCT_snn_res.0.2 
Idents(pbmc_merge_lym_pass4) <- "SCT_snn_res.0.2"
# get number of cells in each cluster
table(Idents(pbmc_merge_lym_pass4))
#0     1     2     3     4     5     6 
#12439  9879  8579  7159  2354  1111   250

# Get cell names from the clusters that I think are lymphoid cells
pbmc_merge_T_cells <- subset(pbmc_merge_lym_pass4, idents = t_cells)
table(pbmc_merge_T_cells@meta.data$SCT_snn_res.0.2)
#    0     1     2     3     4     5     6 
# 12439     0  8579  7159     0  1111     0 


# Drop levels
pbmc_merge_T_cells@meta.data[["SCT_snn_res.0.2"]] <- as.factor(pbmc_merge_T_cells@meta.data[["SCT_snn_res.0.2"]])
pbmc_merge_T_cells@meta.data[["SCT_snn_res.0.2"]] <- droplevels(pbmc_merge_T_cells@meta.data[["SCT_snn_res.0.2"]])

DefaultAssay(pbmc_merge_T_cells) <- "SoupX_RNA"
pbmc_merge_T_cells <- JoinLayers(pbmc_merge_T_cells)
#pbmc_merge_T_cells

# Split layers by RUN
#pbmc_merge_T_cells[["SoupX_RNA"]] <- split(pbmc_merge_T_cells[["SoupX_RNA"]], f = pbmc_merge_T_cells$UPN_Cycle_Day_Sample_Type_Batch)
pbmc_merge_T_cells[["SoupX_RNA"]] <- split(pbmc_merge_T_cells[["SoupX_RNA"]], f = pbmc_merge_T_cells$Batch)
#pbmc_merge_T_cells


# SCTransform is performed per sample, since they are split into layers
pbmc_merge_T_cells <- SCTransform(pbmc_merge_T_cells, assay = "SoupX_RNA", 
                                  vst.flavor = "v2",
                                  vars.to.regress = c("S.Score", "G2M.Score"),
                                  return.only.var.genes = FALSE,
                                  variable.features.n = 2500)


# Get variable features minus TRBV genes
orig_var_feat_tcells <- pbmc_merge_T_cells@assays[["SCT"]]@var.features
filt_var_feat_tcells <- grep("^TRBV", orig_var_feat_tcells, invert = TRUE, value = TRUE)


pbmc_merge_T_cells <- RunPCA(pbmc_merge_T_cells, features = filt_var_feat_tcells)
npcs <- min(get_pcs(pbmc_merge_T_cells))
npcs # 14 (2500)
Sys.time()


min(table(pbmc_merge_T_cells$UPN_Cycle_Day_Sample_Type_Batch))
min(table(pbmc_merge_T_cells$Batch))
#[1] 320

# Integrate by Batch when getting this error: Error in idx[i, ] <- res[[i]][[1]] : 
#number of items to replace is not a multiple of replacement length
# Integrate the layers using the SCT values
pbmc_merge_T_cells <- IntegrateLayers(object = pbmc_merge_T_cells, 
                                      method = RPCAIntegration,
                                      assay = "SCT", # either specify here or run default assay to SCT
                                      orig.reduction = "pca", 
                                      new.reduction = "integrated.rpca",
                                      verbose = FALSE,
                                      normalization.method = "SCT", 
                                      dims = 1:npcs#,
                                      #k.weight = 20
)



# Find neighbors and clusters, and create UMAP
pbmc_merge_T_cells <- FindNeighbors(pbmc_merge_T_cells, dims = 1:npcs, verbose = FALSE, reduction = "integrated.rpca")
pbmc_merge_T_cells <- RunUMAP(pbmc_merge_T_cells, dims = 1:npcs, verbose = FALSE, 
                              reduction = "integrated.rpca", reduction.name = "umap.integrated.rpca")
pbmc_merge_T_cells <- FindClusters(pbmc_merge_T_cells, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), verbose = FALSE)


#------------------------------------------------------------------------------#
## Sub-cluster ----
#------------------------------------------------------------------------------#
# 5 is doulet and part of cluster 3
Idents(pbmc_merge_T_cells) <- "SCT_snn_res.0.2"
pbmc_merge_T_cells_sbclt <- FindSubCluster(pbmc_merge_T_cells, "3", 
                                           resolution = 0.12, graph.name = "SCT_snn")
Idents(pbmc_merge_T_cells_sbclt) <- "sub.cluster"
pbmc_merge_T_cells_sbclt <- FindSubCluster(pbmc_merge_T_cells_sbclt, "3_1", 
                                           resolution = 0.15, graph.name = "SCT_snn")


#==============================================================================#
# Integrate T cells, pass 2 ----
#==============================================================================#
t_cells_pass2 <- c("0", "1", "2_0", "2_1_0", "2_1_1", "3_0", "3_1_0", "4")
# we are using sub.cluster 
Idents(pbmc_merge_T_cells_sbclt) <- "sub.cluster"
# get number of cells in each cluster
table(Idents(pbmc_merge_T_cells_sbclt))
#3_0     1     4     0   2_0     5 3_1_0 3_1_1 2_1_0 2_1_1 
#2782  7425  1886  8692  6050   572   943   165   523   250 

# Get cell names from the clusters that I think are lymphoid cells
pbmc_merge_T_cells_pass2 <- subset(pbmc_merge_T_cells_sbclt, idents = t_cells_pass2)
table(pbmc_merge_T_cells_pass2@meta.data$sub.cluster)
#0     1   2_0 2_1_0 2_1_1   3_0 3_1_0     4 
#8692  7425  6050   523   250  2782   943  1886 


# Drop levels
pbmc_merge_T_cells_pass2@meta.data[["sub.cluster"]] <- as.factor(pbmc_merge_T_cells_pass2@meta.data[["sub.cluster"]])
pbmc_merge_T_cells_pass2@meta.data[["sub.cluster"]] <- droplevels(pbmc_merge_T_cells_pass2@meta.data[["sub.cluster"]])

DefaultAssay(pbmc_merge_T_cells_pass2) <- "SoupX_RNA"
pbmc_merge_T_cells_pass2 <- JoinLayers(pbmc_merge_T_cells_pass2)
#pbmc_merge_T_cells_pass2

# Split layers by RUN
#pbmc_merge_T_cells_pass2[["SoupX_RNA"]] <- split(pbmc_merge_T_cells_pass2[["SoupX_RNA"]], f = pbmc_merge_T_cells_pass2$UPN_Cycle_Day_Sample_Type_Batch)
pbmc_merge_T_cells_pass2[["SoupX_RNA"]] <- split(pbmc_merge_T_cells_pass2[["SoupX_RNA"]], f = pbmc_merge_T_cells_pass2$Batch)


# SCTransform is performed per sample, since they are split into layers
pbmc_merge_T_cells_pass2 <- SCTransform(pbmc_merge_T_cells_pass2, assay = "SoupX_RNA", 
                                        vst.flavor = "v2",
                                        vars.to.regress = c("S.Score", "G2M.Score"),
                                        return.only.var.genes = FALSE,
                                        variable.features.n = 2500)


# Get variable features minus TRBV genes
orig_var_feat_tcells <- pbmc_merge_T_cells_pass2@assays[["SCT"]]@var.features
filt_var_feat_tcells <- grep("^TRBV", orig_var_feat_tcells, invert = TRUE, value = TRUE)
length(orig_var_feat_tcells)
length(filt_var_feat_tcells)

pbmc_merge_T_cells_pass2 <- RunPCA(pbmc_merge_T_cells_pass2, features = filt_var_feat_tcells)
npcs <- min(get_pcs(pbmc_merge_T_cells_pass2))
npcs # 14 (2500)
Sys.time()

min(table(pbmc_merge_T_cells_pass2$UPN_Cycle_Day_Sample_Type_Batch))
min(table(pbmc_merge_T_cells_pass2$Batch))
#[1] 307

# Integrate by Batch when getting this error: Error in idx[i, ] <- res[[i]][[1]] : 
#number of items to replace is not a multiple of replacement length
# Integrate the layers using the SCT values
pbmc_merge_T_cells_pass2 <- IntegrateLayers(object = pbmc_merge_T_cells_pass2, 
                                            method = RPCAIntegration,
                                            assay = "SCT", # either specify here or run default assay to SCT
                                            orig.reduction = "pca", 
                                            new.reduction = "integrated.rpca",
                                            verbose = FALSE,
                                            normalization.method = "SCT", 
                                            dims = 1:npcs#,
                                            #k.weight = 20
)



# Find neighbors and clusters, and create UMAP
pbmc_merge_T_cells_pass2 <- FindNeighbors(pbmc_merge_T_cells_pass2, dims = 1:npcs, verbose = FALSE, reduction = "integrated.rpca")
pbmc_merge_T_cells_pass2 <- RunUMAP(pbmc_merge_T_cells_pass2, dims = 1:npcs, verbose = FALSE, 
                                    reduction = "integrated.rpca", reduction.name = "umap.integrated.rpca")
pbmc_merge_T_cells_pass2 <- FindClusters(pbmc_merge_T_cells_pass2, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), verbose = FALSE)



DefaultAssay(pbmc_merge_T_cells_pass2)
# SoupX_RNA
DefaultAssay(pbmc_merge_T_cells_pass2) <- "RNA"
pbmc_merge_T_cells_pass2 <- JoinLayers(pbmc_merge_T_cells_pass2)
pbmc_merge_T_cells_pass2@meta.data$IL13OPCounts <- pbmc_merge_T_cells_pass2[["RNA"]]$counts["IL13OP",]
pbmc_merge_T_cells_pass2@meta.data$CART <- ifelse((pbmc_merge_T_cells_pass2@meta.data$IL13OPCounts >= 3), "Positive", "Negative")
table(pbmc_merge_T_cells_pass2@meta.data$CART)
DefaultAssay(pbmc_merge_T_cells_pass2) <- "SoupX_RNA"


#------------------------------------------------------------------------------#
## Sub-cluster ----
#------------------------------------------------------------------------------#
# Cluster 2 
Idents(pbmc_merge_T_cells_pass2) <- "SCT_snn_res.0.2"
pbmc_merge_T_cells_pass2_sub <- FindSubCluster(pbmc_merge_T_cells_pass2, 
                                                 "2", resolution = 0.3, graph.name = "SCT_snn")

pbmc_merge_T_cells_pass2_sub@meta.data$ct_clusters <- NA

pbmc_merge_T_cells_pass2_sub@meta.data$ct_clusters[
  pbmc_merge_T_cells_pass2_sub@meta.data$sub.cluster %in% c("0")] <- "T1"
pbmc_merge_T_cells_pass2_sub@meta.data$ct_clusters[
  pbmc_merge_T_cells_pass2_sub@meta.data$sub.cluster %in% c("1")] <- "T2"
pbmc_merge_T_cells_pass2_sub@meta.data$ct_clusters[
  pbmc_merge_T_cells_pass2_sub@meta.data$sub.cluster %in% c("2_0")] <- "T7"
pbmc_merge_T_cells_pass2_sub@meta.data$ct_clusters[
  pbmc_merge_T_cells_pass2_sub@meta.data$sub.cluster %in% c("2_1")] <- "T8"
pbmc_merge_T_cells_pass2_sub@meta.data$ct_clusters[
  pbmc_merge_T_cells_pass2_sub@meta.data$sub.cluster %in% c("2_2")] <- "T9"

pbmc_merge_T_cells_pass2_sub@meta.data$ct_clusters[
  pbmc_merge_T_cells_pass2_sub@meta.data$sub.cluster %in% c("2_3")] <- "T3"
pbmc_merge_T_cells_pass2_sub@meta.data$ct_clusters[
  pbmc_merge_T_cells_pass2_sub@meta.data$sub.cluster %in% c("2_4")] <- "T10"
pbmc_merge_T_cells_pass2_sub@meta.data$ct_clusters[
  pbmc_merge_T_cells_pass2_sub@meta.data$sub.cluster %in% c("2_5")] <- "T4"
pbmc_merge_T_cells_pass2_sub@meta.data$ct_clusters[
  pbmc_merge_T_cells_pass2_sub@meta.data$sub.cluster %in% c("3")] <- "T5"
pbmc_merge_T_cells_pass2_sub@meta.data$ct_clusters[
  pbmc_merge_T_cells_pass2_sub@meta.data$sub.cluster %in% c("4")] <- "T6"


# Then subcluster T3
Idents(pbmc_merge_T_cells_pass2_sub) <- "ct_clusters_num"
pbmc_merge_T_cells_pass2_sub <- FindSubCluster(pbmc_merge_T_cells_pass2_sub, "T3", resolution = 0.2, graph.name = "SCT_snn")

# Then subcluster T1
Idents(pbmc_merge_T_cells_pass2_sub) <- "sub.cluster"
pbmc_merge_T_cells_pass2_sub <- FindSubCluster(pbmc_merge_T_cells_pass2_sub, "T1", resolution = 0.2, graph.name = "SCT_snn")


Idents(pbmc_merge_T_cells_pass2_sub) <- "sub.cluster"
pbmc_merge_T_cells_pass2_sub <- FindSubCluster(pbmc_merge_T_cells_pass2_sub, "T1_1", resolution = 0.15, graph.name = "SCT_snn")


#------------------------------------------------------------------------------#
## Label T cell states ----
#------------------------------------------------------------------------------#
pbmc_merge_T_cells_pass2_sub@meta.data$ct_final <- NA # or any other initialization value
pbmc_merge_T_cells_pass2_sub@meta.data$ct_final[pbmc_merge_T_cells_pass2_sub@meta.data$sub.cluster == "T1_0"] <- "Effector"
pbmc_merge_T_cells_pass2_sub@meta.data$ct_final[pbmc_merge_T_cells_pass2_sub@meta.data$sub.cluster == "T1_1_0"] <- "Effector-memory"
pbmc_merge_T_cells_pass2_sub@meta.data$ct_final[pbmc_merge_T_cells_pass2_sub@meta.data$sub.cluster == "T1_1_1"] <- "Effector"
pbmc_merge_T_cells_pass2_sub@meta.data$ct_final[pbmc_merge_T_cells_pass2_sub@meta.data$sub.cluster == "T1_2"] <- "Effector"
pbmc_merge_T_cells_pass2_sub@meta.data$ct_final[pbmc_merge_T_cells_pass2_sub@meta.data$sub.cluster == "T2"] <- "Naïve"
pbmc_merge_T_cells_pass2_sub@meta.data$ct_final[pbmc_merge_T_cells_pass2_sub@meta.data$sub.cluster == "T3_0"] <- "Treg"
pbmc_merge_T_cells_pass2_sub@meta.data$ct_final[pbmc_merge_T_cells_pass2_sub@meta.data$sub.cluster == "T3_1"] <- "GAPDH+ Glycolysis"
pbmc_merge_T_cells_pass2_sub@meta.data$ct_final[pbmc_merge_T_cells_pass2_sub@meta.data$sub.cluster == "T3_2"] <- "Proliferating"
pbmc_merge_T_cells_pass2_sub@meta.data$ct_final[pbmc_merge_T_cells_pass2_sub@meta.data$sub.cluster == "T4"] <- "Central memory"
pbmc_merge_T_cells_pass2_sub@meta.data$ct_final[pbmc_merge_T_cells_pass2_sub@meta.data$sub.cluster == "T5"] <- "Effector"
pbmc_merge_T_cells_pass2_sub@meta.data$ct_final[pbmc_merge_T_cells_pass2_sub@meta.data$sub.cluster == "T6"] <- "MAIT"
pbmc_merge_T_cells_pass2_sub@meta.data$ct_final[pbmc_merge_T_cells_pass2_sub@meta.data$sub.cluster == "T7"] <- "Central memory"
pbmc_merge_T_cells_pass2_sub@meta.data$ct_final[pbmc_merge_T_cells_pass2_sub@meta.data$sub.cluster == "T8"] <- "Memory"
pbmc_merge_T_cells_pass2_sub@meta.data$ct_final[pbmc_merge_T_cells_pass2_sub@meta.data$sub.cluster == "T9"] <- "Memory"
pbmc_merge_T_cells_pass2_sub@meta.data$ct_final[pbmc_merge_T_cells_pass2_sub@meta.data$sub.cluster == "T10"] <- "Treg"


#------------------------------------------------------------------------------#
## Add cell types and T cell states to the lymphoid object ----
#------------------------------------------------------------------------------#
pbmc_merge_lym_pass4_meta <- pbmc_merge_lym_pass4@meta.data
pbmc_merge_lym_pass4_meta$cellID <- rownames(pbmc_merge_lym_pass4_meta)
pbmc_merge_T_cells_pass2_sub_meta <- pbmc_merge_T_cells_pass2_sub@meta.data %>% dplyr::select(ct_final)
pbmc_merge_T_cells_pass2_sub_meta$cellID <- rownames(pbmc_merge_T_cells_pass2_sub_meta)
colnames(pbmc_merge_T_cells_pass2_sub_meta) <- c("ct_final_T", "cellID")  


pbmc_merge_lym_pass4_meta_join <- left_join(pbmc_merge_lym_pass4_meta, pbmc_merge_T_cells_pass2_sub_meta)
pbmc_merge_lym_pass4@meta.data <- pbmc_merge_lym_pass4_meta_join
rownames(pbmc_merge_lym_pass4@meta.data) <- pbmc_merge_lym_pass4_meta$cellID


pbmc_merge_lym_pass4@meta.data$ct_final_all <- pbmc_merge_lym_pass4@meta.data$ct_final # or any other initialization value
pbmc_merge_lym_pass4@meta.data$ct_final_all[pbmc_merge_lym_pass4@meta.data$ct_final_T == "Effector"] <- "T (Effector)"
pbmc_merge_lym_pass4@meta.data$ct_final_all[pbmc_merge_lym_pass4@meta.data$ct_final_T == "Effector-memory"] <- "T (Effector-memory)"
pbmc_merge_lym_pass4@meta.data$ct_final_all[pbmc_merge_lym_pass4@meta.data$ct_final_T == "Naïve"] <- "T (Naïve)"
pbmc_merge_lym_pass4@meta.data$ct_final_all[pbmc_merge_lym_pass4@meta.data$ct_final_T == "Treg"] <- "Treg"
pbmc_merge_lym_pass4@meta.data$ct_final_all[pbmc_merge_lym_pass4@meta.data$ct_final_T == "GAPDH+ Glycolysis"] <- "T (GAPDH+ Glycolysis)"
pbmc_merge_lym_pass4@meta.data$ct_final_all[pbmc_merge_lym_pass4@meta.data$ct_final_T == "Proliferating"] <- "T (Proliferating)"
pbmc_merge_lym_pass4@meta.data$ct_final_all[pbmc_merge_lym_pass4@meta.data$ct_final_T == "Central memory"] <- "T (Central memory)"
pbmc_merge_lym_pass4@meta.data$ct_final_all[pbmc_merge_lym_pass4@meta.data$ct_final_T == "MAIT"] <- "T (MAIT)"
pbmc_merge_lym_pass4@meta.data$ct_final_all[pbmc_merge_lym_pass4@meta.data$ct_final_T == "Memory"] <- "T (Memory)"


#------------------------------------------------------------------------------#
## Add TCR data to lymphoid object ----
#------------------------------------------------------------------------------#
combined2 <- readRDS("/scratch/aoill/projects/CAR-T/00_new_2025/geo_2026/pbmc_scRepertoire_output.rds")

# Pull out T cells
table(pbmc_merge_lym_pass4@meta.data$ct_final)
t_cells <- c("T (Effector)", "T (Naïve)", "T (MAIT)",
             "T (Effector-memory)", "T (Memory)", "T (GAPDH+ Glycolysis)",
             "T (Central memory)", "Treg", "T (Proliferating)"
)
pbmc_tcells_TCR <- subset(pbmc_merge_lym_pass4, subset = ct_final_all %in% t_cells)

pbmc_tcells_TCR <- combineExpression(combined2, pbmc_tcells,
                               group.by = "UPN_Cycle_Day_Sample_Type",
                               filterNA = FALSE,
                               proportion = FALSE,
                               cloneCall="strict",
                               cloneSize=c(Single=1, Small=5, Medium=20,
                                           Large=100, Hyperexpanded=500)
)
pbmc_tcells_TCR@meta.data$cloneSize[is.na(pbmc_tcells_TCR_TCR@meta.data$cloneSize)] <- "None ( < X <= 0)"



# calculate Frequency_norm
tmp_meta <- pbmc_tcells_TCR@meta.data

# Add a T cell number count for each UPN patient and cycle
tmp_meta <- pbmc_tcells_TCR@meta.data %>%
  #add_count(UPN, name = "UPN_n")
  add_count(UPN_Cycle, name = "UPN_Cycle_n")

# Get normalized TCR frequency
tmp_meta$Frequency_norm <- tmp_meta$clonalFrequency/tmp_meta$UPN_Cycle_n

pbmc_tcells_TCR@meta.data$UPN_Cycle_n <- tmp_meta$UPN_Cycle_n
pbmc_tcells_TCR@meta.data$Frequency_norm <- tmp_meta$Frequency_norm


tmp_meta_sub <- tmp_meta %>%
  dplyr::select(UPN, Sample_Type, Cycle, CART, CTstrict, clonalFrequency,
                Frequency_norm)
rownames(tmp_meta) <- NULL
tmp_meta_sub <- na.omit(tmp_meta_sub)

# There could be duplicated rows in this metadata (clonotypes with freq >1 will 
# come up multiple times)
tmp_meta_sub_unique <- unique(tmp_meta_sub)

tmp_meta <- tmp_meta_sub_unique

# Remove any TCR with a frequency of 1
tmp_meta <- tmp_meta %>% dplyr::filter(clonalFrequency > 1)

# Get cut offs
n <- 20
TCRs_top20 <- tmp_meta[tmp_meta$Frequency_norm > quantile(tmp_meta$Frequency_norm,prob=1-n/100),]
min(TCRs_top20$Frequency_norm)


n <- 5
TCRs_top5 <- tmp_meta[tmp_meta$Frequency_norm > quantile(tmp_meta$Frequency_norm,prob=1-n/100),]
min(TCRs_top5$Frequency_norm)


n <- 1
TCRs_top1 <- tmp_meta[tmp_meta$Frequency_norm > quantile(tmp_meta$Frequency_norm,prob=1-n/100),]
min(TCRs_top1$Frequency_norm)


ggplot(tmp_meta, aes(x=Frequency_norm)) + 
  geom_histogram(color="black", fill="white", bins = 50#, 
                 #aes(y=..density..), position="identity"
  ) +
  #geom_density() +
  scale_x_log10() +
  xlab("Normalized TCR Freq") + ylab("Count") +
  theme_classic() +
  geom_vline(xintercept = min(TCRs_top20$Frequency_norm), linetype="dotted", 
             color = "red") +
  geom_vline(xintercept = min(TCRs_top5$Frequency_norm), linetype="dotted", 
             color = "red") +
  geom_vline(xintercept = min(TCRs_top1$Frequency_norm), linetype="dotted", 
             color = "red")


# Get expanded vs not expanded categories with expanded categories separated by
# how expanded they are
pbmc_tcells_TCR@meta.data <- pbmc_tcells_TCR@meta.data %>% 
  mutate(Frequency_norm_cat = case_when(clonalFrequency == 1 ~ "Not Expanded",
                                        Frequency_norm < min(TCRs_top20$Frequency_norm) ~ "Not Expanded",
                                        (Frequency_norm >= min(TCRs_top20$Frequency_norm) & Frequency_norm < min(TCRs_top5$Frequency_norm)) ~ "Expanded (Top 20% - 5%)",
                                        (Frequency_norm >= min(TCRs_top5$Frequency_norm) & Frequency_norm < min(TCRs_top1$Frequency_norm)) ~ "More Expanded (Top 5% - 1%)",
                                        Frequency_norm >= min(TCRs_top1$Frequency_norm) ~ "Most Expanded (Top 1%)",
                                        
                                        TRUE ~ as.character(NA))) 
table(pbmc_tcells_TCR@meta.data$Frequency_norm_cat)


# Get expanded vs not expanded categories
pbmc_tcells_TCR@meta.data <- pbmc_tcells_TCR@meta.data %>% 
  mutate(Frequency_norm_cat_2 = case_when(clonalFrequency == 1 ~ "Not Expanded",
                                          Frequency_norm < min(TCRs_top20$Frequency_norm) ~ "Not Expanded",
                                          (Frequency_norm >= min(TCRs_top20$Frequency_norm) & Frequency_norm < min(TCRs_top5$Frequency_norm)) ~ "Expanded",
                                          (Frequency_norm >= min(TCRs_top5$Frequency_norm) & Frequency_norm < min(TCRs_top1$Frequency_norm)) ~ "Expanded",
                                          Frequency_norm >= min(TCRs_top1$Frequency_norm) ~ "Expanded",
                                          
                                          TRUE ~ as.character(NA))) 
table(pbmc_tcells_TCR@meta.data$Frequency_norm_cat_2)



# add in T cell TCR frequency back into the lymphoid object (will also need to
# add it to the T cell object TODO)
#colnames(pbmc_tcells_TCR@meta.data)
pbmc_tcells_TCR_tcr_meta <- pbmc_tcells_TCR@meta.data %>%
  dplyr::select(cellID, 
                CTgene, CTnt, CTaa, CTstrict, clonalProportion, clonalFrequency, 
                cloneSize, 
                Frequency_norm, Frequency_norm_cat, Frequency_norm_cat_2)
pbmc_tcells_TCR_meta_rnames <- rownames(pbmc_tcells_TCR@meta.data)

pbmc_merge_lym_pass4_meta <- pbmc_merge_lym_pass4@meta.data
pbmc_merge_lym_pass4_meta_rnames <- rownames(pbmc_merge_lym_pass4@meta.data)


dim(pbmc_merge_lym_pass4_meta)
dim(pbmc_tcells_TCR_tcr_meta)
pbmc_merge_lym_pass4_meta_l_join <- left_join(pbmc_merge_lym_pass4_meta, pbmc_tcells_TCR_tcr_meta)
dim(pbmc_merge_lym_pass4_meta_l_join) 


pbmc_merge_lym_pass4@meta.data <- pbmc_merge_lym_pass4_meta_l_join
rownames(pbmc_merge_lym_pass4@meta.data) <- pbmc_merge_lym_pass4_meta_rnames


# Reorder sample
pbmc_merge_lym_pass4@meta.data$UPN_Cycle <- factor(x = pbmc_merge_lym_pass4@meta.data$UPN_Cycle, 
                                       levels = c("514_2", "514_5", "514_8", "514_11", "514_13",
                                                  "515_4", "515_8",
                                                  "574_1", "574_3", "574_8",  
                                                  "689_1", "689_4", "689_8", "689_12", "689_17", 
                                                  "692_2", "692_4",
                                                  "705_5", "705_8",  
                                                  "716_2", "716_4", "716_8"))


pbmc_merge_lym_pass4@meta.data$Frequency_norm_cat <- factor(pbmc_merge_lym_pass4@meta.data$Frequency_norm_cat, 
                                                levels = c("Most Expanded (Top 1%)", 
                                                           "More Expanded (Top 5% - 1%)", 
                                                           "Expanded (Top 20% - 5%)",
                                                           "Not Expanded",
                                                           "NA"))


#------------------------------------------------------------------------------#
## Add metadata ----
#------------------------------------------------------------------------------#
# Response categories
pbmc_merge_lym_pass4@meta.data$Response_cat <- "NA" # or any other initialization value
pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "514_2"] <- "Baseline"
pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "514_5"] <- "Non-response"
pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "514_8"] <- "Response"
pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "514_11"] <- "Non-response"
pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "514_13"] <- "Non-response"

pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "515_4"] <- "Response"
pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "515_8"] <- "Non-response"

pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "574_1"] <- "Baseline"
pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "574_3"] <- "Non-response"
pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "574_8"] <- "Non-response"

pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "689_1"] <- "Baseline"
pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "689_4"] <- "Response"
pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "689_8"] <- "Response"
pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "689_12"] <- "Non-response"
pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "689_17"] <- "Non-response"

pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "692_2"] <- "Baseline"
pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "692_4"] <- "Response"

pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "705_1"] <- "Baseline"
pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "705_5"] <- "Non-response"
pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "705_8"] <- "Non-response"

pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "716_2"] <- "Baseline"
pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "716_4"] <- "Response"
pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "716_8"] <- "Response"

pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "625_1"] <- "Baseline" 
pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "625_3"] <- "Non-response" 
pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "625_4"] <- "Non-response" 

pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "626_1"] <- "Baseline" 
pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "626_4"] <- "Response" 
pbmc_merge_lym_pass4@meta.data$Response_cat[pbmc_merge_lym_pass4@meta.data$UPN_Cycle == "626_7"] <- "Non-response" 

table(pbmc_merge_lym_pass4@meta.data$Response_cat)


# add disease info
pbmc_merge_lym_pass4@meta.data$tumor_type <- NA # or any other initialization value
pbmc_merge_lym_pass4@meta.data$tumor_type[pbmc_merge_lym_pass4@meta.data$UPN == "514"] <- "Ependymoma"
pbmc_merge_lym_pass4@meta.data$tumor_type[pbmc_merge_lym_pass4@meta.data$UPN == "574"] <- "Ependymoma"
pbmc_merge_lym_pass4@meta.data$tumor_type[pbmc_merge_lym_pass4@meta.data$UPN == "625"] <- "Ependymoma" 
pbmc_merge_lym_pass4@meta.data$tumor_type[pbmc_merge_lym_pass4@meta.data$UPN == "689"] <- "Ependymoma"

pbmc_merge_lym_pass4@meta.data$tumor_type[pbmc_merge_lym_pass4@meta.data$UPN == "515"] <- "DMG"
pbmc_merge_lym_pass4@meta.data$tumor_type[pbmc_merge_lym_pass4@meta.data$UPN == "626"] <- "DMG" 
pbmc_merge_lym_pass4@meta.data$tumor_type[pbmc_merge_lym_pass4@meta.data$UPN == "692"] <- "DMG"
pbmc_merge_lym_pass4@meta.data$tumor_type[pbmc_merge_lym_pass4@meta.data$UPN == "716"] <- "DMG"

pbmc_merge_lym_pass4@meta.data$tumor_type[pbmc_merge_lym_pass4@meta.data$UPN == "705"] <- "Glioblastoma"

head(pbmc_merge_lym_pass4@meta.data)


# add lymphodepletion status
pbmc_merge_lym_pass4@meta.data$Lymphodepletion <- NA # or any other initialization value
pbmc_merge_lym_pass4@meta.data$Lymphodepletion[pbmc_merge_lym_pass4@meta.data$UPN == "514"] <- "Non-Lymphodepleted"
pbmc_merge_lym_pass4@meta.data$Lymphodepletion[pbmc_merge_lym_pass4@meta.data$UPN == "515"] <- "Non-Lymphodepleted"
pbmc_merge_lym_pass4@meta.data$Lymphodepletion[pbmc_merge_lym_pass4@meta.data$UPN == "574"] <- "Non-Lymphodepleted" 
pbmc_merge_lym_pass4@meta.data$Lymphodepletion[pbmc_merge_lym_pass4@meta.data$UPN == "625"] <- "Lymphodepleted" 
pbmc_merge_lym_pass4@meta.data$Lymphodepletion[pbmc_merge_lym_pass4@meta.data$UPN == "626"] <- "Lymphodepleted" 
pbmc_merge_lym_pass4@meta.data$Lymphodepletion[pbmc_merge_lym_pass4@meta.data$UPN == "689"] <- "Lymphodepleted" 
pbmc_merge_lym_pass4@meta.data$Lymphodepletion[pbmc_merge_lym_pass4@meta.data$UPN == "705"] <- "Lymphodepleted" 
pbmc_merge_lym_pass4@meta.data$Lymphodepletion[pbmc_merge_lym_pass4@meta.data$UPN == "716"] <- "Lymphodepleted" 
pbmc_merge_lym_pass4@meta.data$Lymphodepletion[pbmc_merge_lym_pass4@meta.data$UPN == "692"] <- "Lymphodepleted" 


#------------------------------------------------------------------------------#
## Save ----
#------------------------------------------------------------------------------#
saveRDS(pbmc_merge_lym_pass4, "/scratch/aoill/projects/CAR-T/00_new_2025/geo_2026/pbmc_lymphoid_seurat_obj.rds")


#==============================================================================#
# Integrate Myeloid, pass 2 ----
#==============================================================================#
# Myeloid pass 2: Remove cluster 6 from pbmc_merge_sub_mye and re-integrate 
# pbmc_merge_sub_mye, SCT_snn_res.0.2
myeloid_pass2 <- c("0", "1", "2", "3", "4", "5")

# pbmc_merge_sub_mye, SCT_snn_res.0.2
myeloid_pass2 <- c("0", "1", "2", "3", "4", "5")


Idents(pbmc_merge_sub_mye) <- "SCT_snn_res.0.2"
# get number of cells in each cluster
table(Idents(pbmc_merge_sub_mye))
#0     1     2     3     4     5     6 
#40718 30436 11231  8525  6077  2375   152 

# Get cell names from the clusters that I think are lymphoid cells
pbmc_merge_mye_pass2 <- subset(pbmc_merge_sub_mye, idents = myeloid_pass2)
table(pbmc_merge_mye_pass2@meta.data$SCT_snn_res.0.2)
#0     1     2     3     4     5     6 
#40718 30436 11231  8525  6077  2375     0 

# Drop levels
pbmc_merge_mye_pass2@meta.data[["SCT_snn_res.0.2"]] <- as.factor(pbmc_merge_mye_pass2@meta.data[["SCT_snn_res.0.2"]])
pbmc_merge_mye_pass2@meta.data[["SCT_snn_res.0.2"]] <- droplevels(pbmc_merge_mye_pass2@meta.data[["SCT_snn_res.0.2"]])

DefaultAssay(pbmc_merge_mye_pass2) <- "SoupX_RNA"
pbmc_merge_mye_pass2 <- JoinLayers(pbmc_merge_mye_pass2)
#pbmc_merge_mye_pass2

# Split layers by RUN
pbmc_merge_mye_pass2[["SoupX_RNA"]] <- split(pbmc_merge_mye_pass2[["SoupX_RNA"]], f = pbmc_merge_mye_pass2$UPN_Cycle_Day_Sample_Type_Batch)
#pbmc_merge_mye_pass2


# SCTransform is performed per sample, since they are split into layers
pbmc_merge_mye_pass2 <- SCTransform(pbmc_merge_mye_pass2, assay = "SoupX_RNA", 
                                    vst.flavor = "v2",
                                    vars.to.regress = c("S.Score", "G2M.Score"),
                                    return.only.var.genes = FALSE,
                                    variable.features.n = 2500)

pbmc_merge_mye_pass2 <- RunPCA(pbmc_merge_mye_pass2)
npcs <- min(get_pcs(pbmc_merge_mye_pass2))
npcs # 16 (2500)
Sys.time()

min(table(pbmc_merge_mye_pass2$UPN_Cycle_Day_Sample_Type_Batch))
#[1] 32

# Integrate the layers using the SCT values
pbmc_merge_mye_pass2 <- IntegrateLayers(object = pbmc_merge_mye_pass2, 
                                        method = RPCAIntegration,
                                        assay = "SCT", # either specify here or run default assay to SCT
                                        orig.reduction = "pca", 
                                        new.reduction = "integrated.rpca",
                                        verbose = FALSE,
                                        normalization.method = "SCT", 
                                        dims = 1:npcs,
                                        k.weight = 25
)



# Find neighbors and clusters, and create UMAP
pbmc_merge_mye_pass2 <- FindNeighbors(pbmc_merge_mye_pass2, dims = 1:npcs, verbose = FALSE, reduction = "integrated.rpca")
pbmc_merge_mye_pass2 <- RunUMAP(pbmc_merge_mye_pass2, dims = 1:npcs, verbose = FALSE, 
                                reduction = "integrated.rpca", reduction.name = "umap.integrated.rpca")
pbmc_merge_mye_pass2 <- FindClusters(pbmc_merge_mye_pass2, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), verbose = FALSE)


#------------------------------------------------------------------------------#
## Add cell type labels ----
#------------------------------------------------------------------------------#
pbmc_merge_mye_pass2_join@meta.data$ct_clusters <- NA

pbmc_merge_mye_pass2_join@meta.data$ct_clusters[
  pbmc_merge_mye_pass2_join@meta.data$SCT_snn_res.0.2 %in% c("0")] <- "Classical Monocyte"
pbmc_merge_mye_pass2_join@meta.data$ct_clusters[
  pbmc_merge_mye_pass2_join@meta.data$SCT_snn_res.0.2 %in% c("1")] <- "S100 high, HLA class II low Monocyte"
pbmc_merge_mye_pass2_join@meta.data$ct_clusters[
  pbmc_merge_mye_pass2_join@meta.data$SCT_snn_res.0.2 %in% c("2")] <- "Intermediate Monocyte"
pbmc_merge_mye_pass2_join@meta.data$ct_clusters[
  pbmc_merge_mye_pass2_join@meta.data$SCT_snn_res.0.2 %in% c("3")] <- "Classical Monocyte"
pbmc_merge_mye_pass2_join@meta.data$ct_clusters[
  pbmc_merge_mye_pass2_join@meta.data$SCT_snn_res.0.2 %in% c("4")] <- "Non-classical Monocyte"
pbmc_merge_mye_pass2_join@meta.data$ct_clusters[
  pbmc_merge_mye_pass2_join@meta.data$SCT_snn_res.0.2 %in% c("5")] <- "Classical Monocyte"
pbmc_merge_mye_pass2_join@meta.data$ct_clusters[
  pbmc_merge_mye_pass2_join@meta.data$SCT_snn_res.0.2 %in% c("6")] <- "cDC2"


#------------------------------------------------------------------------------#
## Add metadata ----
#------------------------------------------------------------------------------#
# Response categories
pbmc_merge_mye_pass2_join@meta.data$Response_cat <- "NA" # or any other initialization value
pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "514_2"] <- "Baseline"
pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "514_5"] <- "Non-response"
pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "514_8"] <- "Response"
pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "514_11"] <- "Non-response"
pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "514_13"] <- "Non-response"

pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "515_4"] <- "Response"
pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "515_8"] <- "Non-response"

pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "574_1"] <- "Baseline"
pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "574_3"] <- "Non-response"
pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "574_8"] <- "Non-response"

pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "689_1"] <- "Baseline"
pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "689_4"] <- "Response"
pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "689_8"] <- "Response"
pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "689_12"] <- "Non-response"
pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "689_17"] <- "Non-response"

pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "692_2"] <- "Baseline"
pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "692_4"] <- "Response"

pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "705_1"] <- "Baseline"
pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "705_5"] <- "Non-response"
pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "705_8"] <- "Non-response"

pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "716_2"] <- "Baseline"
pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "716_4"] <- "Response"
pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "716_8"] <- "Response"

pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "625_1"] <- "Baseline" 
pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "625_3"] <- "Non-response" 
pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "625_4"] <- "Non-response" 

pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "626_1"] <- "Baseline" 
pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "626_4"] <- "Response" 
pbmc_merge_mye_pass2_join@meta.data$Response_cat[pbmc_merge_mye_pass2_join@meta.data$UPN_Cycle == "626_7"] <- "Non-response" 

table(pbmc_merge_mye_pass2_join@meta.data$Response_cat)


# add disease info
pbmc_merge_mye_pass2_join@meta.data$tumor_type <- NA # or any other initialization value
pbmc_merge_mye_pass2_join@meta.data$tumor_type[pbmc_merge_mye_pass2_join@meta.data$UPN == "514"] <- "Ependymoma"
pbmc_merge_mye_pass2_join@meta.data$tumor_type[pbmc_merge_mye_pass2_join@meta.data$UPN == "574"] <- "Ependymoma"
pbmc_merge_mye_pass2_join@meta.data$tumor_type[pbmc_merge_mye_pass2_join@meta.data$UPN == "625"] <- "Ependymoma" 
pbmc_merge_mye_pass2_join@meta.data$tumor_type[pbmc_merge_mye_pass2_join@meta.data$UPN == "689"] <- "Ependymoma"

pbmc_merge_mye_pass2_join@meta.data$tumor_type[pbmc_merge_mye_pass2_join@meta.data$UPN == "515"] <- "DMG"
pbmc_merge_mye_pass2_join@meta.data$tumor_type[pbmc_merge_mye_pass2_join@meta.data$UPN == "626"] <- "DMG" 
pbmc_merge_mye_pass2_join@meta.data$tumor_type[pbmc_merge_mye_pass2_join@meta.data$UPN == "692"] <- "DMG"
pbmc_merge_mye_pass2_join@meta.data$tumor_type[pbmc_merge_mye_pass2_join@meta.data$UPN == "716"] <- "DMG"

pbmc_merge_mye_pass2_join@meta.data$tumor_type[pbmc_merge_mye_pass2_join@meta.data$UPN == "705"] <- "Glioblastoma"

head(pbmc_merge_mye_pass2_join@meta.data)


# add lymphodepletion status
pbmc_merge_mye_pass2_join@meta.data$Lymphodepletion <- NA # or any other initialization value
pbmc_merge_mye_pass2_join@meta.data$Lymphodepletion[pbmc_merge_mye_pass2_join@meta.data$UPN == "514"] <- "Non-Lymphodepleted"
pbmc_merge_mye_pass2_join@meta.data$Lymphodepletion[pbmc_merge_mye_pass2_join@meta.data$UPN == "515"] <- "Non-Lymphodepleted"
pbmc_merge_mye_pass2_join@meta.data$Lymphodepletion[pbmc_merge_mye_pass2_join@meta.data$UPN == "574"] <- "Non-Lymphodepleted" 
pbmc_merge_mye_pass2_join@meta.data$Lymphodepletion[pbmc_merge_mye_pass2_join@meta.data$UPN == "625"] <- "Lymphodepleted" 
pbmc_merge_mye_pass2_join@meta.data$Lymphodepletion[pbmc_merge_mye_pass2_join@meta.data$UPN == "626"] <- "Lymphodepleted" 
pbmc_merge_mye_pass2_join@meta.data$Lymphodepletion[pbmc_merge_mye_pass2_join@meta.data$UPN == "689"] <- "Lymphodepleted" 
pbmc_merge_mye_pass2_join@meta.data$Lymphodepletion[pbmc_merge_mye_pass2_join@meta.data$UPN == "705"] <- "Lymphodepleted" 
pbmc_merge_mye_pass2_join@meta.data$Lymphodepletion[pbmc_merge_mye_pass2_join@meta.data$UPN == "716"] <- "Lymphodepleted" 
pbmc_merge_mye_pass2_join@meta.data$Lymphodepletion[pbmc_merge_mye_pass2_join@meta.data$UPN == "692"] <- "Lymphodepleted" 


#------------------------------------------------------------------------------#
## Save ----
#------------------------------------------------------------------------------#
saveRDS(pbmc_merge_mye_pass2_join, "/scratch/aoill/projects/CAR-T/00_new_2025/geo_2026/pbmc_myeloid_seurat_obj.rds")

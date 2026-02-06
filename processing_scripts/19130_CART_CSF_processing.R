#==============================================================================#
# Author : Angela M. Oill, aoill@tgen.org
# Date: 2025/04/14
# Project: Pediatric CAR-T 
# Description: Preprocessing and integration of CSF samples
#==============================================================================#


#==============================================================================#
# Load Libraries, helper modules and functions ----
#==============================================================================#
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(DropletUtils)
library(scRepertoire)
library(scGSVA)
library(RColorBrewer)
library(scRepertoire)
library(stringr)


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
sample_metadata_not_multiplexed_csf <- sample_metadata %>% 
  filter(Multiplexed == "No") %>%
  filter(Sample_Type == "CSF")

sample_metadata_multiplexed_csf <- sample_metadata %>% 
  filter(Multiplexed == "Yes") %>%
  filter(Sample_Type == "CSF")


#==============================================================================#
# Read in Seurat data as a list ----
#==============================================================================#
#--------------------#
## Original data ----
#--------------------#
# This data wasnt multiplexed. This is a list of seurat objects
seurat_list <- prep_seurat_list(
  sample_metadata_not_multiplexed_csf, 
  batch_ID = "Cell_Prefix", 
  cellRanger_path = "CellRanger_path", 
  cell_ID_prefix = "Cell_Prefix", 
  run_soupX = FALSE) # Run SoupX after filtering


#--------------------#
## New data ----
#--------------------#
# This data is multiplexed and was demultiplexed using demuxlet. I have a file 
# for each batch with which cells are singlets/doublets/ambiguous and which 
# sample it belongs to (demuxlet.best). I will need to join this to the metadata
# based on the cell ID then only keep singlets and split the objects by sample.
# After this proceed like 
metadata = sample_metadata_multiplexed_csf
batch_ID = "Cell_Prefix"
cellRanger_path = "CellRanger_path"
cell_ID_prefix = "Cell_Prefix"

sample_list <- unique(metadata[[batch_ID]])

#i = "Batch11"
sample_seurat_list <- lapply(sample_list, function(i){
  message(i)
  print(paste("Reading in ", 
              length(unique(paste0(metadata[which(metadata[[batch_ID]]==i),][[cellRanger_path]], "/outs/filtered_feature_bc_matrix/"))), " file.",
              sep = ""
  )
  )
  print(unique(paste0(metadata[which(metadata[[batch_ID]]==i),][[cellRanger_path]], "/outs/filtered_feature_bc_matrix/")))
  sample_10x_data <- Read10X(unique(paste0(metadata[which(metadata[[batch_ID]]==i),][[cellRanger_path]], "/outs/filtered_feature_bc_matrix/")))
  seurat_object <- CreateSeuratObject(counts = sample_10x_data)
  seurat_object <- PercentageFeatureSet(seurat_object, pattern = "^MT-", col.name = "percent.mt_RNA", assay = "RNA")
  seurat_object <- PercentageFeatureSet(seurat_object, pattern = "^RP[SL]|^MRP[SL]", col.name = "percent.ribo_RNA", assay = "RNA")
  seurat_object@meta.data$orig_cellID <- rownames(seurat_object@meta.data)
  
  # Add demuxlet info here (maybe add?)
  #dmx_fn <- paste("/scratch/aoill/projects/CAR-T/00_new_2025/demuxlet_outs/", tolower(i), "/demuxlet.best", sep = "")
  #dmx_results <- read.table(dmx_fn, header = T)
  
  # add souporcell results here
  soc_fn <- paste("/scratch/aoill/projects/CAR-T/00_new_2025/souporcell/", tolower(i), "_results.tsv", sep = "")
  soc_results <- read.table(soc_fn, header = T)

  
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
seurat_list_all <- c(seurat_list, sample_seurat_list)


#==============================================================================#
# Filter ----
#==============================================================================#

#----------------------------------------#
## Merge Seurat and visualize quality ----
#----------------------------------------#
seurat_merge <- merge(x = seurat_list_all[[1]], y = seurat_list_all[2:length(seurat_list_all)])
seurat_merge@meta.data$New_Ident <- "SeuratProject"
Idents(seurat_merge) <- "New_Ident"

# CSF %mt log(nFeature_RNA)
smoothScatter(seurat_merge@meta.data$percent.mt, log(seurat_merge@meta.data$nFeature_RNA),
              xlab = "% MT", ylab = "log(nFeature_RNA)",
              main = "CSF")
abline(h = log(650), v = 10)
text(15,log(750), "nFeature_RNA = 650,\npercent.mt = 10", adj = c(0, -.1))


# CSF %mt log(nCount_RNA)
smoothScatter(seurat_merge@meta.data$percent.mt, log(seurat_merge@meta.data$nCount_RNA),
              xlab = "% MT", ylab = "log(nCount_RNA)",
              main = "CSF")
abline(h = log(1200), v = 10)
text(15,log(1300), "nCount_RNA = 1200,\npercent.mt = 10", adj = c(0, -.1))

# I think the same filter as before works

#-------------#
## Filter ----
#-------------#
seurat_list_filtered <- filter_manual(seurat_list_all, 
                                      pt_mt = 10, 
                                      nFeature = 650,
                                      nCount = 1200)

seurat_list_separated <- list()
for (i in 1:length(seurat_list_filtered)) {
  print(i)
  
  pedi_batch <- seurat_list_filtered[[i]]
  pedi_batch@meta.data$Sample_ID <- paste(pedi_batch@meta.data$UPN, 
                                          pedi_batch@meta.data$Sample_Type,
                                          pedi_batch@meta.data$Cycle, 
                                          pedi_batch@meta.data$Day, 
                                          pedi_batch@meta.data$Batch,
                                          sep = "_")
  pedi_batch <- Seurat::SplitObject(pedi_batch, split.by = "Sample_ID")
  seurat_list_separated <- c(seurat_list_separated, pedi_batch)
}


for (i in 1:length(seurat_list_separated)) {
  print(paste0("Converting sample ", names(seurat_list_separated[i])))
  obj.sub <- seurat_list_separated[[i]]
  DropletUtils::write10xCounts(path = paste0("/scratch/aoill/projects/CAR-T/00_new_2025/soupX/demultiplexed_",names(seurat_list_separated[i])), 
                               x = obj.sub[["RNA"]]@layers$counts, 
                               barcodes = colnames(obj.sub[["RNA"]]), # cell names
                               gene.id = rownames(obj.sub[["RNA"]]), # Gene identifiers, one per row of X
                               type = "sparse", version="3", overwrite = T)
}


#-----------------------------------------------------#
## Get number of cells before and after filtering ----
#-----------------------------------------------------#
#seurat_merge <- merge(x = seurat_list[[1]], y = seurat_list[2:length(seurat_list)])
seurat_merge <- merge(x = seurat_list_all[[1]], y = seurat_list_all[2:length(seurat_list_all)])
seurat_merge@meta.data$UPN_Cycle_Day_Sample_Type_Batch <- paste0(seurat_merge@meta.data$UPN, "_",
                                                                 seurat_merge@meta.data$Cycle,  "_", 
                                                                 seurat_merge@meta.data$Day, "_",
                                                                 seurat_merge@meta.data$Sample_Type, "_",
                                                                 seurat_merge@meta.data$Batch)
cell_numbers <- as.data.frame(table(seurat_merge@meta.data$UPN_Cycle_Day_Sample_Type_Batch))
colnames(cell_numbers) <- c("Sample", "Unfiltered")

seurat_merge <- merge(x = seurat_list_filtered[[1]], y = seurat_list_filtered[2:length(seurat_list_filtered)])
seurat_merge@meta.data$UPN_Cycle_Day_Sample_Type_Batch <- paste0(seurat_merge@meta.data$UPN, "_",
                                                                 seurat_merge@meta.data$Cycle,  "_", 
                                                                 seurat_merge@meta.data$Day, "_",
                                                                 seurat_merge@meta.data$Sample_Type, "_",
                                                                 seurat_merge@meta.data$Batch)
cell_numbers_filtered <- as.data.frame(table(seurat_merge@meta.data$UPN_Cycle_Day_Sample_Type_Batch))
colnames(cell_numbers_filtered) <- c("Sample", "Filtered")

cell_numbers_all <- join(cell_numbers, cell_numbers_filtered)
cell_numbers_all$prop_keep <- cell_numbers_all$Filtered/cell_numbers_all$Unfiltered

cell_numbers_all <- cell_numbers_all %>% separate_wider_delim(Sample, delim = "_", names = c("UPN", "Cycle", "Day", "Sample_Type", "Batch"))
sum(cell_numbers_all$Unfiltered)
# 175829
sum(cell_numbers_all$Filtered)
# 143409


# Cell numbers of new samples
cell_numbers_all %>% filter(Batch %in% c("11", "2", "3", "4", "5", "6")) %>% 
  dplyr::select(Unfiltered) %>% sum()
# 111449
cell_numbers_all %>% filter(Batch %in% c("11", "2", "3", "4", "5", "6")) %>% 
  dplyr::select(Filtered) %>% sum()
# 104509

print(cell_numbers_all, n = nrow(cell_numbers_all))


#==============================================================================#
# Run SoupX ----
#==============================================================================#
# Extract product samples from these above objects
sample_list <- names(seurat_list_separated)

seurat_list_separated_SoupX <- sapply(sample_list,  function(i){
  print(i)
  # Read in count and droplet data
  d10x_toc <- Read10X(paste0("/scratch/aoill/projects/CAR-T/00_new_2025/soupX/demultiplexed_", i))
  
  # Need to read in batch specific empty droplet file
  upn_id <- str_split(i, "_")[[1]][1]
  sample_type_id <- str_split(i, "_")[[1]][2]
  cycle_type_id <- str_split(i, "_")[[1]][3]
  day_type_id <- str_split(i, "_")[[1]][4]
  batch_id <- str_split(i, "_")[[1]][5]
  # batches 37-39 CSF samples were not multiplexed
  if (batch_id %in% c("37", "38", "39")) {
    batch_path <- sample_metadata_not_multiplexed_csf %>% 
      filter(UPN == upn_id) %>%
      filter(Sample_Type == sample_type_id) %>%
      filter(Cycle == cycle_type_id) %>%
      filter(Day == day_type_id) %>%
      filter(Batch == batch_id) %>%
      pull(CellRanger_path) %>% unique()
    d10x_tod <- Read10X(paste0(batch_path, "/outs/raw_feature_bc_matrix/"))
  } else { # this is demultiplexed data
    batch_path <- sample_metadata_multiplexed_csf %>% filter(Batch == batch_id) %>%
      pull(CellRanger_path) %>% unique()
    d10x_tod <- Read10X(paste0(batch_path, "/outs/raw_feature_bc_matrix/"))
  }
  
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
#print(sample_list)


# Add metadata
#batch_ID = "FID_GEXFB"
cellRanger_path = "CellRanger_path"
for (i in names(seurat_list_separated_SoupX)) {
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
      seurat_list_separated_SoupX[[i]]@meta.data[[m]] = meta.data[[m]]
    }
  }
  seurat_list_separated_SoupX[[i]]@meta.data$Sample_ID <- paste(seurat_list_separated_SoupX[[i]]@meta.data$UPN, 
                                                                seurat_list_separated_SoupX[[i]]@meta.data$Sample_Type,
                                                                seurat_list_separated_SoupX[[i]]@meta.data$Cycle, 
                                                                seurat_list_separated_SoupX[[i]]@meta.data$Day, 
                                                                seurat_list_separated_SoupX[[i]]@meta.data$Batch,
                                                                sep = "_")
  
}
# probably also want to remove TCRpath


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
#DefaultAssay(seurat_list_separated_SoupX[[1]])
for (i in 1:length(seurat_list_separated_SoupX)) {
  DefaultAssay(seurat_list_separated_SoupX[[i]]) = "SoupX_RNA"
  seurat_list_separated_SoupX[[i]] = SCTransform(seurat_list_separated_SoupX[[i]], 
                                                 assay = "SoupX_RNA", 
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
csf_merge_rpca <- merge(x = seurat_list_separated_SoupX_batch_DoubletFinder_noRBSMT[[1]], 
                        y = seurat_list_separated_SoupX_batch_DoubletFinder_noRBSMT[2:length(seurat_list_separated_SoupX_batch_DoubletFinder_noRBSMT)])

# REMOVE DOUBLETS HERE
csf_merge_rpca_singlets <- subset(csf_merge_rpca, subset = doublet_finder == "Singlet")
table(csf_merge_rpca@meta.data$doublet_finder)
table(csf_merge_rpca_singlets@meta.data$doublet_finder)
csf_merge_rpca <- csf_merge_rpca_singlets

DefaultAssay(csf_merge_rpca) <- "SoupX_RNA"
csf_merge_rpca <- JoinLayers(csf_merge_rpca)
#csf_merge_rpca

# Split layers by RUN
csf_merge_rpca[["SoupX_RNA"]] <- split(csf_merge_rpca[["SoupX_RNA"]], f = csf_merge_rpca$UPN_Cycle_Day_Sample_Type_Batch)
#csf_merge_rpca

Sys.time()

# SCTransform is performed per sample, since they are split into layers
csf_merge_rpca <- SCTransform(csf_merge_rpca, assay = "SoupX_RNA", 
                              vst.flavor = "v2",
                              vars.to.regress = c("S.Score", "G2M.Score"),
                              return.only.var.genes = FALSE,
                              variable.features.n = 2500)
csf_merge_rpca <- RunPCA(csf_merge_rpca)
npcs <- min(get_pcs(csf_merge_rpca))
npcs # 18 (2500)
Sys.time()

# Integrate the layers using the SCT values
csf_merge_rpca <- IntegrateLayers(object = csf_merge_rpca, 
                                  method = RPCAIntegration,
                                  assay = "SCT", # either specify here or run default assay to SCT
                                  orig.reduction = "pca", 
                                  new.reduction = "integrated.rpca",
                                  verbose = FALSE,
                                  normalization.method = "SCT", 
                                  dims = 1:npcs)

# Find neighbors and clusters, and create UMAP
csf_merge_rpca <- FindNeighbors(csf_merge_rpca, dims = 1:npcs, verbose = FALSE, reduction = "integrated.rpca")
csf_merge_rpca <- RunUMAP(csf_merge_rpca, dims = 1:npcs, verbose = FALSE, 
                          reduction = "integrated.rpca", reduction.name = "umap.integrated.rpca")
csf_merge_rpca <- FindClusters(csf_merge_rpca, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), verbose = FALSE)
csf_merge_rpca@meta.data$UPN_Cycle <- paste(csf_merge_rpca@meta.data$UPN, csf_merge_rpca@meta.data$Cycle, sep = "_")


#==============================================================================#
# Subcluster ----
#==============================================================================#
# Subset cluster 7 (half looks like B cells, other half looks like DCs)
csf_merge_rpca_sbclt <- FindSubCluster(csf_merge_rpca, "7", resolution = 0.02, graph.name = "SCT_snn")

#==============================================================================#
# Lymphoid Integration ----
#==============================================================================#
# There are multiple passes of integration where I removed clusters that looked
# like doublets 

#------------------------------------------------------------------------------#
## Pass 1 ----
#------------------------------------------------------------------------------#
# Lymphoid: 0, 2, 3, 4, 5, 6, 7_1, 8, 11, 12, 15
DefaultAssay(csf_merge_rpca_sbclt) <- "SoupX_RNA"
csf_merge_rpca_sbclt_join <- JoinLayers(csf_merge_rpca_sbclt)

# we are using sub.cluster 
# csf_merge_rpca_sbclt_join
lymphoid <- c("0", "2", "3", "4", "5", "6", "7_1", "8", "11", "12", "15")
Idents(csf_merge_rpca_sbclt_join) <- "sub.cluster"
# get number of cells in each cluster
table(Idents(csf_merge_rpca_sbclt_join))
#0     8     2     5     4    15   7_1     9    11     3    14    13   7_0    10     1    12     6 
#35543  5137 12465  5376  8994   801  1498  3756  2263 12172  1051  1098  3641  2678 20128  1399  5187 

# Get cell names from the clusters that I think are T cells
csf_merge_rpca_lym <- subset(csf_merge_rpca_sbclt_join, idents = lymphoid)
table(csf_merge_rpca_lym@meta.data$sub.cluster)
#0    11    12    15     2     3     4     5     6   7_1     8 
#35543  2263  1399   801 12465 12172  8994  5376  5187  1498  5137 

# Drop levels
csf_merge_rpca_lym@meta.data[["sub.cluster"]] <- as.factor(csf_merge_rpca_lym@meta.data[["sub.cluster"]])
csf_merge_rpca_lym@meta.data[["sub.cluster"]] <- droplevels(csf_merge_rpca_lym@meta.data[["sub.cluster"]])

DefaultAssay(csf_merge_rpca_lym) <- "SoupX_RNA"
csf_merge_rpca_lym <- JoinLayers(csf_merge_rpca_lym)
#csf_merge_rpca_lym

# Split layers by RUN
csf_merge_rpca_lym[["SoupX_RNA"]] <- split(csf_merge_rpca_lym[["SoupX_RNA"]], f = csf_merge_rpca_lym$UPN_Cycle_Day_Sample_Type_Batch)
#csf_merge_rpca_lym


# SCTransform is performed per sample, since they are split into layers
csf_merge_rpca_lym <- SCTransform(csf_merge_rpca_lym, assay = "SoupX_RNA", 
                                  vst.flavor = "v2",
                                  vars.to.regress = c("S.Score", "G2M.Score"),
                                  return.only.var.genes = FALSE,
                                  variable.features.n = 2500)

csf_merge_rpca_lym <- RunPCA(csf_merge_rpca_lym)
npcs <- min(get_pcs(csf_merge_rpca_lym))
npcs # 14 (2500)
Sys.time()

# Integrate the layers using the SCT values
csf_merge_rpca_lym <- IntegrateLayers(object = csf_merge_rpca_lym, 
                                      method = RPCAIntegration,
                                      assay = "SCT", # either specify here or run default assay to SCT
                                      orig.reduction = "pca", 
                                      new.reduction = "integrated.rpca",
                                      verbose = FALSE,
                                      normalization.method = "SCT", 
                                      dims = 1:npcs)


# Find neighbors and clusters, and create UMAP
csf_merge_rpca_lym <- FindNeighbors(csf_merge_rpca_lym, dims = 1:npcs, verbose = FALSE, reduction = "integrated.rpca")
csf_merge_rpca_lym <- RunUMAP(csf_merge_rpca_lym, dims = 1:npcs, verbose = FALSE, 
                              reduction = "integrated.rpca", reduction.name = "umap.integrated.rpca")
csf_merge_rpca_lym <- FindClusters(csf_merge_rpca_lym, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), verbose = FALSE)


#------------------------------------------------------------------------------#
### Sub-cluster ----
#------------------------------------------------------------------------------#
# Possible doublet
Idents(csf_merge_rpca_lym) <- "SCT_snn_res.0.1"
csf_merge_rpca_lym_TCR_sbclt <- FindSubCluster(csf_merge_rpca_lym, "4", resolution = 0.1, graph.name = "SCT_snn")


#------------------------------------------------------------------------------#
## Pass 2 ----
#------------------------------------------------------------------------------#
# Remove cluster 4_0
lymphoid_pass2 <- c("0", "1", "2", "3", "4_1", "5", "6", "7", "8", "9")
Idents(csf_merge_rpca_lym_TCR_sbclt) <- "sub.cluster"
# get number of cells in each cluster
table(Idents(csf_merge_rpca_lym_TCR_sbclt))
#0     2     1     5     3     7     8   4_1     6   4_0     9 
#45768  8615  9315  5221  8193  2002  1521  1678  3430  3578  1514 

# Get cell names from the clusters that I think are T cells
csf_merge_rpca_lym_pass2 <- subset(csf_merge_rpca_lym_TCR_sbclt, idents = lymphoid_pass2)
table(csf_merge_rpca_lym_pass2@meta.data$sub.cluster)
#0     1     2     3   4_1     5     6     7     8     9 
#45768  9315  8615  8193  1678  5221  3430  2002  1521  1514 

# Drop levels
csf_merge_rpca_lym_pass2@meta.data[["sub.cluster"]] <- as.factor(csf_merge_rpca_lym_pass2@meta.data[["sub.cluster"]])
csf_merge_rpca_lym_pass2@meta.data[["sub.cluster"]] <- droplevels(csf_merge_rpca_lym_pass2@meta.data[["sub.cluster"]])

DefaultAssay(csf_merge_rpca_lym_pass2) <- "SoupX_RNA"
csf_merge_rpca_lym_pass2 <- JoinLayers(csf_merge_rpca_lym_pass2)
#csf_merge_rpca_lym_pass2

# Split layers by RUN
csf_merge_rpca_lym_pass2[["SoupX_RNA"]] <- split(csf_merge_rpca_lym_pass2[["SoupX_RNA"]], f = csf_merge_rpca_lym_pass2$UPN_Cycle_Day_Sample_Type_Batch)
#csf_merge_rpca_lym_pass2


# SCTransform is performed per sample, since they are split into layers
csf_merge_rpca_lym_pass2 <- SCTransform(csf_merge_rpca_lym_pass2, assay = "SoupX_RNA", 
                                  vst.flavor = "v2",
                                  vars.to.regress = c("S.Score", "G2M.Score"),
                                  return.only.var.genes = FALSE,
                                  variable.features.n = 2500)

csf_merge_rpca_lym_pass2 <- RunPCA(csf_merge_rpca_lym_pass2)
npcs <- min(get_pcs(csf_merge_rpca_lym_pass2))
npcs # 13 (2500)
Sys.time()

# Integrate the layers using the SCT values
csf_merge_rpca_lym_pass2 <- IntegrateLayers(object = csf_merge_rpca_lym_pass2, 
                                      method = RPCAIntegration,
                                      assay = "SCT", # either specify here or run default assay to SCT
                                      orig.reduction = "pca", 
                                      new.reduction = "integrated.rpca",
                                      verbose = FALSE,
                                      normalization.method = "SCT", 
                                      dims = 1:npcs)


# Find neighbors and clusters, and create UMAP
csf_merge_rpca_lym_pass2 <- FindNeighbors(csf_merge_rpca_lym_pass2, dims = 1:npcs, verbose = FALSE, reduction = "integrated.rpca")
csf_merge_rpca_lym_pass2 <- RunUMAP(csf_merge_rpca_lym_pass2, dims = 1:npcs, verbose = FALSE, 
                              reduction = "integrated.rpca", reduction.name = "umap.integrated.rpca")
csf_merge_rpca_lym_pass2 <- FindClusters(csf_merge_rpca_lym_pass2, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), verbose = FALSE)


#------------------------------------------------------------------------------#
### Sub-cluster ----
#------------------------------------------------------------------------------#
Idents(csf_merge_rpca_lym_pass2) <- "SCT_snn_res.0.1"
csf_merge_rpca_lym_pass2_sbclt_1 <- FindSubCluster(csf_merge_rpca_lym_pass2, "4", resolution = 0.1, graph.name = "SCT_snn")
Idents(csf_merge_rpca_lym_pass2_sbclt_1) <- "sub.cluster"
csf_merge_rpca_lym_pass2_sbclt_2 <- FindSubCluster(csf_merge_rpca_lym_pass2_sbclt_1, "4_1", resolution = 0.11, graph.name = "SCT_snn")
# Cluster 4_1_0 looks like a doublet


#------------------------------------------------------------------------------#
## Pass 3 ----
#------------------------------------------------------------------------------#
# Pass 3 integration with 4_1_0 removed and not including TRBV genes in variable features
# we are using sub.cluster 
# csf_merge_rpca_lym_pass2_sbclt_2
lymphoid_pass3 <- c("0", "1", "2", "3", "4_0", "4_1_1", "4_2", "5", "6", "7")
Idents(csf_merge_rpca_lym_pass2_sbclt_2) <- "sub.cluster"
# get number of cells in each cluster
table(Idents(csf_merge_rpca_lym_pass2_sbclt_2))
#0     3     2     1     6 4_1_0 4_1_1   4_0     5   4_2     7 
#46972  7977  8308 14118  1576  1210   589  2092  2379   585  1451 

# Get cell names from the clusters that I think are T cells
csf_merge_rpca_lym_pass3 <- subset(csf_merge_rpca_lym_pass2_sbclt_2, idents = lymphoid_pass3)
table(csf_merge_rpca_lym_pass3@meta.data$sub.cluster)
#0     1     2     3   4_0 4_1_1   4_2     5     6     7 
#46972 14118  8308  7977  2092   589   585  2379  1576  1451 

# Drop levels
csf_merge_rpca_lym_pass3@meta.data[["sub.cluster"]] <- as.factor(csf_merge_rpca_lym_pass3@meta.data[["sub.cluster"]])
csf_merge_rpca_lym_pass3@meta.data[["sub.cluster"]] <- droplevels(csf_merge_rpca_lym_pass3@meta.data[["sub.cluster"]])

DefaultAssay(csf_merge_rpca_lym_pass3) <- "SoupX_RNA"
csf_merge_rpca_lym_pass3 <- JoinLayers(csf_merge_rpca_lym_pass3)
#csf_merge_rpca_lym_pass3

# Split layers by RUN
csf_merge_rpca_lym_pass3[["SoupX_RNA"]] <- split(csf_merge_rpca_lym_pass3[["SoupX_RNA"]], f = csf_merge_rpca_lym_pass3$UPN_Cycle_Day_Sample_Type_Batch)
#csf_merge_rpca_lym_pass3


# SCTransform is performed per sample, since they are split into layers
csf_merge_rpca_lym_pass3 <- SCTransform(csf_merge_rpca_lym_pass3, assay = "SoupX_RNA", 
                                        vst.flavor = "v2",
                                        vars.to.regress = c("S.Score", "G2M.Score"),
                                        return.only.var.genes = FALSE,
                                        variable.features.n = 2500)


# Get variable features minus TRBV genes
orig_var_feat_tcells <- csf_merge_rpca_lym_pass3@assays[["SCT"]]@var.features
filt_var_feat_tcells <- grep("^TRBV", orig_var_feat_tcells, invert = TRUE, value = TRUE)

csf_merge_rpca_lym_pass3 <- RunPCA(csf_merge_rpca_lym_pass3, features = filt_var_feat_tcells)
npcs <- min(get_pcs(csf_merge_rpca_lym_pass3))
npcs # 19 (2500)
Sys.time()

# Integrate the layers using the SCT values
csf_merge_rpca_lym_pass3 <- IntegrateLayers(object = csf_merge_rpca_lym_pass3, 
                                            method = RPCAIntegration,
                                            assay = "SCT", # either specify here or run default assay to SCT
                                            orig.reduction = "pca", 
                                            new.reduction = "integrated.rpca",
                                            verbose = FALSE,
                                            normalization.method = "SCT", 
                                            dims = 1:npcs)


# Find neighbors and clusters, and create UMAP
csf_merge_rpca_lym_pass3 <- FindNeighbors(csf_merge_rpca_lym_pass3, dims = 1:npcs, verbose = FALSE, reduction = "integrated.rpca")
csf_merge_rpca_lym_pass3 <- RunUMAP(csf_merge_rpca_lym_pass3, dims = 1:npcs, verbose = FALSE, 
                                    reduction = "integrated.rpca", reduction.name = "umap.integrated.rpca")
csf_merge_rpca_lym_pass3 <- FindClusters(csf_merge_rpca_lym_pass3, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), verbose = FALSE)


#------------------------------------------------------------------------------#
### Sub-cluster ----
#------------------------------------------------------------------------------#
csf_merge_rpca_lym_pass3_join <- JoinLayers(csf_merge_rpca_lym_pass3)
csf_merge_rpca_sbclt <- FindSubCluster(csf_merge_rpca_lym_pass3_join, "6", resolution = 0.1, graph.name = "SCT_snn")


#------------------------------------------------------------------------------#
### Add cell type lables ----
#------------------------------------------------------------------------------#
csf_merge_rpca_sbclt@meta.data$ct_final <- NA # or any other initialization value
csf_merge_rpca_sbclt@meta.data$ct_final[csf_merge_rpca_sbclt@meta.data$sub.cluster == "0"] <- "T cells"
csf_merge_rpca_sbclt@meta.data$ct_final[csf_merge_rpca_sbclt@meta.data$sub.cluster == "1"] <- "T cells"
csf_merge_rpca_sbclt@meta.data$ct_final[csf_merge_rpca_sbclt@meta.data$sub.cluster == "2"] <- "T cells"
csf_merge_rpca_sbclt@meta.data$ct_final[csf_merge_rpca_sbclt@meta.data$sub.cluster == "3"] <- "T cells"
csf_merge_rpca_sbclt@meta.data$ct_final[csf_merge_rpca_sbclt@meta.data$sub.cluster == "4"] <- "T cells"
csf_merge_rpca_sbclt@meta.data$ct_final[csf_merge_rpca_sbclt@meta.data$sub.cluster == "7"] <- "T cells"
csf_merge_rpca_sbclt@meta.data$ct_final[csf_merge_rpca_sbclt@meta.data$sub.cluster == "8"] <- "T cells"
csf_merge_rpca_sbclt@meta.data$ct_final[csf_merge_rpca_sbclt@meta.data$sub.cluster == "9"] <- "T cells"
csf_merge_rpca_sbclt@meta.data$ct_final[csf_merge_rpca_sbclt@meta.data$sub.cluster == "5"] <- "NK"
csf_merge_rpca_sbclt@meta.data$ct_final[csf_merge_rpca_sbclt@meta.data$sub.cluster == "6_1"] <- "B"
csf_merge_rpca_sbclt@meta.data$ct_final[csf_merge_rpca_sbclt@meta.data$sub.cluster == "6_0"] <- "pDC"


#------------------------------------------------------------------------------#
### Add some metadata ----
#------------------------------------------------------------------------------#
### Radiographic response category (Response_cat) ----
csf_merge_rpca_sbclt@meta.data$Response_cat <- NA # or any other initialization value
csf_merge_rpca_sbclt@meta.data$Response_cat[csf_merge_rpca_sbclt@meta.data$UPN_Cycle == "514_2"] <- "Baseline"
csf_merge_rpca_sbclt@meta.data$Response_cat[csf_merge_rpca_sbclt@meta.data$UPN_Cycle == "514_5"] <- "Non-response"
csf_merge_rpca_sbclt@meta.data$Response_cat[csf_merge_rpca_sbclt@meta.data$UPN_Cycle == "514_8"] <- "Response"
csf_merge_rpca_sbclt@meta.data$Response_cat[csf_merge_rpca_sbclt@meta.data$UPN_Cycle == "514_11"] <- "Non-response"
csf_merge_rpca_sbclt@meta.data$Response_cat[csf_merge_rpca_sbclt@meta.data$UPN_Cycle == "514_13"] <- "Non-response"

csf_merge_rpca_sbclt@meta.data$Response_cat[csf_merge_rpca_sbclt@meta.data$UPN_Cycle == "515_4"] <- "Response"
csf_merge_rpca_sbclt@meta.data$Response_cat[csf_merge_rpca_sbclt@meta.data$UPN_Cycle == "515_8"] <- "Non-response"

csf_merge_rpca_sbclt@meta.data$Response_cat[csf_merge_rpca_sbclt@meta.data$UPN_Cycle == "574_1"] <- "Baseline"
csf_merge_rpca_sbclt@meta.data$Response_cat[csf_merge_rpca_sbclt@meta.data$UPN_Cycle == "574_3"] <- "Non-response"
csf_merge_rpca_sbclt@meta.data$Response_cat[csf_merge_rpca_sbclt@meta.data$UPN_Cycle == "574_8"] <- "Non-response"

csf_merge_rpca_sbclt@meta.data$Response_cat[csf_merge_rpca_sbclt@meta.data$UPN_Cycle == "689_1"] <- "Baseline"
csf_merge_rpca_sbclt@meta.data$Response_cat[csf_merge_rpca_sbclt@meta.data$UPN_Cycle == "689_4"] <- "Response"
csf_merge_rpca_sbclt@meta.data$Response_cat[csf_merge_rpca_sbclt@meta.data$UPN_Cycle == "689_8"] <- "Response"
csf_merge_rpca_sbclt@meta.data$Response_cat[csf_merge_rpca_sbclt@meta.data$UPN_Cycle == "689_12"] <- "Non-response"
csf_merge_rpca_sbclt@meta.data$Response_cat[csf_merge_rpca_sbclt@meta.data$UPN_Cycle == "689_17"] <- "Non-response"

csf_merge_rpca_sbclt@meta.data$Response_cat[csf_merge_rpca_sbclt@meta.data$UPN_Cycle == "692_2"] <- "Baseline"
csf_merge_rpca_sbclt@meta.data$Response_cat[csf_merge_rpca_sbclt@meta.data$UPN_Cycle == "692_4"] <- "Response"

csf_merge_rpca_sbclt@meta.data$Response_cat[csf_merge_rpca_sbclt@meta.data$UPN_Cycle == "705_1"] <- "Baseline"
csf_merge_rpca_sbclt@meta.data$Response_cat[csf_merge_rpca_sbclt@meta.data$UPN_Cycle == "705_5"] <- "Non-response"
csf_merge_rpca_sbclt@meta.data$Response_cat[csf_merge_rpca_sbclt@meta.data$UPN_Cycle == "705_8"] <- "Non-response"

csf_merge_rpca_sbclt@meta.data$Response_cat[csf_merge_rpca_sbclt@meta.data$UPN_Cycle == "716_2"] <- "Baseline"
csf_merge_rpca_sbclt@meta.data$Response_cat[csf_merge_rpca_sbclt@meta.data$UPN_Cycle == "716_4"] <- "Response"
csf_merge_rpca_sbclt@meta.data$Response_cat[csf_merge_rpca_sbclt@meta.data$UPN_Cycle == "716_8"] <- "Response"


#------------------------------------------------------------------------------#
### Add disease info ----
#------------------------------------------------------------------------------#
csf_merge_rpca_sbclt@meta.data$tumor_type <- NA # or any other initialization value
csf_merge_rpca_sbclt@meta.data$tumor_type[csf_merge_rpca_sbclt@meta.data$UPN == "514"] <- "Ependymoma"
csf_merge_rpca_sbclt@meta.data$tumor_type[csf_merge_rpca_sbclt@meta.data$UPN == "574"] <- "Ependymoma"
csf_merge_rpca_sbclt@meta.data$tumor_type[csf_merge_rpca_sbclt@meta.data$UPN == "689"] <- "Ependymoma"

csf_merge_rpca_sbclt@meta.data$tumor_type[csf_merge_rpca_sbclt@meta.data$UPN == "515"] <- "DMG"
csf_merge_rpca_sbclt@meta.data$tumor_type[csf_merge_rpca_sbclt@meta.data$UPN == "692"] <- "DMG"
csf_merge_rpca_sbclt@meta.data$tumor_type[csf_merge_rpca_sbclt@meta.data$UPN == "716"] <- "DMG"

csf_merge_rpca_sbclt@meta.data$tumor_type[csf_merge_rpca_sbclt@meta.data$UPN == "705"] <- "Glioblastoma"


#------------------------------------------------------------------------------#
### Add lymphodepletion status ----
#------------------------------------------------------------------------------#
csf_merge_rpca_sbclt@meta.data$Lymphodepletion <- NA # or any other initialization value
csf_merge_rpca_sbclt@meta.data$Lymphodepletion[csf_merge_rpca_sbclt@meta.data$UPN == "514"] <- "Non-Lymphodepleted"
csf_merge_rpca_sbclt@meta.data$Lymphodepletion[csf_merge_rpca_sbclt@meta.data$UPN == "515"] <- "Non-Lymphodepleted"
csf_merge_rpca_sbclt@meta.data$Lymphodepletion[csf_merge_rpca_sbclt@meta.data$UPN == "574"] <- "Non-Lymphodepleted" 
csf_merge_rpca_sbclt@meta.data$Lymphodepletion[csf_merge_rpca_sbclt@meta.data$UPN == "689"] <- "Lymphodepleted" 
csf_merge_rpca_sbclt@meta.data$Lymphodepletion[csf_merge_rpca_sbclt@meta.data$UPN == "705"] <- "Lymphodepleted" 
csf_merge_rpca_sbclt@meta.data$Lymphodepletion[csf_merge_rpca_sbclt@meta.data$UPN == "716"] <- "Lymphodepleted" 
csf_merge_rpca_sbclt@meta.data$Lymphodepletion[csf_merge_rpca_sbclt@meta.data$UPN == "692"] <- "Lymphodepleted" 


#------------------------------------------------------------------------------#
### Add CAR positivity status ----
#------------------------------------------------------------------------------#
# At least 3 reads covering construct (IL13OP) in T cells
DefaultAssay(csf_merge_rpca_sbclt) <- "RNA"
csf_merge_rpca_sbclt <- JoinLayers(csf_merge_rpca_sbclt)
csf_merge_rpca_sbclt@meta.data$IL13OPCounts <- csf_merge_rpca_sbclt[["RNA"]]$counts["IL13OP",]
csf_merge_rpca_sbclt@meta.data$CART <- ifelse((csf_merge_rpca_sbclt@meta.data$IL13OPCounts >= 3 & csf_merge_rpca_sbclt@meta.data$ct_final %in% c("T cells")), "Positive", "Negative")


#------------------------------------------------------------------------------#
## T cell integration ----
#------------------------------------------------------------------------------#
t_cell_clusters <- c("0", "1", "2", "3", "4", "7", "8", "9")
Idents(csf_merge_rpca_sbclt) <- "sub.cluster"
# get number of cells in each cluster
table(Idents(csf_merge_rpca_sbclt))
#0     3     4     5     2     1   6_1   6_0     8     9     7 
#24302  9054  8582  5497  9415 21288   859  1905  2341   345  2459

# Get cell names from the clusters that I think are T cells
csf_merge_rpca_t_cells <- subset(csf_merge_rpca_sbclt, idents = t_cell_clusters)
table(csf_merge_rpca_t_cells@meta.data$sub.cluster)
#0     1     2     3     4     7     8     9 
#24302 21288  9415  9054  8582  2459  2341   345

# Drop levels
csf_merge_rpca_t_cells@meta.data[["sub.cluster"]] <- as.factor(csf_merge_rpca_t_cells@meta.data[["sub.cluster"]])
csf_merge_rpca_t_cells@meta.data[["sub.cluster"]] <- droplevels(csf_merge_rpca_t_cells@meta.data[["sub.cluster"]])

DefaultAssay(csf_merge_rpca_t_cells) <- "SoupX_RNA"
csf_merge_rpca_t_cells <- JoinLayers(csf_merge_rpca_t_cells)
#csf_merge_rpca_t_cells

# Split layers by RUN
csf_merge_rpca_t_cells[["SoupX_RNA"]] <- split(csf_merge_rpca_t_cells[["SoupX_RNA"]], f = csf_merge_rpca_t_cells$UPN_Cycle_Day_Sample_Type_Batch)


# SCTransform is performed per sample, since they are split into layers
csf_merge_rpca_t_cells <- SCTransform(csf_merge_rpca_t_cells, assay = "SoupX_RNA", 
                                      vst.flavor = "v2",
                                      vars.to.regress = c("S.Score", "G2M.Score"),
                                      return.only.var.genes = FALSE,
                                      variable.features.n = 2500)


# Get variable features minus TRBV genes
orig_var_feat_tcells <- csf_merge_rpca_t_cells@assays[["SCT"]]@var.features
filt_var_feat_tcells <- grep("^TRBV", orig_var_feat_tcells, invert = TRUE, value = TRUE)

csf_merge_rpca_t_cells <- RunPCA(csf_merge_rpca_t_cells, features = filt_var_feat_tcells)
npcs <- min(get_pcs(csf_merge_rpca_t_cells))
npcs # 19 (2500)
Sys.time()

min(table(csf_merge_rpca_t_cells@meta.data$UPN_Cycle_Day_Sample_Type_Batch))

# Integrate the layers using the SCT values
csf_merge_rpca_t_cells <- IntegrateLayers(object = csf_merge_rpca_t_cells, 
                                          method = RPCAIntegration,
                                          assay = "SCT", # either specify here or run default assay to SCT
                                          orig.reduction = "pca", 
                                          new.reduction = "integrated.rpca",
                                          verbose = FALSE,
                                          normalization.method = "SCT", 
                                          dims = 1:npcs)


# Find neighbors and clusters, and create UMAP
csf_merge_rpca_t_cells <- FindNeighbors(csf_merge_rpca_t_cells, dims = 1:npcs, verbose = FALSE, reduction = "integrated.rpca")
csf_merge_rpca_t_cells <- RunUMAP(csf_merge_rpca_t_cells, dims = 1:npcs, verbose = FALSE, 
                                  reduction = "integrated.rpca", reduction.name = "umap.integrated.rpca")
csf_merge_rpca_t_cells <- FindClusters(csf_merge_rpca_t_cells, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), verbose = FALSE)


#------------------------------------------------------------------------------#
### Subcluster ----
#------------------------------------------------------------------------------#
# 0 and 3
Idents(csf_merge_rpca_t_cells) <- "SCT_snn_res.0.2"
csf_merge_rpca_t_cells_sbclt <- FindSubCluster(csf_merge_rpca_t_cells, "0", resolution = 0.11, graph.name = "SCT_snn")
csf_merge_rpca_t_cells_sbclt_join <- JoinLayers(csf_merge_rpca_t_cells_sbclt)
Idents(csf_merge_rpca_t_cells_sbclt_join) <- "sub.cluster"
csf_merge_rpca_t_cells_sbclt_join <- FindSubCluster(csf_merge_rpca_t_cells_sbclt_join, "3", resolution = 0.07, graph.name = "SCT_snn")


#------------------------------------------------------------------------------#
### Add T cell states to T cell object ----
#------------------------------------------------------------------------------#
csf_merge_rpca_t_cells_sbclt_join@meta.data$ct_final <- NA # or any other initialization value
csf_merge_rpca_t_cells_sbclt_join@meta.data$ct_final[csf_merge_rpca_t_cells_sbclt_join@meta.data$sub.cluster == "0_0"] <- "Memory"
csf_merge_rpca_t_cells_sbclt_join@meta.data$ct_final[csf_merge_rpca_t_cells_sbclt_join@meta.data$sub.cluster == "1"] <- "GAPDH+ Glycolysis"
csf_merge_rpca_t_cells_sbclt_join@meta.data$ct_final[csf_merge_rpca_t_cells_sbclt_join@meta.data$sub.cluster == "2"] <- "Effector-memory"
csf_merge_rpca_t_cells_sbclt_join@meta.data$ct_final[csf_merge_rpca_t_cells_sbclt_join@meta.data$sub.cluster == "3_0"] <- "Undifferentiated"
csf_merge_rpca_t_cells_sbclt_join@meta.data$ct_final[csf_merge_rpca_t_cells_sbclt_join@meta.data$sub.cluster == "3_2"] <- "Undifferentiated"
csf_merge_rpca_t_cells_sbclt_join@meta.data$ct_final[csf_merge_rpca_t_cells_sbclt_join@meta.data$sub.cluster == "4"] <- "Cytotoxic"
csf_merge_rpca_t_cells_sbclt_join@meta.data$ct_final[csf_merge_rpca_t_cells_sbclt_join@meta.data$sub.cluster == "6"] <- "Heat-shock"
csf_merge_rpca_t_cells_sbclt_join@meta.data$ct_final[csf_merge_rpca_t_cells_sbclt_join@meta.data$sub.cluster == "7"] <- "Proliferating"
csf_merge_rpca_t_cells_sbclt_join@meta.data$ct_final[csf_merge_rpca_t_cells_sbclt_join@meta.data$sub.cluster == "0_1"] <- "Undifferentiated"
csf_merge_rpca_t_cells_sbclt_join@meta.data$ct_final[csf_merge_rpca_t_cells_sbclt_join@meta.data$sub.cluster == "3_1"] <- "Undifferentiated"
csf_merge_rpca_t_cells_sbclt_join@meta.data$ct_final[csf_merge_rpca_t_cells_sbclt_join@meta.data$sub.cluster == "5"] <- "Treg"


#------------------------------------------------------------------------------#
## Add cell types and T cell states to the lymphoid object ----
#------------------------------------------------------------------------------#
csf_lym_meta <- csf_merge_rpca_sbclt@meta.data
csf_lym_meta$cellID <- rownames(csf_lym_meta)
csf_t_meta <- csf_merge_rpca_t_cells_sbclt_join@meta.data %>% dplyr::select(ct_final)
csf_t_meta$cellID <- rownames(csf_t_meta)
colnames(csf_t_meta) <- c("ct_final_T", "cellID")  


csf_lym_meta_join <- left_join(csf_lym_meta, csf_t_meta)
csf_merge_rpca_sbclt@meta.data <- csf_lym_meta_join
rownames(csf_merge_rpca_sbclt@meta.data) <- csf_lym_meta$cellID


csf_merge_rpca_sbclt@meta.data$ct_final_all <- csf_merge_rpca_sbclt@meta.data$ct_final # or any other initialization value
csf_merge_rpca_sbclt@meta.data$ct_final_all[csf_merge_rpca_sbclt@meta.data$ct_final_T == "Memory"] <- "T (Memory)"
csf_merge_rpca_sbclt@meta.data$ct_final_all[csf_merge_rpca_sbclt@meta.data$ct_final_T == "Undifferentiated"] <- "T (Undifferentiated)"
csf_merge_rpca_sbclt@meta.data$ct_final_all[csf_merge_rpca_sbclt@meta.data$ct_final_T == "Effector-memory"] <- "T (Effector-memory)"
csf_merge_rpca_sbclt@meta.data$ct_final_all[csf_merge_rpca_sbclt@meta.data$ct_final_T == "Cytotoxic"] <- "T (Cytotoxic)"
csf_merge_rpca_sbclt@meta.data$ct_final_all[csf_merge_rpca_sbclt@meta.data$ct_final_T == "Treg"] <- "Treg"
csf_merge_rpca_sbclt@meta.data$ct_final_all[csf_merge_rpca_sbclt@meta.data$ct_final_T == "GAPDH+ Glycolysis"] <- "T (GAPDH+ Glycolysis)"
csf_merge_rpca_sbclt@meta.data$ct_final_all[csf_merge_rpca_sbclt@meta.data$ct_final_T == "Heat-shock"] <- "T (Heat-shock)"
csf_merge_rpca_sbclt@meta.data$ct_final_all[csf_merge_rpca_sbclt@meta.data$ct_final_T == "Proliferating"] <- "T (Proliferating)"


#------------------------------------------------------------------------------#
## Add TCR data to lymphoid object ----
#------------------------------------------------------------------------------#
combined2 <- readRDS("/scratch/aoill/projects/CAR-T/00_new_2025/geo_2026/csf_scRepertoire_output.rds")


csf_lym_TCR_join_TCR <- combineExpression(combined2, csf_merge_rpca_sbclt,
                               group.by = "UPN_Cycle_Day_Sample_Type",
                               filterNA = FALSE,
                               proportion = FALSE,
                               cloneCall="strict",
                               cloneSize=c(Single=1, Small=5, Medium=20,
                                           Large=100, Hyperexpanded=500)
)
csf_lym_TCR_join_TCR@meta.data$cloneSize[is.na(csf_lym_TCR_join_TCR@meta.data$cloneSize)] <- "None ( < X <= 0)"



# resubset the lymphoid object to just the T cells to get the normalized frequency
t_cells <- c("T (Memory)", "T (Undifferentiated)", "T (Effector-memory)", 
             "T (Cytotoxic)", "Treg", "T (GAPDH+ Glycolysis)",
             "T (Heat-shock)", "T (Proliferating)")
csf_tcells <- subset(csf_lym_TCR_join_TCR, subset = ct_final_all %in% t_cells)


# calculate Frequency_norm
tmp_meta <- csf_tcells@meta.data

# Add a T cell number count for each UPN patient and cycle
tmp_meta <- csf_tcells@meta.data %>%
  #add_count(UPN, name = "UPN_n")
  add_count(UPN_Cycle, name = "UPN_Cycle_n")


# Get normalized TCR frequency
tmp_meta$Frequency_norm <- tmp_meta$clonalFrequency/tmp_meta$UPN_Cycle_n

csf_tcells@meta.data$UPN_Cycle_n <- tmp_meta$UPN_Cycle_n
csf_tcells@meta.data$Frequency_norm <- tmp_meta$Frequency_norm


tmp_meta_sub <- tmp_meta %>%
  dplyr::select(UPN, Sample_Type, Cycle, CART, CTstrict, clonalFrequency, #UPN_Cycle_n
                Frequency_norm)
rownames(tmp_meta) <- NULL
tmp_meta_sub <- na.omit(tmp_meta_sub)

# There could be duplicated rows in this metadata (clonotypes with freq >1 will 
# come up multiple times)
tmp_meta_sub_unique <- unique(tmp_meta_sub)

tmp_meta <- tmp_meta_sub_unique

# Remove any TCR with a frequency of 1
tmp_meta <- tmp_meta %>% dplyr::filter(clonalFrequency > 1)

# I need to double check that this is conceptually correct
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
csf_tcells@meta.data <- csf_tcells@meta.data %>% 
  mutate(Frequency_norm_cat = case_when(clonalFrequency == 1 ~ "Not Expanded",
                                        Frequency_norm < min(TCRs_top20$Frequency_norm) ~ "Not Expanded",
                                        (Frequency_norm >= min(TCRs_top20$Frequency_norm) & Frequency_norm < min(TCRs_top5$Frequency_norm)) ~ "Expanded (Top 20% - 5%)",
                                        (Frequency_norm >= min(TCRs_top5$Frequency_norm) & Frequency_norm < min(TCRs_top1$Frequency_norm)) ~ "More Expanded (Top 5% - 1%)",
                                        Frequency_norm >= min(TCRs_top1$Frequency_norm) ~ "Most Expanded (Top 1%)",
                                        
                                        TRUE ~ as.character(NA))) 
table(csf_tcells@meta.data$Frequency_norm_cat)


# Get expanded vs not expanded categories
csf_tcells@meta.data <- csf_tcells@meta.data %>% 
  mutate(Frequency_norm_cat_2 = case_when(clonalFrequency == 1 ~ "Not Expanded",
                                          Frequency_norm < min(TCRs_top20$Frequency_norm) ~ "Not Expanded",
                                          (Frequency_norm >= min(TCRs_top20$Frequency_norm) & Frequency_norm < min(TCRs_top5$Frequency_norm)) ~ "Expanded",
                                          (Frequency_norm >= min(TCRs_top5$Frequency_norm) & Frequency_norm < min(TCRs_top1$Frequency_norm)) ~ "Expanded",
                                          Frequency_norm >= min(TCRs_top1$Frequency_norm) ~ "Expanded",
                                          
                                          TRUE ~ as.character(NA))) 
table(csf_tcells@meta.data$Frequency_norm_cat_2)



# add in T cell TCR frequency back into the lymphoid object (will also need to
# add it to the T cell object TODO)
#colnames(csf_tcells@meta.data)
csf_tcells_tcr_meta <- csf_tcells@meta.data %>%
  dplyr::select(cellID, 
                CTgene, CTnt, CTaa, CTstrict, clonalProportion, clonalFrequency, 
                cloneSize, 
                Frequency_norm, Frequency_norm_cat, Frequency_norm_cat_2)
csf_tcells_meta_rnames <- rownames(csf_tcells@meta.data)

csf_lym_TCR_join_meta <- csf_lym_TCR_join@meta.data
csf_lym_TCR_join_meta_rnames <- rownames(csf_lym_TCR_join@meta.data)



dim(csf_lym_TCR_join_meta)
dim(csf_tcells_tcr_meta)
csf_lym_TCR_join_meta_l_join <- left_join(csf_lym_TCR_join_meta, csf_tcells_tcr_meta)
dim(csf_lym_TCR_join_meta_l_join) 


csf_lym_TCR_join@meta.data <- csf_lym_TCR_join_meta_l_join
rownames(csf_lym_TCR_join@meta.data) <- csf_lym_TCR_join_meta_rnames


# SAVE #
saveRDS(csf_lym_TCR_join, "/scratch/aoill/projects/CAR-T/00_new_2025/geo_2026/csf_lymphoid_seurat_obj.rds")


#==============================================================================#
# Myeloid Integration ----
#==============================================================================#

#------------------------------------------------------------------------------#
## Pass 1 ----
#------------------------------------------------------------------------------#
# we are using sub.cluster 
# csf_merge_rpca_sbclt_join
myeloid <- c("1", "7_0", "9", "10", "13", "14")
Idents(csf_merge_rpca_sbclt_join) <- "sub.cluster"
# get number of cells in each cluster
table(Idents(csf_merge_rpca_sbclt_join))
#0     8     2     5     4    15   7_1     9    11     3    14    13   7_0    10     1    12     6 
#35543  5137 12465  5376  8994   801  1498  3756  2263 12172  1051  1098  3641  2678 20128  1399  5187 

# Get cell names from the clusters that I think are T cells
csf_merge_rpca_mye <- subset(csf_merge_rpca_sbclt_join, idents = myeloid)
table(csf_merge_rpca_mye@meta.data$sub.cluster)
#1    10    13    14   7_0     9 
#20128  2678  1098  1051  3641  3756 

# Drop levels
csf_merge_rpca_mye@meta.data[["sub.cluster"]] <- as.factor(csf_merge_rpca_mye@meta.data[["sub.cluster"]])
csf_merge_rpca_mye@meta.data[["sub.cluster"]] <- droplevels(csf_merge_rpca_mye@meta.data[["sub.cluster"]])

DefaultAssay(csf_merge_rpca_mye) <- "SoupX_RNA"
csf_merge_rpca_mye <- JoinLayers(csf_merge_rpca_mye)
#csf_merge_rpca_mye

# Split layers by RUN
csf_merge_rpca_mye[["SoupX_RNA"]] <- split(csf_merge_rpca_mye[["SoupX_RNA"]], f = csf_merge_rpca_mye$UPN_Cycle_Day_Sample_Type_Batch)
#csf_merge_rpca_mye


# SCTransform is performed per sample, since they are split into layers
csf_merge_rpca_mye <- SCTransform(csf_merge_rpca_mye, assay = "SoupX_RNA", 
                                  vst.flavor = "v2",
                                  vars.to.regress = c("S.Score", "G2M.Score"),
                                  return.only.var.genes = FALSE,
                                  variable.features.n = 2500)

csf_merge_rpca_mye <- RunPCA(csf_merge_rpca_mye)
npcs <- min(get_pcs(csf_merge_rpca_mye))
npcs # 14 (2500)
Sys.time()

min(table(csf_merge_rpca_mye@meta.data$UPN_Cycle_Day_Sample_Type_Batch))
# 112

# Integrate the layers using the SCT values
csf_merge_rpca_mye <- IntegrateLayers(object = csf_merge_rpca_mye, 
                                      method = RPCAIntegration,
                                      assay = "SCT", # either specify here or run default assay to SCT
                                      orig.reduction = "pca", 
                                      new.reduction = "integrated.rpca",
                                      verbose = FALSE,
                                      normalization.method = "SCT", 
                                      dims = 1:npcs,
                                      k.weight = 90)


# Find neighbors and clusters, and create UMAP
csf_merge_rpca_mye <- FindNeighbors(csf_merge_rpca_mye, dims = 1:npcs, verbose = FALSE, reduction = "integrated.rpca")
csf_merge_rpca_mye <- RunUMAP(csf_merge_rpca_mye, dims = 1:npcs, verbose = FALSE, 
                              reduction = "integrated.rpca", reduction.name = "umap.integrated.rpca")
csf_merge_rpca_mye <- FindClusters(csf_merge_rpca_mye, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), verbose = FALSE)


#------------------------------------------------------------------------------#
### Subcluster ----
#------------------------------------------------------------------------------#
# Cluster 6
Idents(csf_merge_rpca_mye_TCR) <- "SCT_snn_res.0.2"
csf_merge_rpca_mye_TCR_sbclt <- FindSubCluster(csf_merge_rpca_mye_TCR, "6", resolution = 0.1, graph.name = "SCT_snn")


#------------------------------------------------------------------------------#
## Pass 2 ----
#------------------------------------------------------------------------------#
# Myeloid pass 2 with clusters 3, 6_1, and 9 removed
myeloid_pass2 <- c("0", "1", "2", "4", "5", "6_0", "7", "8")
Idents(csf_merge_rpca_mye_TCR_sbclt) <- "sub.cluster"
# get number of cells in each cluster
table(Idents(csf_merge_rpca_mye_TCR_sbclt))
#2   6_1     4     1     8     0   6_0     3     7     5     9 
#2876   457  2351  4397   609 15925   678  2658   674  1700    27

# Get cell names from the clusters that I think are T cells
csf_merge_rpca_mye_pass2 <- subset(csf_merge_rpca_mye_TCR_sbclt, idents = myeloid_pass2)
table(csf_merge_rpca_mye_pass2@meta.data$sub.cluster)
#0     1     2     4     5   6_0     7     8 
#15925  4397  2876  2351  1700   678   674   609 

# Drop levels
csf_merge_rpca_mye_pass2@meta.data[["sub.cluster"]] <- as.factor(csf_merge_rpca_mye_pass2@meta.data[["sub.cluster"]])
csf_merge_rpca_mye_pass2@meta.data[["sub.cluster"]] <- droplevels(csf_merge_rpca_mye_pass2@meta.data[["sub.cluster"]])

DefaultAssay(csf_merge_rpca_mye_pass2) <- "SoupX_RNA"
csf_merge_rpca_mye_pass2 <- JoinLayers(csf_merge_rpca_mye_pass2)
#csf_merge_rpca_mye_pass2

# Split layers by RUN
csf_merge_rpca_mye_pass2[["SoupX_RNA"]] <- split(csf_merge_rpca_mye_pass2[["SoupX_RNA"]], f = csf_merge_rpca_mye_pass2$UPN_Cycle_Day_Sample_Type_Batch)
#csf_merge_rpca_mye_pass2


# SCTransform is performed per sample, since they are split into layers
csf_merge_rpca_mye_pass2 <- SCTransform(csf_merge_rpca_mye_pass2, assay = "SoupX_RNA", 
                                        vst.flavor = "v2",
                                        vars.to.regress = c("S.Score", "G2M.Score"),
                                        return.only.var.genes = FALSE,
                                        variable.features.n = 2500)

csf_merge_rpca_mye_pass2 <- RunPCA(csf_merge_rpca_mye_pass2)
npcs <- min(get_pcs(csf_merge_rpca_mye_pass2))
npcs # 14 (2500)
Sys.time()

min(table(csf_merge_rpca_mye_pass2@meta.data$UPN_Cycle_Day_Sample_Type_Batch))
# 104

# Integrate the layers using the SCT values
csf_merge_rpca_mye_pass2 <- IntegrateLayers(object = csf_merge_rpca_mye_pass2, 
                                            method = RPCAIntegration,
                                            assay = "SCT", # either specify here or run default assay to SCT
                                            orig.reduction = "pca", 
                                            new.reduction = "integrated.rpca",
                                            verbose = FALSE,
                                            normalization.method = "SCT", 
                                            dims = 1:npcs,
                                            k.weight = 70)


# Find neighbors and clusters, and create UMAP
csf_merge_rpca_mye_pass2 <- FindNeighbors(csf_merge_rpca_mye_pass2, dims = 1:npcs, verbose = FALSE, reduction = "integrated.rpca")
csf_merge_rpca_mye_pass2 <- RunUMAP(csf_merge_rpca_mye_pass2, dims = 1:npcs, verbose = FALSE, 
                                    reduction = "integrated.rpca", reduction.name = "umap.integrated.rpca")
csf_merge_rpca_mye_pass2 <- FindClusters(csf_merge_rpca_mye_pass2, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), verbose = FALSE)


#------------------------------------------------------------------------------#
### Subcluster ----
#------------------------------------------------------------------------------#
DefaultAssay(csf_merge_rpca_mye_pass2) <- "SoupX_RNA"
csf_merge_rpca_mye_pass2_sbclt_join <- JoinLayers(csf_merge_rpca_mye_pass2)
Idents(csf_merge_rpca_mye_pass2_sbclt_join) <- "SCT_snn_res.0.2"
csf_merge_rpca_mye_pass2_sbclt_join_sbclt <- FindSubCluster(csf_merge_rpca_mye_pass2_sbclt_join, "6", resolution = 0.2, graph.name = "SCT_snn")
Idents(csf_merge_rpca_mye_pass2_sbclt_join_sbclt) <- "sub.cluster"
csf_merge_rpca_mye_pass2_sbclt_join_sbclt <- FindSubCluster(csf_merge_rpca_mye_pass2_sbclt_join_sbclt, "7", resolution = 0.1, graph.name = "SCT_snn")


#------------------------------------------------------------------------------#
## Pass 3 ----
#------------------------------------------------------------------------------#
# Myeloid pass 3 removing 6_1 and 7_1 
# csf_merge_rpca_mye_pass2_sbclt_join_sbclt
myeloid_pass3 <- c("0", "1", "2", "3", "4", "5", "6_0", "7_0", "7_2")
Idents(csf_merge_rpca_mye_pass2_sbclt_join_sbclt) <- "sub.cluster"
# get number of cells in each cluster
table(Idents(csf_merge_rpca_mye_pass2_sbclt_join_sbclt))
#2     3     1   7_0     0   7_1   6_0   6_1     5     4   7_2 
#4858  2177  6893   367 11711   198   428   259   688  1568    63 

csf_merge_rpca_mye_pass3 <- subset(csf_merge_rpca_mye_pass2_sbclt_join_sbclt, idents = myeloid_pass3)
table(csf_merge_rpca_mye_pass3@meta.data$sub.cluster)
#0     1     2     3     4     5   6_0   7_0   7_2 
#11711  6893  4858  2177  1568   688   428   367    63 

# Drop levels
csf_merge_rpca_mye_pass3@meta.data[["sub.cluster"]] <- as.factor(csf_merge_rpca_mye_pass3@meta.data[["sub.cluster"]])
csf_merge_rpca_mye_pass3@meta.data[["sub.cluster"]] <- droplevels(csf_merge_rpca_mye_pass3@meta.data[["sub.cluster"]])

DefaultAssay(csf_merge_rpca_mye_pass3) <- "SoupX_RNA"
csf_merge_rpca_mye_pass3 <- JoinLayers(csf_merge_rpca_mye_pass3)
#csf_merge_rpca_mye_pass3

# Split layers by RUN
csf_merge_rpca_mye_pass3[["SoupX_RNA"]] <- split(csf_merge_rpca_mye_pass3[["SoupX_RNA"]], f = csf_merge_rpca_mye_pass3$UPN_Cycle_Day_Sample_Type_Batch)
#csf_merge_rpca_mye_pass3


# SCTransform is performed per sample, since they are split into layers
csf_merge_rpca_mye_pass3 <- SCTransform(csf_merge_rpca_mye_pass3, assay = "SoupX_RNA", 
                                        vst.flavor = "v2",
                                        vars.to.regress = c("S.Score", "G2M.Score"),
                                        return.only.var.genes = FALSE,
                                        variable.features.n = 2500)

csf_merge_rpca_mye_pass3 <- RunPCA(csf_merge_rpca_mye_pass3)
npcs <- min(get_pcs(csf_merge_rpca_mye_pass3))
npcs # 14 (2500)
Sys.time()

min(table(csf_merge_rpca_mye_pass3@meta.data$UPN_Cycle_Day_Sample_Type_Batch))
# 96

# Integrate the layers using the SCT values
csf_merge_rpca_mye_pass3 <- IntegrateLayers(object = csf_merge_rpca_mye_pass3, 
                                            method = RPCAIntegration,
                                            assay = "SCT", # either specify here or run default assay to SCT
                                            orig.reduction = "pca", 
                                            new.reduction = "integrated.rpca",
                                            verbose = FALSE,
                                            normalization.method = "SCT", 
                                            dims = 1:npcs,
                                            k.weight = 70)


# Find neighbors and clusters, and create UMAP
csf_merge_rpca_mye_pass3 <- FindNeighbors(csf_merge_rpca_mye_pass3, dims = 1:npcs, verbose = FALSE, reduction = "integrated.rpca")
csf_merge_rpca_mye_pass3 <- RunUMAP(csf_merge_rpca_mye_pass3, dims = 1:npcs, verbose = FALSE, 
                                    reduction = "integrated.rpca", reduction.name = "umap.integrated.rpca")
csf_merge_rpca_mye_pass3 <- FindClusters(csf_merge_rpca_mye_pass3, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), verbose = FALSE)


### Add cell type lables ----
csf_merge_rpca_mye_pass3_join <- JoinLayers(csf_merge_rpca_mye_pass3)

csf_merge_rpca_mye_pass3_join@meta.data$ct_final <- NA

csf_merge_rpca_mye_pass3_join@meta.data$ct_final[
  csf_merge_rpca_mye_pass3_join@meta.data$SCT_snn_res.0.1 %in% c("0")] <- "Classical Monocyte"
csf_merge_rpca_mye_pass3_join@meta.data$ct_final[
  csf_merge_rpca_mye_pass3_join@meta.data$SCT_snn_res.0.1 %in% c("1")] <- "Mo/Mac"
csf_merge_rpca_mye_pass3_join@meta.data$ct_final[
  csf_merge_rpca_mye_pass3_join@meta.data$SCT_snn_res.0.1 %in% c("2")] <- "cDC2"
csf_merge_rpca_mye_pass3_join@meta.data$ct_final[
  csf_merge_rpca_mye_pass3_join@meta.data$SCT_snn_res.0.1 %in% c("3")] <- "Non-classical Monocyte"
csf_merge_rpca_mye_pass3_join@meta.data$ct_final[
  csf_merge_rpca_mye_pass3_join@meta.data$SCT_snn_res.0.1 %in% c("4")] <- "Macrophage"
csf_merge_rpca_mye_pass3_join@meta.data$ct_final[
  csf_merge_rpca_mye_pass3_join@meta.data$SCT_snn_res.0.1 %in% c("5")] <- "migDC"
csf_merge_rpca_mye_pass3_join@meta.data$ct_final[
  csf_merge_rpca_mye_pass3_join@meta.data$SCT_snn_res.0.1 %in% c("6")] <- "cDC1"


table(csf_merge_rpca_mye_pass3_join@meta.data$ct_final)

# save 
saveRDS(csf_merge_rpca_mye_pass3_join, "/scratch/aoill/projects/CAR-T/00_new_2025/geo_2026/csf_myeloid_seurat_obj.rds")



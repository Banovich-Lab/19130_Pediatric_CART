#==============================================================================#
# Authors: 
#     Angela M. Oill, aoill@tgen.org
#     Rishika Mudunuri, rmudunuri@tgen.org 
# Date: 2025
# Project: Pediatric CAR-T 
# Description: Preprocessing and integration of product samples
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
library(scGSVA)
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
sample_metadata_multiplexed <- sample_metadata %>% 
  filter(Multiplexed == "Yes") %>% filter(Batch == "37" | Batch == "38" |
                                            Batch == "39" | Batch == "40")

sample_metadata_not_multiplexed <- sample_metadata %>% 
  filter(Multiplexed == "No") %>% filter(Batch == "37" | Batch == "38" |
                                           Batch == "39" | Batch == "40")

#==============================================================================#
# Demultiplex ----
#==============================================================================#
seurat_list_demultiplexed <- prep_seurat_list_multiplexed(
  sample_metadata_multiplexed, 
  batch_ID = "Cell_Prefix", 
  cellRanger_path = "CellRanger_path", 
  cell_ID_prefix = "Cell_Prefix", 
  CellHashing_Ab = "CellHashing_Ab")


# Keep only singlets
seurat_list_demultiplexed_singlets <- filter_doublets_cellhashing(seurat_list_demultiplexed)

# Keep only product samples
product_seurat_list_demultiplexed_singlets <- list()
for (i in 1:length(seurat_list_demultiplexed_singlets)) {
  print(i)
  sub_i <- subset(seurat_list_demultiplexed_singlets[[i]], subset = Sample_Type == "Product") #No 'Product' in Batch 40
  product_seurat_list_demultiplexed_singlets <- c(product_seurat_list_demultiplexed_singlets, sub_i)
}
names(product_seurat_list_demultiplexed_singlets) <- c("Batch37_filtered", "Batch38_filtered", "Batch39_filtered")


#==============================================================================#
# Read in product sample that wasn't multiplexed (626) ----
#==============================================================================#
sample_metadata_not_multiplexed_626 <- sample_metadata_not_multiplexed %>% filter(FID_GEXFB == "F05565")
seurat_list_626 <- prep_seurat_list(sample_metadata_not_multiplexed_626, batch_ID = "FID_GEXFB", cellRanger_path = "CellRanger_path",
                                    cell_ID_prefix = "Cell_Prefix", run_soupX = FALSE)

product_seurat_list <- c(product_seurat_list_demultiplexed_singlets, seurat_list_626)


#==============================================================================#
# Filter ----
#==============================================================================#
product_seurat_list_filtered <- filter_manual(product_seurat_list, 
                                              pt_mt = 10, 
                                              nFeature = 1300, 
                                              nCount = 2500)

# Save these samples as 10X files
seurat_list_separated <- list()
for (i in 1:length(product_seurat_list_filtered)) {
  print(i)
  
  pedi_batch <- product_seurat_list_filtered[[i]]
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
  DropletUtils::write10xCounts(path = paste0("/scratch/rmudunuri/cart_project/Product/Fourth_version/soupX/demultiplexed_",names(seurat_list_separated[i])), 
                               x = obj.sub[["RNA"]]@layers$counts, 
                               barcodes = colnames(obj.sub[["RNA"]]), # cell names
                               gene.id = rownames(obj.sub[["RNA"]]), # Gene identifiers, one per row of X
                               type = "sparse", version="3", overwrite = T)
}


#==============================================================================#
# Run SoupX ----
#==============================================================================#
product_samples <- c("515_Product_NA_NA_37", "574_Product_NA_NA_38", "514_Product_NA_NA_38",
                     "625_Product_NA_NA_39", "626_Product_NA_NA_40")

#-----------------------------------------#
## Run SoupX on demultiplexed samples ----
#-----------------------------------------#
# Extract product samples from these above objects
sample_list <- product_samples
#sample_list <- names(seurat_list_separated)
product_seurat_list_separated_SoupX <- sapply(sample_list,  function(i){
  print(i)
  # Read in count and droplet data
  d10x_toc <- Read10X(paste0("/scratch/rmudunuri/cart_project/Product/Fourth_version/soupX/demultiplexed_", i))
  
  # Need to read in batch specific empty droplet file
  batch_id <- str_split(i, "_")[[1]][5]
  if ("37" %in% batch_id) {
    batch_path <- sample_metadata_multiplexed %>% filter(Batch == 37) %>%
      pull(CellRanger_path) %>% unique()
    d10x_tod <- Read10X(paste0(batch_path, "/outs/raw_feature_bc_matrix/"))
  } else if ("38" %in% batch_id) {
    batch_path <- sample_metadata_multiplexed %>% filter(Batch == 38) %>%
      pull(CellRanger_path) %>% unique()
    d10x_tod <- Read10X(paste0(batch_path, "/outs/raw_feature_bc_matrix/"))  
  } else if ("39" %in% batch_id) {
    batch_path <- sample_metadata_multiplexed %>% filter(Batch == 39) %>%
      pull(CellRanger_path) %>% unique()
    d10x_tod <- Read10X(paste0(batch_path, "/outs/raw_feature_bc_matrix/"))
  } else if ("40" %in% batch_id) {
    batch_path <- sample_metadata_not_multiplexed %>% filter(Batch == 40) %>% 
      filter(Sample_Type == "Product") %>% 
      pull(CellRanger_path) %>% unique()
    d10x_tod <- Read10X(paste0(batch_path, "/outs/raw_feature_bc_matrix/"))
  }
  
  if ("40" %in% batch_id) {
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
    
  } else {
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
    
  }
})
names(product_seurat_list_separated_SoupX) <- sample_list

# Add metadata
batch_ID = "FID_GEXFB"
cellRanger_path = "CellRanger_path"
for (i in names(product_seurat_list_separated_SoupX)) {
  upn_id <- str_split(i, "_")[[1]][1]
  sample_type_id <- str_split(i, "_")[[1]][2]
  meta.data <- sample_metadata %>% filter(UPN == upn_id) %>% filter(Sample_Type == sample_type_id)
  # Add sample metadata
  for (m in colnames(meta.data)){
    if (m == cellRanger_path) {
      next
    } else {
      product_seurat_list_separated_SoupX[[i]]@meta.data[[m]] = meta.data[[m]]
    }
  }
  product_seurat_list_separated_SoupX[[i]]@meta.data$Sample_ID <- paste(product_seurat_list_separated_SoupX[[i]]@meta.data$UPN, 
                                                                        product_seurat_list_separated_SoupX[[i]]@meta.data$Sample_Type,
                                                                        product_seurat_list_separated_SoupX[[i]]@meta.data$Cycle, 
                                                                        product_seurat_list_separated_SoupX[[i]]@meta.data$Day, 
                                                                        product_seurat_list_separated_SoupX[[i]]@meta.data$Batch,
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

product_seurat_list_separated_SoupX <- add_cell_cycle_score_2(product_seurat_list_separated_SoupX, rna_assay = "SoupX_RNA")


#-----------------------------------------#
## Re-normalize with cell cycle scores ----
#-----------------------------------------#
for (i in 1:length(product_seurat_list_separated_SoupX)) {
  DefaultAssay(product_seurat_list_separated_SoupX[[i]]) = "SoupX_RNA"
  product_seurat_list_separated_SoupX[[i]] = SCTransform(product_seurat_list_separated_SoupX[[i]], 
                                                         method = "glmGamPoi", 
                                                         vars.to.regress = c("S.Score", "G2M.Score"),
                                                         vst.flavor = "v2",
                                                         verbose = T)
}


#==============================================================================#
# Run DoubletFinder by batch ----
#==============================================================================#
# merge and re split the object by batch instead of each unique sample (UPN,cycle,day,sample type, batch)
product_seurat_list_separated_SoupX_merge <- merge(x = product_seurat_list_separated_SoupX[[1]], y = product_seurat_list_separated_SoupX[2:length(product_seurat_list_separated_SoupX)])
product_seurat_list_separated_SoupX_batch <- Seurat::SplitObject(product_seurat_list_separated_SoupX_merge, split.by = "Batch")


dblrate <- 0.15 # doublet rate

for (i in 1:length(product_seurat_list_separated_SoupX_batch)) {
  print(paste("Starting analysis on ",names(product_seurat_list_separated_SoupX_batch)[i], sep = ""))
  # Pre-process Seurat object (sctransform)
  DefaultAssay(product_seurat_list_separated_SoupX_batch[[i]]) <- "SoupX_RNA"
  product_seurat_list_separated_SoupX_batch[[i]] <- SCTransform(product_seurat_list_separated_SoupX_batch[[i]], 
                                                                assay = "SoupX_RNA", 
                                                                method = "glmGamPoi", 
                                                                vars.to.regress = c("S.Score", "G2M.Score"),
                                                                vst.flavor = "v2",
                                                                verbose = F)
  product_seurat_list_separated_SoupX_batch[[i]] <- RunPCA(product_seurat_list_separated_SoupX_batch[[i]])
  npcs <- min(get_pcs(product_seurat_list_separated_SoupX_batch[[i]]))
  product_seurat_list_separated_SoupX_batch[[i]] <- RunUMAP(product_seurat_list_separated_SoupX_batch[[i]], dims = 1:npcs)
  product_seurat_list_separated_SoupX_batch[[i]] <- FindNeighbors(product_seurat_list_separated_SoupX_batch[[i]], dims = 1:npcs, verbose = F)
  product_seurat_list_separated_SoupX_batch[[i]] <- FindClusters(product_seurat_list_separated_SoupX_batch[[i]], resolution = 1, verbose = F)
  
  
  # pK Identification (no ground-truth) 
  # For some reason, after running the preprocessing steps RNA and SoupX RNA 
  # gets split, not sure why this is happening but I am joining it all here
  # before running steps for doubletfinder
  DefaultAssay(product_seurat_list_separated_SoupX_batch[[i]]) <- "SoupX_RNA"
  product_seurat_list_separated_SoupX_batch[[i]] <- JoinLayers(product_seurat_list_separated_SoupX_batch[[i]])
  DefaultAssay(product_seurat_list_separated_SoupX_batch[[i]]) <- "RNA"
  product_seurat_list_separated_SoupX_batch[[i]] <- JoinLayers(product_seurat_list_separated_SoupX_batch[[i]])
  DefaultAssay(product_seurat_list_separated_SoupX_batch[[i]]) <- "SCT"
  sweep.res.list_seu_test <- paramSweep(product_seurat_list_separated_SoupX_batch[[i]], PCs = 1:npcs, sct = TRUE)
  sweep.stats_seu_test <- summarizeSweep(sweep.res.list_seu_test, GT = FALSE)
  bcmvn_seu_test <- find.pK(sweep.stats_seu_test)
  
  # Select optimal pK 
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn_seu_test[which.max(bcmvn_seu_test$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  
  # Homotypic Doublet Proportion Estimate 
  homotypic.prop <- modelHomotypic(product_seurat_list_separated_SoupX_batch[[i]]@meta.data$seurat_clusters)
  nExp_poi <- round(dblrate*nrow(product_seurat_list_separated_SoupX_batch[[i]]@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  
  # Run DoubletFinder
  product_seurat_list_separated_SoupX_batch[[i]] <- doubletFinder(seu = product_seurat_list_separated_SoupX_batch[[i]], 
                                                                  PCs = 1:npcs,
                                                                  pK = optimal.pk, 
                                                                  nExp = nExp_poi, 
                                                                  sct = TRUE
  )
  # Run DoubletFinder again w/ pANN
  pANN_col_pos <- ncol(product_seurat_list_separated_SoupX_batch[[i]]@meta.data) - 1
  product_seurat_list_separated_SoupX_batch[[i]] <- doubletFinder(seu = product_seurat_list_separated_SoupX_batch[[i]], 
                                                                  PCs = 1:npcs, 
                                                                  pK = optimal.pk, 
                                                                  nExp = nExp_poi.adj, 
                                                                  reuse.pANN = colnames(product_seurat_list_separated_SoupX_batch[[i]]@meta.data)[pANN_col_pos],
                                                                  sct = TRUE
  )
}

# Change colnames
for (i in 1:length(product_seurat_list_separated_SoupX_batch)) {
  pANN_col_pos <- ncol(product_seurat_list_separated_SoupX_batch[[i]]@meta.data) - 2
  names(product_seurat_list_separated_SoupX_batch[[i]]@meta.data)[pANN_col_pos] <- "pANN_col"
  names(product_seurat_list_separated_SoupX_batch[[i]]@meta.data)[pANN_col_pos+1] <- "DF.class_col"
  names(product_seurat_list_separated_SoupX_batch[[i]]@meta.data)[pANN_col_pos+2] <- "doublet_finder"
  
}



for (i in 1:length(product_seurat_list_separated_SoupX_batch)) {
  DefaultAssay(product_seurat_list_separated_SoupX_batch[[i]]) = "SoupX_RNA"
  product_seurat_list_separated_SoupX_batch[[i]] = SCTransform(product_seurat_list_separated_SoupX_batch[[i]], 
                                                               assay = "SoupX_RNA", 
                                                               method = "glmGamPoi", 
                                                               vars.to.regress = c("S.Score", "G2M.Score"),
                                                               vst.flavor = "v2",
                                                               verbose = T)
}


# then re-merge and re-split by how I will integrate
product_seurat_list_separated_SoupX_batch_merge <- merge(x = product_seurat_list_separated_SoupX_batch[[1]], 
                                                         y = product_seurat_list_separated_SoupX_batch[2:length(product_seurat_list_separated_SoupX_batch)]#,
                                                         #add.cell.ids = names(seurat_list_separated_SoupX_batch)
)

product_seurat_list_separated_SoupX_batch_merge@meta.data$UPN_Cycle_Day_Sample_Type_Batch <- paste(
  product_seurat_list_separated_SoupX_batch_merge@meta.data$UPN,
  product_seurat_list_separated_SoupX_batch_merge@meta.data$Cycle, 
  product_seurat_list_separated_SoupX_batch_merge@meta.data$Day, 
  product_seurat_list_separated_SoupX_batch_merge@meta.data$Sample_Type, 
  product_seurat_list_separated_SoupX_batch_merge@meta.data$Batch, 
  sep = "_")
product_seurat_list_separated_SoupX_batch_DoubletFinder <- Seurat::SplitObject(product_seurat_list_separated_SoupX_batch_merge, split.by = "UPN_Cycle_Day_Sample_Type_Batch")


#==============================================================================#
# Remove RB and MT genes ----
#==============================================================================#
# No doublet finder, new cell cycle scoring
product_seurat_list_separated_SoupX_batch_DoubletFinder_noRBSMT <- list()
message("Removing ribosomal and mitochondrial genes from Seurat object")
for (i in 1:length(product_seurat_list_separated_SoupX_batch_DoubletFinder)) {
  RBMTgenes <- grep(pattern = "^RP[SL]|^MRP[SL]|^MT-", x = rownames(product_seurat_list_separated_SoupX_batch_DoubletFinder[[i]]@assays$SoupX_RNA$counts), value = TRUE, invert = TRUE)
  product_seurat_list_separated_SoupX_batch_DoubletFinder_noRBSMT[[i]] = DietSeurat(product_seurat_list_separated_SoupX_batch_DoubletFinder[[i]], features = RBMTgenes, counts = TRUE, data = FALSE, scale.data = FALSE)
  
}
names(product_seurat_list_separated_SoupX_batch_DoubletFinder_noRBSMT) <- names(product_seurat_list_separated_SoupX_batch_DoubletFinder)


#==============================================================================#
# Integration ----
#==============================================================================#

#------------------------------------------------------------------------------#
## Integrate Product ----
#------------------------------------------------------------------------------#
# https://satijalab.org/seurat/articles/seurat5_integration
product_merge <- merge(x = product_seurat_list_separated_SoupX_batch_DoubletFinder_noRBSMT[[1]], 
                       y = product_seurat_list_separated_SoupX_batch_DoubletFinder_noRBSMT[2:length(product_seurat_list_separated_SoupX_batch_DoubletFinder_noRBSMT)])

# REMOVE DOUBLETS HERE
product_merge_singlets <- subset(product_merge, subset = doublet_finder == "Singlet")
table(product_merge@meta.data$doublet_finder)
table(product_merge_singlets@meta.data$doublet_finder)
product_merge <- product_merge_singlets

DefaultAssay(product_merge) <- "SoupX_RNA"
product_merge <- JoinLayers(product_merge)
product_merge

# Split layers by RUN
product_merge[["SoupX_RNA"]] <- split(product_merge[["SoupX_RNA"]], f = product_merge$UPN_Cycle_Day_Sample_Type_Batch)
Sys.time()


# SCTransform is performed per sample, since they are split into layers
product_merge <- SCTransform(product_merge, assay = "SoupX_RNA", 
                             vst.flavor = "v2",
                             vars.to.regress = c("S.Score", "G2M.Score"),
                             return.only.var.genes = FALSE,
                             variable.features.n = 2500)
product_merge <- RunPCA(product_merge)
npcs <- min(get_pcs(product_merge))
npcs # 11 (2500)
Sys.time()

# Integrate the layers using the SCT values
product_merge <- IntegrateLayers(object = product_merge, 
                                 method = RPCAIntegration,
                                 orig.reduction = "pca", 
                                 new.reduction = "integrated.rpca",
                                 verbose = FALSE,
                                 normalization.method = "SCT", 
                                 dims = 1:npcs)

# Find neighbors and clusters, and create UMAP (rerun after removing LYZ cluster)
product_merge <- FindNeighbors(product_merge, dims = 1:npcs, verbose = FALSE, reduction = "integrated.rpca")
product_merge <- RunUMAP(product_merge, dims = 1:npcs, verbose = FALSE, 
                         reduction = "integrated.rpca", reduction.name = "umap.integrated.rpca")
product_merge <- FindClusters(product_merge, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), verbose = FALSE)

#saveRDS(product_merge, "/scratch/rmudunuri/cart_project/Product/Fourth_version/product_merge_rpca.rds")

#Looking at LYZ cluster
FeaturePlot(product_merge, features = c("LYZ", "MS4A1", "CD19", "CD79A"), ncol = 4)

#Removing LYZ Cluster
product_merge_no_LYZ <- subset(product_merge, subset = SCT_snn_res.0.6 != "11")

#------------------------------------------------------------------------------#
## Re-integrating after dropping LYZ cluster ----
#------------------------------------------------------------------------------#
# Drop levels
product_merge_no_LYZ@meta.data[["SCT_snn_res.0.6"]] <- as.factor(product_merge_no_LYZ@meta.data[["SCT_snn_res.0.6"]])
product_merge_no_LYZ@meta.data[["SCT_snn_res.0.6"]] <- droplevels(product_merge_no_LYZ@meta.data[["SCT_snn_res.0.6"]])

DefaultAssay(product_merge_no_LYZ) <- "SoupX_RNA"
product_merge_no_LYZ <- JoinLayers(product_merge_no_LYZ)
product_merge_no_LYZ

# Split layers by RUN
product_merge_no_LYZ[["SoupX_RNA"]] <- split(product_merge_no_LYZ[["SoupX_RNA"]], f = product_merge_no_LYZ$UPN_Cycle_Day_Sample_Type_Batch)
product_merge_no_LYZ

product_merge_no_LYZ <- SCTransform(product_merge_no_LYZ, assay = "SoupX_RNA", 
                                    vst.flavor = "v2",
                                    vars.to.regress = c("S.Score", "G2M.Score"),
                                    return.only.var.genes = FALSE,
                                    variable.features.n = 2500)

product_merge_no_LYZ <- RunPCA(product_merge_no_LYZ)
npcs <- min(get_pcs(product_merge_no_LYZ))
npcs # 11 (2500)
Sys.time()

# Integrate the layers using the SCT values
product_merge_no_LYZ <- IntegrateLayers(object = product_merge_no_LYZ, 
                                        method = RPCAIntegration,
                                        orig.reduction = "pca", 
                                        new.reduction = "integrated.rpca",
                                        verbose = FALSE,
                                        normalization.method = "SCT", 
                                        dims = 1:npcs)

# Find neighbors and clusters, and create UMAP (rerun after removing LYZ cluster)
product_merge_no_LYZ <- FindNeighbors(product_merge_no_LYZ, dims = 1:npcs, verbose = FALSE, reduction = "integrated.rpca")
product_merge_no_LYZ <- RunUMAP(product_merge_no_LYZ, dims = 1:npcs, verbose = FALSE, 
                                reduction = "integrated.rpca", reduction.name = "umap.integrated.rpca")
product_merge_no_LYZ <- FindClusters(product_merge_no_LYZ, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), verbose = FALSE)

DimPlot(product_merge_no_LYZ, group.by = "SCT_snn_res.0.1", reduction = "umap.integrated.rpca")

#Looking at B cells cluster
FeaturePlot(product_merge_no_LYZ, features = c("LYZ", "MS4A1", "CD19", "CD79A"))
VlnPlot(product_merge_no_LYZ, features = c("LYZ", "MS4A1", "CD19", "CD79A"), pt.size = 0, ncol = 1)


#------------------------------------------------------------------------------#
## Re-integrating after dropping LYZ cluster without TRBV----
#------------------------------------------------------------------------------#
product_merge_no_LYZ_TRBV <- subset(product_merge, subset = SCT_snn_res.0.6 != "11")

# Drop levels
product_merge_no_LYZ_TRBV@meta.data[["SCT_snn_res.0.6"]] <- as.factor(product_merge_no_LYZ_TRBV@meta.data[["SCT_snn_res.0.6"]])
product_merge_no_LYZ_TRBV@meta.data[["SCT_snn_res.0.6"]] <- droplevels(product_merge_no_LYZ_TRBV@meta.data[["SCT_snn_res.0.6"]])

DefaultAssay(product_merge_no_LYZ_TRBV) <- "SoupX_RNA"
product_merge_no_LYZ_TRBV <- JoinLayers(product_merge_no_LYZ_TRBV)
product_merge_no_LYZ_TRBV

# Split layers by RUN
product_merge_no_LYZ_TRBV[["SoupX_RNA"]] <- split(product_merge_no_LYZ_TRBV[["SoupX_RNA"]], f = product_merge_no_LYZ_TRBV$UPN_Cycle_Day_Sample_Type_Batch)
product_merge_no_LYZ_TRBV

product_merge_no_LYZ_TRBV <- SCTransform(product_merge_no_LYZ_TRBV, assay = "SoupX_RNA", 
                                         vst.flavor = "v2",
                                         vars.to.regress = c("S.Score", "G2M.Score"),
                                         return.only.var.genes = FALSE,
                                         variable.features.n = 2500)
# Get variable features minus TRBV genes
orig_var_feat_no_LYZ_TRBV <- product_merge_no_LYZ_TRBV@assays[["SCT"]]@var.features
filtered_var_feat_no_LYZ_TRBV <- grep("^TRBV", orig_var_feat_no_LYZ_TRBV, invert = TRUE, value = TRUE)

# Use filtered list of variable features for PCA 
product_merge_no_LYZ_TRBV <- RunPCA(product_merge_no_LYZ_TRBV, features = filtered_var_feat_no_LYZ_TRBV)
npcs <- min(get_pcs(product_merge_no_LYZ_TRBV))
npcs # 12 (2458)
Sys.time()

# Integrate the layers using the SCT values and filtered variable features
product_merge_no_LYZ_TRBV <- IntegrateLayers(object = product_merge_no_LYZ_TRBV, 
                                             method = RPCAIntegration,
                                             orig.reduction = "pca", 
                                             new.reduction = "integrated.rpca",
                                             verbose = FALSE,
                                             normalization.method = "SCT", 
                                             dims = 1:npcs,
                                             features = filtered_var_feat_no_LYZ_TRBV
)

# Find neighbors and clusters, and create UMAP (rerun after removing LYZ cluster)
product_merge_no_LYZ_TRBV <- FindNeighbors(product_merge_no_LYZ_TRBV, dims = 1:npcs, verbose = FALSE, reduction = "integrated.rpca")
product_merge_no_LYZ_TRBV <- RunUMAP(product_merge_no_LYZ_TRBV, dims = 1:npcs, verbose = FALSE, 
                                     reduction = "integrated.rpca", reduction.name = "umap.integrated.rpca")
product_merge_no_LYZ_TRBV <- FindClusters(product_merge_no_LYZ_TRBV, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), verbose = FALSE)



product_merge_clean <- product_merge_no_LYZ_TRBV



#------------------------------------------------------------------------------#
## Add CART information to Immune object ----
#------------------------------------------------------------------------------#
DefaultAssay(product_merge_clean) <- "RNA"
product_merge_clean <- JoinLayers(product_merge_clean)

product_merge_clean@meta.data$IL13OPCounts <- product_merge_clean[["RNA"]]$counts["IL13OP",]
product_merge_clean@meta.data$CART <- ifelse(
  (product_merge_clean@meta.data$IL13OPCounts >= 3), "Positive", "Negative")

table(product_merge_clean@meta.data$CART)

DimPlot(product_merge_clean, reduction = "umap.integrated.rpca", group.by = "CART", cols = c("gray", "red"))



#==============================================================================#
# Add TCR data to T cell object ----
#==============================================================================#

#------------------------------------------------------------------------------#
## Add TCR data determined from scRepertoire ----
#------------------------------------------------------------------------------#

## Load the contig data ----
## Product
# One sample was not multiplexed
F05629 <- read.csv("/tgen_labs/banovich/SingleCell/CellRanger/5_0_0/Ensemble_98/PipelineData/Projects/BCTCSF/CellRangerOuts/VDJ/BCTCSF_0116_1_PB_Whole_C3_X5TCR_F05629_HNJG2DSX5/outs/filtered_contig_annotations.csv")

## PBMC and Product
# This data is multiplexed so I will have to use a function in 
# scRepetiore to create the list for each sample (patient/cycle/day combinations)

# Batch 37
F05037 <- read.csv("/tgen_labs/banovich/SingleCell/CellRanger/5_0_0/Ensemble_98/PipelineData/Projects/BCTCSF/CellRangerOuts/VDJ/BCTCSF_0103_1_PB_Whole_C1_X5TCR_F05037_HVFC5DSX3/outs/filtered_contig_annotations.csv")

# Batch 38
F05041 <- read.csv("/tgen_labs/banovich/SingleCell/CellRanger/5_0_0/Ensemble_98/PipelineData/Projects/BCTCSF/CellRangerOuts/VDJ/BCTCSF_0105_1_PB_Whole_C1_X5TCR_F05041_HVFC5DSX3/outs/filtered_contig_annotations.csv")

# Batch 39
F05625 <- read.csv("/tgen_labs/banovich/SingleCell/CellRanger/5_0_0/Ensemble_98/PipelineData/Projects/BCTCSF/CellRangerOuts/VDJ/BCTCSF_0114_1_PB_Whole_C1_X5TCR_F05625_HNJG2DSX5/outs/filtered_contig_annotations.csv")

# Batch 40
F05626 <- read.csv("/tgen_labs/banovich/SingleCell/CellRanger/5_0_0/Ensemble_98/PipelineData/Projects/BCTCSF/CellRangerOuts/VDJ/BCTCSF_0117_1_PB_Whole_C1_X5TCR_F05626_HNJG2DSX5/outs/filtered_contig_annotations.csv")

#------------------------------------------------------------------------------#
## Process the contig data ----
#------------------------------------------------------------------------------#
product_list <- SplitObject(product_merge_clean_join, split.by = "Batch")

# Strip the barcode extra labeling in product_list since the imported contig's 
# cell IDs are not formatted with a prefix (like Batch37_BARCODEID)
product_list_cells_renamed <- list()
for (i in names(product_list)){
  print(i)
  new_cell_IDs_tmp <- t(as.data.frame(str_split(rownames(product_list[[i]]@meta.data), "_")))
  new_cell_IDs <-  as.character(new_cell_IDs_tmp[,2])
  renamed_assay <- RenameCells(
    product_list[[i]],
    new.names = new_cell_IDs
  )
  product_list_cells_renamed <- c(product_list_cells_renamed, renamed_assay)
}
names(product_list_cells_renamed) <- names(product_list)

# Get the contig list for the multiplexed batches
contig_list_37 <- createHTOContigList(F05037, product_list_cells_renamed[["37"]], group.by = c("UPN", "Sample_Type"))
contig_list_38 <- createHTOContigList(F05041, product_list_cells_renamed[["38"]], group.by = c("UPN", "Sample_Type"))
contig_list_39 <- createHTOContigList(F05625, product_list_cells_renamed[["39"]], group.by = c("UPN", "Sample_Type"))
contig_list_40 <- createHTOContigList(F05629, product_list_cells_renamed[["40"]], group.by = c("UPN", "Sample_Type"))

# Merge all contigs and contig lists into one contig list
contig_list_product <- append(contig_list_37, contig_list_38)
contig_list_product <- append(contig_list_product, contig_list_39)
contig_list_product <- append(contig_list_product, contig_list_40)

names(contig_list_product)
# "515.Product" "574.Product" "514.Product" "625.Product" "626.Product"
names(contig_list_product) <- c("F05037_515_Product", "F05041_574_Product", 
                                "F05041_514_Product", "F05625_625_Product", 
                                "F05629_626_Product")

#------------------------------------------------------------------------------#
## Combine the contigs ----
#------------------------------------------------------------------------------#
# Will want to use the parameter removeNA set to true to remove any cell barcode 
# with an NA value in at least one of the chains
combined <- combineTCR(contig_list_product, 
                       samples =  c("F05037_515_Product", "F05041_574_Product", 
                                    "F05041_514_Product", "F05625_625_Product", 
                                    "F05629_626_Product"),
                       filterMulti = FALSE, removeNA = TRUE)

## Add TCR info to the Seurat objects ----

### Add a new barcode column in the TCR list to match the Seurat objects
# I need to replace the barcode column in combined which matching barcode to 
# seurat object 
# Path to google sheet with metadata for your samples
#metapath <- "https://docs.google.com/spreadsheets/d/1fmEuOFXm893mfS2T1zpEsMrN7MpFc_Aym24ZODbw2bA/edit#gid=0"
#sheet_name <- "Sheet1"

sample_metadata <- read_csv("/home/rmudunuri/cart_project/19130_all_batches_metadata.csv")

combined2 <- lapply(names(combined), function(i){
  j <- combined[[i]]
  k <- unique(j$sample)
  l <- strsplit(k, "_")[[1]][1] # this is the TCRseq_ID in my metadata file
  # add new barcode column
  # also make an original barcode column (without anything in front of barcode)
  j$tcr_barcode <- j$barcode
  j$orig_barcode <- sapply(strsplit(j$barcode, "_"), `[`, 4)
  cell_id_prefix <- sample_metadata %>% filter(TCRseq_ID == l) %>% pull(Cell_Prefix) %>% unique()
  j$barcode <- paste0(cell_id_prefix, "_", j$orig_barcode)
  j$TCR_ID <- sapply(strsplit(j$sample, "_"), `[`, 1)
  j$UPN <- sapply(strsplit(j$sample, "_"), `[`, 2)
  j$Cycle <- "NA"
  j$Day <- "NA"
  j$Sample_Type <- sapply(strsplit(j$sample, "_"), `[`, 3)
  j$UPN_Cycle_Day_Sample_Type <- paste(j$UPN, j$Cycle, j$Day, j$Sample_Type, sep = "_")
  j$TCR_ID_UPN_Cycle_Day_Sample_Type <- paste(j$TCR_ID, j$UPN, j$Cycle, j$Day, j$Sample_Type, sep = "_")
  j
})
names(combined2) <- c("F05037_515_Product", "F05041_574_Product", 
                      "F05041_514_Product", "F05625_625_Product", 
                      "F05629_626_Product")

saveRDS(combined2, "/scratch/aoill/projects/CAR-T/00_new_2025/geo_2026/product_scRepertiore_output.rds")


#==============================================================================#
# Add cell typing info ----
#==============================================================================#
product_fnl_TCR@meta.data$ct_final_all <- NA # or any other initialization value
product_fnl_TCR@meta.data$ct_final_all[product_fnl_TCR@meta.data$SCT_snn_res.0.3 == "0"] <- "Undifferentiated"
product_fnl_TCR@meta.data$ct_final_all[product_fnl_TCR@meta.data$SCT_snn_res.0.3 == "1"] <- "Effector"
product_fnl_TCR@meta.data$ct_final_all[product_fnl_TCR@meta.data$SCT_snn_res.0.3 == "2"] <- "Undifferentiated"
product_fnl_TCR@meta.data$ct_final_all[product_fnl_TCR@meta.data$SCT_snn_res.0.3 == "3"] <- "Activated"
product_fnl_TCR@meta.data$ct_final_all[product_fnl_TCR@meta.data$SCT_snn_res.0.3 == "4"] <- "Effector"
product_fnl_TCR@meta.data$ct_final_all[product_fnl_TCR@meta.data$SCT_snn_res.0.3 == "5"] <- "Proliferating"
product_fnl_TCR@meta.data$ct_final_all[product_fnl_TCR@meta.data$SCT_snn_res.0.3 == "6"] <- "Proliferating"
product_fnl_TCR@meta.data$ct_final_all[product_fnl_TCR@meta.data$SCT_snn_res.0.3 == "7"] <- "CAR high"


#==============================================================================#
# Add more metadata ----
#==============================================================================#
# add disease info
product_fnl_TCR@meta.data$tumor_type <- NA # or any other initialization value
product_fnl_TCR@meta.data$tumor_type[product_fnl_TCR@meta.data$UPN == "514"] <- "Ependymoma"
product_fnl_TCR@meta.data$tumor_type[product_fnl_TCR@meta.data$UPN == "574"] <- "Ependymoma"
product_fnl_TCR@meta.data$tumor_type[product_fnl_TCR@meta.data$UPN == "625"] <- "Ependymoma" 

product_fnl_TCR@meta.data$tumor_type[product_fnl_TCR@meta.data$UPN == "515"] <- "DMG"
product_fnl_TCR@meta.data$tumor_type[product_fnl_TCR@meta.data$UPN == "626"] <- "DMG" 


saveRDS(product_fnl_TCR, "/scratch/rmudunuri/cart_project/Product/Fourth_version/product_fnl_TCR.rds")

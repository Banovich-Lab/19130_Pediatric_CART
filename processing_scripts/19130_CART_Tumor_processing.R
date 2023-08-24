#==============================================================================#
# Author : Angela M. Oill, aoill@tgen.org
# Date: 2023/06/09
# Project: Pediatric CAR-T 
# Description: Preprocessing and integration of tumor samples
#==============================================================================#

#==============================================================================#
# Load Libraries, helper modules and functions ----
#==============================================================================#
#library(Matrix)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(DropletUtils)

source("/home/aoill/projects/SingleCellBestPractices/scripts/preprocessing_qc_module.R")
source("/home/aoill/projects/SingleCellBestPractices/scripts/helper_functions_module.R")


#==============================================================================#
# Set variables ----
#==============================================================================#
set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)

# Path to google sheet with metadata for your samples
metapath <- "https://docs.google.com/spreadsheets/d/1gRL53qgRBApRHxznwTK_17I1cRlAL0GGf8nMikOJle0/edit#gid=0"
sheet_name <- "Sheet1"

# Output directory path
out_dir <- "/scratch/aoill/projects/CAR-T/00_final/"

#==============================================================================#
# Prep metadata ----
#==============================================================================#

#--------------------#
## Read metadata ----
#--------------------#
sample_metadata <- load_metadata(metadata_path = metapath, data_type = "googlesheet", sheet_name = sheet_name)

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
  DropletUtils::write10xCounts(path = paste0("/scratch/aoill/projects/CAR-T/00_final/soupX/demultiplexed_",names(seurat_list_separated[i])), 
                               x = obj.sub[["RNA"]]@data, 
                               barcodes = colnames(obj.sub[["RNA"]]@data), # cell names
                               gene.id = rownames(obj.sub[["RNA"]]@data), # Gene identifiers, one per row of X
                               type = "sparse", version="3", overwrite = T)
}


#-----------------------------------------------------#
## Get number of cells before and after filtering ----
#-----------------------------------------------------#
seurat_merge <- merge(x = seurat_list[[1]], y = seurat_list[2:length(seurat_list)])
seurat_merge@meta.data$UPN_Cycle_Day_Sample_Type_Batch <- paste0(seurat_merge@meta.data$UPN, "_",
                                                                 #seurat_merge@meta.data$Cycle,  "_", 
                                                                 #seurat_merge@meta.data$Day, "_",
                                                                 seurat_merge@meta.data$Sample_Type_More, "_",
                                                                 seurat_merge@meta.data$Batch)
table(seurat_merge@meta.data$UPN_Cycle_Day_Sample_Type_Batch)

seurat_merge <- merge(x = seurat_list_filtered[[1]], y = seurat_list_filtered[2:length(seurat_list_filtered)])
seurat_merge@meta.data$UPN_Cycle_Day_Sample_Type_Batch <- paste0(seurat_merge@meta.data$UPN, "_",
                                                                 #seurat_merge@meta.data$Cycle,  "_", 
                                                                 #seurat_merge@meta.data$Day, "_",
                                                                 seurat_merge@meta.data$Sample_Type_More, "_",
                                                                 seurat_merge@meta.data$Batch)
table(seurat_merge@meta.data$UPN_Cycle_Day_Sample_Type_Batch)


#==============================================================================#
# Run SoupX ----
#==============================================================================#
# Extract product samples from these above objects
sample_list <- names(seurat_list_separated)

seurat_list_separated_SoupX <- sapply(sample_list,  function(i){
  print(i)
  # Read in count and droplet data
  d10x_toc <- Read10X(paste0("/scratch/aoill/projects/CAR-T/00_final/soupX/demultiplexed_", i))
  
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
    sample_seurat_list[[i]] <- SCTransform(sample_seurat_list[[i]], 
                                           vst.flavor = "v2", verbose = F)
    DefaultAssay(sample_seurat_list[[i]]) <- "SCT"
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
# Remove RB and MT genes ----
#==============================================================================#
# Remove ribosomal and mt genes
seurat_list_separated_SoupX_noRBSMT <- list()
message("Removing ribosomal and mitochondrial genes from Seurat object")
for (i in 1:length(seurat_list_separated_SoupX)) {
  RBMTgenes <- grep(pattern = "^RP[SL]|^MRP[SL]|^MT-", x = rownames(seurat_list_separated_SoupX[[i]]@assays$RNA@data), value = TRUE, invert = TRUE)
  seurat_list_separated_SoupX_noRBSMT[[i]] = subset(seurat_list_separated_SoupX[[i]], features = RBMTgenes)
  
}
names(seurat_list_separated_SoupX_noRBSMT) <- names(seurat_list_separated_SoupX)


#==============================================================================#
# Integration ----
#==============================================================================#

#------------------------------------------------------------------------------#
## Add some more metadata ----
#------------------------------------------------------------------------------#
for (i in 1:length(seurat_list_separated_SoupX_noRBSMT)) {
  seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$UPN_Cycle_Day_Sample_Type_Batch <- paste0(seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$UPN, "_",
                                                                                               #seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$Cycle,  "_", 
                                                                                               #seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$Day, "_",
                                                                                               "NA_", 
                                                                                               "NA_",
                                                                                               seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$Sample_Type_More, "_",
                                                                                               seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$Batch)

}



#------------------------------------------------------------------------------#
## Integrate Tumor ----
#------------------------------------------------------------------------------#
seurat_list <- seurat_list_separated_SoupX_noRBSMT

seurat_list <- lapply(seurat_list, function(xx){
  # Rerunning SCTransform
  DefaultAssay(xx) <- "SoupX_RNA"
  xx <- SCTransform(xx, 
                    method = "glmGamPoi",
                    vars.to.regress = c("S.Score", "G2M.Score"), 
                    vst.flavor = "v2",
                    verbose = T)
  
})

integrated_obj <- basic_sct_rpca_integration(seurat_list, npcs=50, k_weight=100, numfeatures = 2500)



integrated_obj@meta.data$UPN_Sample_Type_More_Batch <- paste(integrated_obj@meta.data$UPN,
                                                             integrated_obj@meta.data$Sample_Type_More,
                                                             integrated_obj@meta.data$Batch)

## SAVE ##
saveRDS(integrated_obj, "/scratch/aoill/projects/CAR-T/00_final/rds_files/Tumor_integrated_20230512.rds")
#integrated_obj <- readRDS("/scratch/aoill/projects/CAR-T/00_final/rds_files/Tumor_integrated_20230512.rds")


#==============================================================================#
# Add immune sub-setting reintegrating, and T cell re clustering (TO ADD) ----
#==============================================================================#
# TODO

#==============================================================================#
# NEW Subset immune cells and re-integrate ----
#==============================================================================#
DimPlot(integrated_obj, group.by = "integrated_sct_snn_res.0.2", 
        reduction = "integrated_sct_umap",
        label = T) 

# Immune clusters
immune_clusters <- c("3", "7", "9")

# we are using resolution 0.2 
Idents(integrated_obj) <- "integrated_sct_snn_res.0.2"
# get number of cells in each cluster
table(Idents(integrated_obj))
# Get cell names from the clusters that I think are T cells
integrated_obj_immune_clusters <- subset(integrated_obj, idents = immune_clusters)
table(integrated_obj_immune_clusters@meta.data$integrated_sct_snn_res.0.2)

# Not sure if I need to do this but drop levels
integrated_obj_immune_clusters@meta.data[["integrated_sct_snn_res.0.2"]] <- droplevels(integrated_obj_immune_clusters@meta.data[["integrated_sct_snn_res.0.2"]])
table(integrated_obj_immune_clusters@meta.data$integrated_sct_snn_res.0.2)

# Split object by Cell_Prefix, since this is how they were split originally 
integrated_obj_immune_clusters_batches_list <- SplitObject(integrated_obj_immune_clusters, split.by = "Cell_Prefix")
names(integrated_obj_immune_clusters_batches_list) <- names(seurat_list_separated_SoupX_noRBSMT)
names(integrated_obj_immune_clusters_batches_list)

# Remove the CD45 negative samples
integrated_obj_immune_clusters_batches_list_sub <- c(integrated_obj_immune_clusters_batches_list["574_Tumor_Cells_CD45_positive_39"],
                                                     integrated_obj_immune_clusters_batches_list["625_Tumor_Cells_myelin_negative_40"],
                                                     integrated_obj_immune_clusters_batches_list["625_Tumor_Cells_CD45_positive_NA"])
# Re-run SCTransform
# Has to iterate by each item in list
integrated_obj_immune_clusters_batches_list_sub <- lapply(integrated_obj_immune_clusters_batches_list_sub, function(xx){
  # Rerunning SCTransform
  DefaultAssay(xx) <- "SoupX_RNA"
  xx <- SCTransform(xx, 
                    method = "glmGamPoi",
                    vars.to.regress = c("S.Score", "G2M.Score"), 
                    vst.flavor = "v2",
                    verbose = T)
  
})


# Integrate
immune_integrated_obj <- basic_sct_rpca_integration(integrated_obj_immune_clusters_batches_list_sub, npcs=50, k_weight=100, numfeatures = 2500)


p1 <- DimPlot(immune_integrated_obj, group.by = "integrated_sct_snn_res.0.2", 
              reduction = "integrated_sct_umap",
              label = T) 
p2 <- DimPlot(immune_integrated_obj, group.by = "UPN", reduction = "integrated_sct_umap", label = F) 
p2 <- DimPlot(immune_integrated_obj, group.by = "UPN_Tumor_info", reduction = "integrated_sct_umap", label = F) 
p3 <- DimPlot(immune_integrated_obj, group.by = "Phase", reduction = "integrated_sct_umap", label = F) 

p4 <- DimPlot(immune_integrated_obj, group.by = "CART") + 
  scale_color_manual(values=c("grey", "red")) +
  ggtitle("CAR-T (+/-)") +
  theme(legend.position="right") 

ggarrange(p1, p2,
          p3, p4,
          ncol = 2,
          nrow = 2, align = "hv")


# Sub-cluster 0
Idents(immune_integrated_obj) <- "integrated_sct_snn_res.0.2"
immune_integrated_obj_subcluster_0 <- FindSubCluster(immune_integrated_obj, "0", resolution = .1, graph.name = "integrated_sct_snn")
DimPlot(immune_integrated_obj_subcluster_0, group.by = "sub.cluster", label = T)

VlnPlot(immune_integrated_obj_subcluster_0, features = feature_genes, 
        group.by = "sub.cluster",
        pt.size = 0,
        ncol = 2)

Idents(immune_integrated_obj_subcluster_0) <- "sub.cluster"
immune_integrated_obj_subcluster_0_1 <- FindSubCluster(immune_integrated_obj_subcluster_0, "1", resolution = .3, graph.name = "integrated_sct_snn")
DimPlot(immune_integrated_obj_subcluster_0_1, group.by = "sub.cluster", label = T)

# Draft Label cell types
Idents(immune_integrated_obj_subcluster_0_1) <- "sub.cluster"

immune_integrated_obj_subcluster_0_1@meta.data$CT <- immune_integrated_obj_subcluster_0_1@meta.data$sub.cluster

immune_integrated_obj_subcluster_0_1@meta.data$CT[
  immune_integrated_obj_subcluster_0_1@meta.data$sub.cluster %in% c("2")] <- "T cells"
immune_integrated_obj_subcluster_0_1@meta.data$CT[
  immune_integrated_obj_subcluster_0_1@meta.data$sub.cluster %in% c("4")] <- "Monocytes"
immune_integrated_obj_subcluster_0_1@meta.data$CT[
  immune_integrated_obj_subcluster_0_1@meta.data$sub.cluster %in% c("0_0", "3", "6", "1_1")] <- "Macrophages"
immune_integrated_obj_subcluster_0_1@meta.data$CT[
  immune_integrated_obj_subcluster_0_1@meta.data$sub.cluster %in% c("0_1", "5")] <- "cDCs"
immune_integrated_obj_subcluster_0_1@meta.data$CT[
  immune_integrated_obj_subcluster_0_1@meta.data$sub.cluster %in% c("7")] <- "pDC"
immune_integrated_obj_subcluster_0_1@meta.data$CT[
  immune_integrated_obj_subcluster_0_1@meta.data$sub.cluster %in% c("1_0", "1_2")] <- "Myeloid"


table(immune_integrated_obj_subcluster_0_1@meta.data$CT)

DimPlot(immune_integrated_obj_subcluster_0_1, group.by = "CT", label = T)


## SAVE ##
#immune_integrated_obj_subcluster_0_1 <- readRDS("/scratch/aoill/projects/CAR-T/00_final/rds_files/Tumor_integrated_immune_cells_20230515.rds")



#==============================================================================#
# Add metadata (TO ADD) ----
#==============================================================================#
# Add IL13OP info to metadata
DefaultAssay(integrated_obj) <- "RNA"

integrated_obj@meta.data$IL13OPCounts <- integrated_obj@assays$RNA@counts["IL13OP",]
integrated_obj@meta.data$CART <- ifelse(integrated_obj@meta.data$IL13OPCounts >= 3, "Positive", "Negative")
unique(integrated_obj@meta.data$CART)
table(integrated_obj@meta.data$CART)
#Negative 
# 6119
# Negative 
# 16687
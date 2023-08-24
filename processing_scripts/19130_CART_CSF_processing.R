#==============================================================================#
# Author : Angela M. Oill, aoill@tgen.org
# Date: 2023/06/09
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
metapath <- "https://docs.google.com/spreadsheets/d/1gRL53qgRBApRHxznwTK_17I1cRlAL0GGf8nMikOJle0/edit#gid=0"
sheet_name <- "Sheet1"


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
sample_metadata_not_multiplexed_csf <- sample_metadata %>% 
  filter(Multiplexed == "No") %>%
  filter(Sample_Type == "CSF")


#==============================================================================#
# Read in Seurat data as a list ----
#==============================================================================#
seurat_list <- prep_seurat_list(
  sample_metadata_not_multiplexed_csf, 
  batch_ID = "Cell_Prefix", 
  cellRanger_path = "CellRanger_path", 
  cell_ID_prefix = "Cell_Prefix", 
  run_soupX = FALSE) # Run SoupX after filtering


#==============================================================================#
# Filter ----
#==============================================================================#
seurat_list_filtered <- filter_manual(seurat_list, 
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
                                                                 seurat_merge@meta.data$Cycle,  "_", 
                                                                 seurat_merge@meta.data$Day, "_",
                                                                 seurat_merge@meta.data$Sample_Type, "_",
                                                                 seurat_merge@meta.data$Batch)
table(seurat_merge@meta.data$UPN_Cycle_Day_Sample_Type_Batch)

seurat_merge <- merge(x = seurat_list_filtered[[1]], y = seurat_list_filtered[2:length(seurat_list_filtered)])
seurat_merge@meta.data$UPN_Cycle_Day_Sample_Type_Batch <- paste0(seurat_merge@meta.data$UPN, "_",
                                                                 seurat_merge@meta.data$Cycle,  "_", 
                                                                 seurat_merge@meta.data$Day, "_",
                                                                 seurat_merge@meta.data$Sample_Type, "_",
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
  cycle_id <- str_split(i, "_")[[1]][3]
  day_id <- str_split(i, "_")[[1]][4]
  batch_id <- str_split(i, "_")[[1]][5]
  sample_metadata_not_multiplexed_csf_i_CR_path <- sample_metadata_not_multiplexed_csf %>%
    filter(UPN == upn_id) %>% filter(Cycle == cycle_id) %>% 
    filter(Day == day_id) %>% filter(Batch == batch_id) %>% 
    pull(CellRanger_path)
  
  d10x_tod <- Read10X(paste0(sample_metadata_not_multiplexed_csf_i_CR_path, "/outs/raw_feature_bc_matrix/"))
  
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

# Add metadata
batch_ID = "FID_GEXFB"
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
    sample_seurat_list[[i]] <- NormalizeData(sample_seurat_list[[i]], verbose = FALSE)
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
                                                                                               seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$Cycle,  "_", 
                                                                                               seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$Day, "_",
                                                                                               seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$Sample_Type, "_",
                                                                                               seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$Batch)
  
  seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$Cycle_Day <- paste0("Cycle", seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$Cycle,  "_Day", 
                                                                         seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$Day)
}



#------------------------------------------------------------------------------#
## Integrate CSF ----
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


#==============================================================================#
# T cell integration ----
#==============================================================================#
# Subset potential T cells (T/NK) cells and re-integrate
T_NK_clusters <- c("0", "2", "3", "4", "5", "7", "8", "10", "11", "13")

# we are using resolution 0.2 
Idents(integrated_obj) <- "integrated_sct_snn_res.0.2"
# get number of cells in each cluster
table(Idents(integrated_obj))
# Get cell names from the clusters that I think are T cells
integrated_obj_T_NK <- subset(integrated_obj, idents = T_NK_clusters)
table(integrated_obj_T_NK@meta.data$integrated_sct_snn_res.0.2)

# Drop levels
integrated_obj_T_NK@meta.data[["integrated_sct_snn_res.0.2"]] <- droplevels(integrated_obj_T_NK@meta.data[["integrated_sct_snn_res.0.2"]])

# Split object by FID_GEXFB, since this is how they were split originally 
integrated_obj_T_NK_batches_list <- SplitObject(integrated_obj_T_NK, split.by = "FID_GEXFB")

# Re-run SCTransform
integrated_obj_T_NK_batches_list <- lapply(integrated_obj_T_NK_batches_list, function(xx){
  # Rerunning SCTransform
  DefaultAssay(xx) <- "SoupX_RNA"
  xx <- SCTransform(xx, 
                    method = "glmGamPoi",
                    vars.to.regress = c("S.Score", "G2M.Score"), 
                    vst.flavor = "v2",
                    verbose = T)
  
})

# Integrate
T_NK_integrated_obj <- basic_sct_rpca_integration(integrated_obj_T_NK_batches_list, npcs=50, k_weight=100, numfeatures = 2500)
#T_NK_integrated_obj <- readRDS("/scratch/aoill/projects/CAR-T/00_final/rds_files/csf_T_NK_integrated_20230329.rds")


## Remove NK cluster (4 an 11) and re-cluster ----
Idents(T_NK_integrated_obj) <- "integrated_sct_snn_res.0.2"
T_only_clusters <- c("0", "1", "2", "3", "5", "6", "7", "8", "9", "10", "12")
T_only_obj <- subset(T_NK_integrated_obj, idents = T_only_clusters)


# Re-cluster
npcs <- min(get_pcs(T_only_obj, reduction_name = "integrated_sct_pca"))
npcs
integrated_object_T_only <- RunUMAP(T_only_obj,
                                    reduction = "integrated_sct_pca",
                                    reduction.name = "integrated_sct_umap",
                                    dims = 1:npcs,
                                    return.model = TRUE)
integrated_object_T_only <- FindNeighbors(integrated_object_T_only,
                                          reduction = "integrated_sct_pca",
                                          dims = 1:npcs,
                                          graph.name = c("integrated_sct_nn",
                                                         "integrated_sct_snn"))
integrated_object_T_only <- FindClusters(integrated_object_T_only,
                                         resolution = c(0.1,0.2,0.3,0.5,0.8,1),
                                         graph.name = "integrated_sct_snn")


# Remove cluster 8 (cluster with LYZ expression)
Idents(integrated_object_T_only) <- "integrated_sct_snn_res.0.2"
table(Idents(integrated_object_T_only))
cells_to_keep <- c("0", "1", "2", "3", "4", "5", "6", "7", "9", "10")
integrated_object_T_only_filter <- subset(integrated_object_T_only, idents = cells_to_keep)
table(Idents(integrated_object_T_only_filter))

# Re-cluster
npcs <- min(get_pcs(integrated_object_T_only_filter, reduction_name = "integrated_sct_pca"))
npcs
obj_T_fltr_rclstr <- RunUMAP(integrated_object_T_only_filter,
                             reduction = "integrated_sct_pca",
                             reduction.name = "integrated_sct_umap",
                             dims = 1:npcs,
                             return.model = TRUE)
obj_T_fltr_rclstr <- FindNeighbors(obj_T_fltr_rclstr,
                                   reduction = "integrated_sct_pca",
                                   dims = 1:npcs,
                                   graph.name = c("integrated_sct_nn",
                                                  "integrated_sct_snn"))
obj_T_fltr_rclstr <- FindClusters(obj_T_fltr_rclstr,
                                  resolution = c(0.1,0.2,0.3,0.5,0.8,1),
                                  graph.name = "integrated_sct_snn")


# Get CD4 and CD8 clusters
DefaultAssay(obj_T_fltr_rclstr) <- "SCT"
feature_genes <- c("CD4", "CD8A",
                   "FOXP3", "IL2RA")
FeaturePlot(obj_T_fltr_rclstr, 
            features = feature_genes,
            reduction = "integrated_sct_umap", 
            ncol = 2) 

# sub-cluster 4
Idents(obj_T_fltr_rclstr) <- "integrated_sct_snn_res.0.2"
obj_T_fltr_rclstr_sbclstr_4 <- FindSubCluster(obj_T_fltr_rclstr, "4", resolution = 0.15, graph.name = "integrated_sct_snn")

# sub-cluster 2
Idents(obj_T_fltr_rclstr_sbclstr_4) <- "sub.cluster"
obj_T_fltr_rclstr_sbclstr_4_2 <- FindSubCluster(obj_T_fltr_rclstr_sbclstr_4, "2", resolution = 0.1, graph.name = "integrated_sct_snn")


## Add T ct labels 
# CD4: 0, 2_1, 4_0, 7, 8
# CD8: 1, 2_0, 3, 4_1, 5, 6, 9

# Add CD4/CD8 to the meta-data (col name ct_T)
Idents(obj_T_fltr_rclstr_sbclstr_4_2) <- "sub.cluster"
old_cluster_IDs <- obj_T_fltr_rclstr_sbclstr_4_2@meta.data$sub.cluster
levels(obj_T_fltr_rclstr_sbclstr_4_2)

new_cluster_ids <- c("CD4", "CD4", "CD4", "CD8", "CD8", "CD8", "CD4", "CD8", 
                     "CD4", "CD8", "CD8", "CD8")

names(new_cluster_ids) <- levels(obj_T_fltr_rclstr_sbclstr_4_2)
obj_T_fltr_rclstr_sbclstr_4_2 <- RenameIdents(obj_T_fltr_rclstr_sbclstr_4_2, new_cluster_ids)
DimPlot(obj_T_fltr_rclstr_sbclstr_4_2,  
        reduction = "integrated_sct_umap", label = TRUE, pt.size = 0.5) 

obj_T_fltr_rclstr_sbclstr_4_2@meta.data$ct_T <- Idents(obj_T_fltr_rclstr_sbclstr_4_2)


Idents(obj_T_fltr_rclstr_sbclstr_4_2) <- "sub.cluster"


#==============================================================================#
# Add metadata and other information ----
#==============================================================================#

#------------------------------------------------------------------------------#
## T cell object ----
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### Add response information T cell object ----
#------------------------------------------------------------------------------#
obj_T_fltr_rclstr_sbclstr_4_2@meta.data$UPN_Cycle <- paste(obj_T_fltr_rclstr_sbclstr_4_2@meta.data$UPN, obj_T_fltr_rclstr_sbclstr_4_2@meta.data$Cycle, sep = "_")
levels(as.factor(obj_T_fltr_rclstr_sbclstr_4_2@meta.data$UPN_Cycle))

obj_T_fltr_rclstr_sbclstr_4_2@meta.data$ResponseNon <- NA # or any other initialization value
obj_T_fltr_rclstr_sbclstr_4_2@meta.data$ResponseNon[obj_T_fltr_rclstr_sbclstr_4_2@meta.data$UPN_Cycle == "514_5"] <- "Non-Response" 
obj_T_fltr_rclstr_sbclstr_4_2@meta.data$ResponseNon[obj_T_fltr_rclstr_sbclstr_4_2@meta.data$UPN_Cycle == "514_8"] <- "Response"
obj_T_fltr_rclstr_sbclstr_4_2@meta.data$ResponseNon[obj_T_fltr_rclstr_sbclstr_4_2@meta.data$UPN_Cycle == "515_4"] <- "Response"
obj_T_fltr_rclstr_sbclstr_4_2@meta.data$ResponseNon[obj_T_fltr_rclstr_sbclstr_4_2@meta.data$UPN_Cycle == "574_8"] <- "Non-Response" 
obj_T_fltr_rclstr_sbclstr_4_2@meta.data$ResponseNon[obj_T_fltr_rclstr_sbclstr_4_2@meta.data$UPN_Cycle == "514_11"] <- "Non-Response" 
obj_T_fltr_rclstr_sbclstr_4_2@meta.data$ResponseNon[obj_T_fltr_rclstr_sbclstr_4_2@meta.data$UPN_Cycle == "515_8"] <- "Non-Response" 
obj_T_fltr_rclstr_sbclstr_4_2@meta.data$ResponseNon[obj_T_fltr_rclstr_sbclstr_4_2@meta.data$UPN_Cycle == "574_3"] <- "Non-Response" 

obj_T_fltr_rclstr_sbclstr_4_2@meta.data$Responses <- NA # or any other initialization value
obj_T_fltr_rclstr_sbclstr_4_2@meta.data$Responses[obj_T_fltr_rclstr_sbclstr_4_2@meta.data$UPN_Cycle == "514_5"] <- "Stable" 
obj_T_fltr_rclstr_sbclstr_4_2@meta.data$Responses[obj_T_fltr_rclstr_sbclstr_4_2@meta.data$UPN_Cycle == "514_8"] <- "Response"
obj_T_fltr_rclstr_sbclstr_4_2@meta.data$Responses[obj_T_fltr_rclstr_sbclstr_4_2@meta.data$UPN_Cycle == "515_4"] <- "Response"
obj_T_fltr_rclstr_sbclstr_4_2@meta.data$Responses[obj_T_fltr_rclstr_sbclstr_4_2@meta.data$UPN_Cycle == "574_8"] <- "Stable"
obj_T_fltr_rclstr_sbclstr_4_2@meta.data$Responses[obj_T_fltr_rclstr_sbclstr_4_2@meta.data$UPN_Cycle == "514_11"] <- "Progression"
obj_T_fltr_rclstr_sbclstr_4_2@meta.data$Responses[obj_T_fltr_rclstr_sbclstr_4_2@meta.data$UPN_Cycle == "515_8"] <- "Progression"
obj_T_fltr_rclstr_sbclstr_4_2@meta.data$Responses[obj_T_fltr_rclstr_sbclstr_4_2@meta.data$UPN_Cycle == "574_3"] <- "Pseudoprogression" 

table(obj_T_fltr_rclstr_sbclstr_4_2@meta.data$ResponseNon)
table(obj_T_fltr_rclstr_sbclstr_4_2@meta.data$Responses)


#------------------------------------------------------------------------------#
### Add T cell state annotations to T cell object ----
#------------------------------------------------------------------------------#
#### Subcluster ----
# sub-cluster 0
Idents(obj_T_fltr_rclstr_sbclstr_4_2) <- "sub.cluster"
obj_T_fltr_rclstr_sbclstr_4_2_0 <- FindSubCluster(obj_T_fltr_rclstr_sbclstr_4_2, "0", resolution = 0.35, graph.name = "integrated_sct_snn")

# sub-cluster 1
Idents(obj_T_fltr_rclstr_sbclstr_4_2_0) <- "sub.cluster"
obj_T_fltr_rclstr_sbclstr_4_2_0_1 <- FindSubCluster(obj_T_fltr_rclstr_sbclstr_4_2_0, "1", resolution = 0.3, graph.name = "integrated_sct_snn")

# sub-cluster 2_0
Idents(obj_T_fltr_rclstr_sbclstr_4_2_0_1) <- "sub.cluster"
obj_T_fltr_rclstr_sbclstr_4_2_0_1_0 <- FindSubCluster(obj_T_fltr_rclstr_sbclstr_4_2_0_1, "2_0", resolution = 0.1, graph.name = "integrated_sct_snn")

# sub-cluster 3
Idents(obj_T_fltr_rclstr_sbclstr_4_2_0_1_0) <- "sub.cluster"
obj_T_fltr_rclstr_sbclstr_4_2_0_1_0_3 <- FindSubCluster(obj_T_fltr_rclstr_sbclstr_4_2_0_1_0, "3", resolution = 0.15, graph.name = "integrated_sct_snn")

# sub-cluster 7
Idents(obj_T_fltr_rclstr_sbclstr_4_2_0_1_0_3) <- "sub.cluster"
obj_T_fltr_rclstr_sbclstr_fnl <- FindSubCluster(obj_T_fltr_rclstr_sbclstr_4_2_0_1_0_3, "7", resolution = 0.1, graph.name = "integrated_sct_snn")


#### Add cell state annotations ----
obj_T_fltr_rclstr_sbclstr_fnl@meta.data$ct_T_states <- obj_T_fltr_rclstr_sbclstr_fnl@meta.data$sub.cluster

obj_T_fltr_rclstr_sbclstr_fnl@meta.data$ct_T_states[
  obj_T_fltr_rclstr_sbclstr_fnl@meta.data$sub.cluster %in% c("9")] <- "CD8 Proliferating"
obj_T_fltr_rclstr_sbclstr_fnl@meta.data$ct_T_states[
  obj_T_fltr_rclstr_sbclstr_fnl@meta.data$sub.cluster %in% c("1_1")] <- "CD8 Naive"
obj_T_fltr_rclstr_sbclstr_fnl@meta.data$ct_T_states[
  obj_T_fltr_rclstr_sbclstr_fnl@meta.data$sub.cluster %in% c("1_0", "1_2", "1_3",
                                                             "2_0_1", "3_0", "6")] <- "CD8 Effector"
obj_T_fltr_rclstr_sbclstr_fnl@meta.data$ct_T_states[
  obj_T_fltr_rclstr_sbclstr_fnl@meta.data$sub.cluster %in% c("2_0_2", "4_1")] <- "CD8 Memory"
obj_T_fltr_rclstr_sbclstr_fnl@meta.data$ct_T_states[
  obj_T_fltr_rclstr_sbclstr_fnl@meta.data$sub.cluster %in% c("3_2", "5")] <- "CD8 Resident Memory-like"
obj_T_fltr_rclstr_sbclstr_fnl@meta.data$ct_T_states[
  obj_T_fltr_rclstr_sbclstr_fnl@meta.data$sub.cluster %in% c("2_0_0", "3_1", "7_0")] <- "CD8 Exhausted"
obj_T_fltr_rclstr_sbclstr_fnl@meta.data$ct_T_states[
  obj_T_fltr_rclstr_sbclstr_fnl@meta.data$sub.cluster %in% c("8")] <- "Treg"
obj_T_fltr_rclstr_sbclstr_fnl@meta.data$ct_T_states[
  obj_T_fltr_rclstr_sbclstr_fnl@meta.data$sub.cluster %in% c("0_0", "0_1", "0_2", "0_3", "4_0")] <- "CD4 Memory"
obj_T_fltr_rclstr_sbclstr_fnl@meta.data$ct_T_states[
  obj_T_fltr_rclstr_sbclstr_fnl@meta.data$sub.cluster %in% c("7_1")] <- "CD4 Activated"
obj_T_fltr_rclstr_sbclstr_fnl@meta.data$ct_T_states[
  obj_T_fltr_rclstr_sbclstr_fnl@meta.data$sub.cluster %in% c("2_1")] <- "CD4 Effector"

table(obj_T_fltr_rclstr_sbclstr_fnl@meta.data$ct_T_states)


#------------------------------------------------------------------------------#
### Add CART information to T cell object ----
#------------------------------------------------------------------------------#
DefaultAssay(obj_T_fltr_rclstr_sbclstr_fnl) <- "RNA"

obj_T_fltr_rclstr_sbclstr_fnl@meta.data$IL13OPCounts <- obj_T_fltr_rclstr_sbclstr_fnl@assays$RNA@counts["IL13OP",]
obj_T_fltr_rclstr_sbclstr_fnl@meta.data$CART <- ifelse(
  (obj_T_fltr_rclstr_sbclstr_fnl@meta.data$IL13OPCounts >= 3 & (obj_T_fltr_rclstr_sbclstr_fnl@meta.data$CT == "CD8+ T" | obj_T_fltr_rclstr_sbclstr_fnl@meta.data$CT == "CD4+ T" | obj_T_fltr_rclstr_sbclstr_fnl@meta.data$CT == "Treg" )), "Positive", "Negative")

table(obj_T_fltr_rclstr_sbclstr_fnl@meta.data$CART)


#------------------------------------------------------------------------------#
### Add TCR data to T cell object ----
#------------------------------------------------------------------------------#

#### Add TCR data determined from scRepertoire ----

##### Load the contig data ----

# TCR info for each sample

## CSF
F05038 <- read.csv("/tgen_labs/banovich/SingleCell/CellRanger/5_0_0/Ensemble_98/PipelineData/Projects/BCTCSF/CellRangerOuts/VDJ/BCTCSF_0097_1_CS_Whole_C4_X5TCR_F05038_HVFC5DSX3/outs/filtered_contig_annotations.csv")
F05039 <- read.csv("/tgen_labs/banovich/SingleCell/CellRanger/5_0_0/Ensemble_98/PipelineData/Projects/BCTCSF/CellRangerOuts/VDJ/BCTCSF_0098_1_CS_Whole_C3_X5TCR_F05039_HVFC5DSX3/outs/filtered_contig_annotations.csv")

F05040 <- read.csv("/tgen_labs/banovich/SingleCell/CellRanger/5_0_0/Ensemble_98/PipelineData/Projects/BCTCSF/CellRangerOuts/VDJ/BCTCSF_0099_1_CS_Whole_C2_X5TCR_F05040_HVFC5DSX3/outs/filtered_contig_annotations.csv")
F05042 <- read.csv("/tgen_labs/banovich/SingleCell/CellRanger/5_0_0/Ensemble_98/PipelineData/Projects/BCTCSF/CellRangerOuts/VDJ/BCTCSF_0097_1_CS_Whole_C6_X5TCR_F05042_HVFC5DSX3/outs/filtered_contig_annotations.csv")
F05043 <- read.csv("/tgen_labs/banovich/SingleCell/CellRanger/5_0_0/Ensemble_98/PipelineData/Projects/BCTCSF/CellRangerOuts/VDJ/BCTCSF_0098_1_CS_Whole_C5_X5TCR_F05043_HVFC5DSX3/outs/filtered_contig_annotations.csv")
F05044 <- read.csv("/tgen_labs/banovich/SingleCell/CellRanger/5_0_0/Ensemble_98/PipelineData/Projects/BCTCSF/CellRangerOuts/VDJ/BCTCSF_0099_1_CS_Whole_C5_X5TCR_F05044_HVFC5DSX3/outs/filtered_contig_annotations.csv")
F05045 <- read.csv("/tgen_labs/banovich/SingleCell/CellRanger/5_0_0/Ensemble_98/PipelineData/Projects/BCTCSF/CellRangerOuts/VDJ/BCTCSF_0097_1_CS_Whole_C7_X5TCR_F05045_HVFC5DSX3/outs/filtered_contig_annotations.csv")
F05046 <- read.csv("/tgen_labs/banovich/SingleCell/CellRanger/5_0_0/Ensemble_98/PipelineData/Projects/BCTCSF/CellRangerOuts/VDJ/BCTCSF_0097_1_CS_Whole_C4_X5TCR_F05046_HVFC5DSX3/outs/filtered_contig_annotations.csv")
F05047 <- read.csv("/tgen_labs/banovich/SingleCell/CellRanger/5_0_0/Ensemble_98/PipelineData/Projects/BCTCSF/CellRangerOuts/VDJ/BCTCSF_0098_1_CS_Whole_C3_X5TCR_F05047_HVFC5DSX3/outs/filtered_contig_annotations.csv")
F05048 <- read.csv("/tgen_labs/banovich/SingleCell/CellRanger/5_0_0/Ensemble_98/PipelineData/Projects/BCTCSF/CellRangerOuts/VDJ/BCTCSF_0099_1_CS_Whole_C2_X5TCR_F05048_HVFC5DSX3/outs/filtered_contig_annotations.csv")


##### Process the contig data ----
seurat_obj <- obj_T_fltr_rclstr_sbclstr_4_2

seurat_list <- SplitObject(seurat_obj, split.by = "FID_GEXFB")

# Strip the barcode extra labeling in seurat_list since the imported contig's 
# cell IDs are not formatted with a prefix (like Batch37_BARCODEID)
seurat_list_cells_renamed <- list()
for (i in names(seurat_list)){
  print(i)
  new_cell_IDs_tmp <- t(as.data.frame(str_split(rownames(seurat_list[[i]]@meta.data), "_")))
  new_cell_IDs <-  as.character(new_cell_IDs_tmp[,2])
  renamed_assay <- RenameCells(
    seurat_list[[i]],
    new.names = new_cell_IDs
  )
  seurat_list_cells_renamed <- c(seurat_list_cells_renamed, renamed_assay)
}
names(seurat_list_cells_renamed) <- names(seurat_list)

# Merge all contigs and contig lists into one contig list
contig_list_csf <- list(F05038, F05039, F05040, F05042, 
                        F05043, F05044, F05045, F05046, 
                        F05047, F05048)


names(contig_list_csf)
names(contig_list_csf) <- c("F05038_514_11_1_CSF", "F05039_515_4_2_CSF", "F05040_574_8_1_CSF", "F05042_514_8_1_CSF", 
                            "F05043_515_8_1_CSF", "F05044_574_3_1_CSF", "F05045_514_5_1_CSF", "F05046_514_11_1_CSF",
                            "F05047_515_4_2_CSF", "F05048_574_8_1_CSF")


##### Combine the contigs ----

# Will want to use the parameter removeNA set to true to remove any cell barcode 
# with an NA value in at least one of the chains
combined <- combineTCR(contig_list_csf, 
                       samples =  c("F05038_514_11_1_CSF", "F05039_515_4_2_CSF", "F05040_574_8_1_CSF", "F05042_514_8_1_CSF", 
                                    "F05043_515_8_1_CSF", "F05044_574_3_1_CSF", "F05045_514_5_1_CSF", "F05046_514_11_1_CSF",
                                    "F05047_515_4_2_CSF", "F05048_574_8_1_CSF"),
                       cells ="T-AB", filterMulti = FALSE, removeNA = TRUE
)



##### Add TCR info to the Seurat objects ----

### Add a new barcode column in the TCR list to match the Seurat objects
# Replace the barcode column in combined which matching barcode to seurat object 
# Path to google sheet with metadata for your samples
metapath <- "https://docs.google.com/spreadsheets/d/1fmEuOFXm893mfS2T1zpEsMrN7MpFc_Aym24ZODbw2bA/edit#gid=0"
sheet_name <- "Sheet1"
sample_metadata <- load_metadata(metadata_path = metapath, data_type = "googlesheet", sheet_name = sheet_name)

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
  j$orig_barcode <- sapply(strsplit(j$barcode, "_"), `[`, 6) # position in the string where the original barcode starts
  cell_id_prefix <- sample_metadata %>% filter(TCRseq_ID == l) %>% pull(Cell_Prefix) %>% unique()
  j$barcode <- paste0(cell_id_prefix, "_", j$orig_barcode)
  j$TCR_ID <- sapply(strsplit(j$sample, "_"), `[`, 1)
  j$UPN <- sapply(strsplit(j$sample, "_"), `[`, 2)
  j$Cycle <- sapply(strsplit(j$sample, "_"), `[`, 3)
  j$Day <- sapply(strsplit(j$sample, "_"), `[`, 4)
  j$Sample_Type <- sapply(strsplit(j$sample, "_"), `[`, 5)
  j$UPN_Cycle_Day_Sample_Type <- paste(j$UPN, j$Cycle, j$Day, j$Sample_Type, sep = "_")
  j$TCR_ID_UPN_Cycle_Day_Sample_Type <- paste(j$TCR_ID, j$UPN, j$Cycle, j$Day, j$Sample_Type, sep = "_")
  j
})
names(combined2) <- c("F05038_514_11_1_CSF", "F05039_515_4_2_CSF", "F05040_574_8_1_CSF", "F05042_514_8_1_CSF", 
                      "F05043_515_8_1_CSF", "F05044_574_3_1_CSF", "F05045_514_5_1_CSF", "F05046_514_11_1_CSF",
                      "F05047_515_4_2_CSF", "F05048_574_8_1_CSF")


#combined2 <- readRDS("/scratch/aoill/projects/CAR-T/00_final/rds_files/csf_T_TCR_clonotype_screp_20230405.rds")

obj_T_fltr_rclstr_sbclstr_fnl_TCR <- combineExpression(combined2, obj_T_fltr_rclstr_sbclstr_fnl,
                                                       group.by = "UPN_Cycle_Day_Sample_Type",
                                                       filterNA = FALSE,
                                                       proportion = FALSE,
                                                       cloneCall="strict",
                                                       cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

slot(obj_T_fltr_rclstr_sbclstr_fnl_TCR, 
     "meta.data")$cloneType <- factor(slot(obj_T_fltr_rclstr_sbclstr_fnl_TCR, "meta.data")$cloneType, 
                                      levels = c("Hyperexpanded (100 < X <= 500)", "Large (20 < X <= 100)", 
                                                 "Medium (5 < X <= 20)", "Small (1 < X <= 5)", 
                                                 "Single (0 < X <= 1)", NA))

## SAVE ##
saveRDS(combined2, "/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/CSF_TCR_clonotype_scRepertoire.rds")


#### Get normalized TCR frequencies and categories ----
# TCRs with a frequency of 1 will be categorized as non-expanded 
# Distribution of TCR frequencies normalized will only include TCRs with a 
# frequency > 1

# calculate Frequency_norm
tmp_meta <- obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data

# Add a T cell number count for each UPN patient and cycle
tmp_meta <- obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data %>%
  #add_count(UPN, name = "UPN_n")
  add_count(UPN_Cycle, name = "UPN_Cycle_n")

# Get normalized TCR frequency
tmp_meta$Frequency_norm <- tmp_meta$Frequency/tmp_meta$UPN_Cycle_n

obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$UPN_Cycle_n <- tmp_meta$UPN_Cycle_n
obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$Frequency_norm <- tmp_meta$Frequency_norm


tmp_meta_sub <- tmp_meta %>%
  dplyr::select(UPN, Sample_Type, Cycle, Cycle_Day, CART, CTstrict, Frequency,
                Frequency_norm)
rownames(tmp_meta) <- NULL
tmp_meta_sub <- na.omit(tmp_meta_sub)

# There could be duplicated rows in this metadata (clonotypes with freq >1 will 
# come up multiple times)
tmp_meta_sub_unique <- unique(tmp_meta_sub)

tmp_meta <- tmp_meta_sub_unique

# Remove any TCR with a frequency of 1
tmp_meta <- tmp_meta %>% dplyr::filter(Frequency > 1)

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
obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data <- obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data %>% 
  mutate(Frequency_norm_cat = case_when(Frequency == 1 ~ "Not Expanded",
                                        Frequency_norm < min(TCRs_top20$Frequency_norm) ~ "Not Expanded",
                                        (Frequency_norm >= min(TCRs_top20$Frequency_norm) & Frequency_norm < min(TCRs_top5$Frequency_norm)) ~ "Expanded (Top 20% - 5%)",
                                        (Frequency_norm >= min(TCRs_top5$Frequency_norm) & Frequency_norm < min(TCRs_top1$Frequency_norm)) ~ "More Expanded (Top 5% - 1%)",
                                        Frequency_norm >= min(TCRs_top1$Frequency_norm) ~ "Most Expanded (Top 1%)",
                                        
                                        TRUE ~ as.character(NA))) 
table(obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$Frequency_norm_cat)


# Get expanded vs not expanded categories
obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data <- obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data %>% 
  mutate(Frequency_norm_cat_2 = case_when(Frequency == 1 ~ "Not Expanded",
                                          Frequency_norm < min(TCRs_top20$Frequency_norm) ~ "Not Expanded",
                                          (Frequency_norm >= min(TCRs_top20$Frequency_norm) & Frequency_norm < min(TCRs_top5$Frequency_norm)) ~ "Expanded",
                                          (Frequency_norm >= min(TCRs_top5$Frequency_norm) & Frequency_norm < min(TCRs_top1$Frequency_norm)) ~ "Expanded",
                                          Frequency_norm >= min(TCRs_top1$Frequency_norm) ~ "Expanded",
                                          
                                          TRUE ~ as.character(NA))) 
table(obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$Frequency_norm_cat_2)


#------------------------------------------------------------------------------#
### Run scGSVA analysis ----
#------------------------------------------------------------------------------#
# To help with cell state annotations we used gene sets for various cell states
# and ran scGSVA to get enrichment scores
# Make sure to load Matrix 1.5.1
library(Matrix)


gs4_deauth()

canonical_markers  <- gs4_get("https://docs.google.com/spreadsheets/d/1S-NFiWaE9TqOUXyklvwevLpahCPHNNC8IJl9t1MPis4/edit#gid=1375531067")
sheet_names(canonical_markers)
canonical_markers <- read_sheet(canonical_markers, sheet = "Sheet3")
head(canonical_markers)
tail(canonical_markers)


canonical_markers_small <- canonical_markers[,which(colnames(canonical_markers) %in% c("RNA", "Gene_sets"))]
canonical_markers_small <- canonical_markers_small[complete.cases(canonical_markers_small),]
colnames(canonical_markers_small) <- c("GeneID", "Annot")
head(canonical_markers_small)

# Genes for module scoring
nkt_markers <- setdiff(canonical_markers[which(canonical_markers$Gene_sets=="NKT"),]$RNA, c(NA))
cd8_naive_markers <- setdiff(canonical_markers[which(canonical_markers$Gene_sets=="CD8_Naive"),]$RNA, c(NA))
cd4_naive_markers <- setdiff(canonical_markers[which(canonical_markers$Gene_sets=="CD4_Naive"),]$RNA, c(NA))
cd8_effector_markers <- setdiff(canonical_markers[which(canonical_markers$Gene_sets=="CD8_Effector"),]$RNA, c(NA))
cd4_effector_markers <- setdiff(canonical_markers[which(canonical_markers$Gene_sets=="CD4_Effector"),]$RNA, c(NA))
cd8_memory_markers <- setdiff(canonical_markers[which(canonical_markers$Gene_sets=="CD8_Memory"),]$RNA, c(NA))
cd4_memory_markers <- setdiff(canonical_markers[which(canonical_markers$Gene_sets=="CD4_Memory"),]$RNA, c(NA))
cd8_effector_memory_markers <- setdiff(canonical_markers[which(canonical_markers$Gene_sets=="CD8_Effector_Memory"),]$RNA, c(NA))
cd8_central_memory_markers <- setdiff(canonical_markers[which(canonical_markers$Gene_sets=="CD8_Central_Memory"),]$RNA, c(NA))
cd8_resident_memory_markers <- setdiff(canonical_markers[which(canonical_markers$Gene_sets=="CD8_Resident_Memory"),]$RNA, c(NA))
cd8_exhausted_markers <- setdiff(canonical_markers[which(canonical_markers$Gene_sets=="CD8_Exhausted"),]$RNA, c(NA))
activated_markers <- setdiff(canonical_markers[which(canonical_markers$Gene_sets=="Activated"),]$RNA, c(NA))
proliferating_markers <- setdiff(canonical_markers[which(canonical_markers$Gene_sets=="Proliferating"),]$RNA, c(NA))


hsko <- buildAnnot(species="human", keytype="SYMBOL", anntype="GO")
hsko@species
hsko@anntype <- "custom"
hsko@keytype

typeof(hsko@annot)
head(hsko@annot)
tmp <- hsko@annot


# Must have three columns: GeneID, Term or dummy, Annotation
canonical_markers_small[,c("Term")] <- canonical_markers_small$Annot
canonical_markers_small <- canonical_markers_small[,c("GeneID", "Term", "Annot")]

hsko@annot <- as.data.frame(canonical_markers_small)
head(hsko@annot)
tail(hsko@annot)
levels(as.factor(hsko@annot$Annot))
typeof(hsko@annot)

DefaultAssay(obj_T_fltr_rclstr_sbclstr_fnl_TCR) 
DefaultAssay(obj_T_fltr_rclstr_sbclstr_fnl_TCR) <- "SoupX_RNA"
res <- scgsva(obj_T_fltr_rclstr_sbclstr_fnl_TCR, hsko) 

res_df <- as.data.frame(res)

rownames(res_df)

identical(rownames(res_df), colnames(obj_T_fltr_rclstr_sbclstr_fnl_TCR))

obj_T_fltr_rclstr_sbclstr_fnl_TCR$NKT <- res_df$NKT
obj_T_fltr_rclstr_sbclstr_fnl_TCR$CD8_Naive <- res_df$CD8_Naive
obj_T_fltr_rclstr_sbclstr_fnl_TCR$CD4_Naive <- res_df$CD4_Naive
obj_T_fltr_rclstr_sbclstr_fnl_TCR$CD8_Effector <- res_df$CD8_Effector
obj_T_fltr_rclstr_sbclstr_fnl_TCR$CD4_Effector <- res_df$CD4_Effector
obj_T_fltr_rclstr_sbclstr_fnl_TCR$CD8_Memory <- res_df$CD8_Memory
obj_T_fltr_rclstr_sbclstr_fnl_TCR$CD4_Memory <- res_df$CD4_Memory
obj_T_fltr_rclstr_sbclstr_fnl_TCR$CD8_Effector_Memory <- res_df$CD8_Effector_Memory
obj_T_fltr_rclstr_sbclstr_fnl_TCR$CD8_Central_Memory <- res_df$CD8_Central_Memory
obj_T_fltr_rclstr_sbclstr_fnl_TCR$CD8_Resident_Memory <- res_df$CD8_Resident_Memory
obj_T_fltr_rclstr_sbclstr_fnl_TCR$CD8_Exhausted <- res_df$CD8_Exhausted
obj_T_fltr_rclstr_sbclstr_fnl_TCR$Activated <- res_df$Activated
obj_T_fltr_rclstr_sbclstr_fnl_TCR$Proliferating <- res_df$Proliferating


# Plot all on same axis 
maxpltval <- max((obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$NKT), (obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$CD8_Naive),
                 (obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$CD4_Naive), (obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$CD8_Effector), 
                 (obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$CD4_Effector), (obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$CD8_Memory), 
                 (obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$CD4_Memory), (obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$CD8_Effector_Memory),
                 (obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$CD8_Central_Memory), (obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$CD8_Resident_Memory),
                 (obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$CD8_Exhausted), (obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$Activated),
                 (obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$Proliferating))
minpltval <- min((obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$NKT), (obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$CD8_Naive),
                 (obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$CD4_Naive), (obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$CD8_Effector), 
                 (obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$CD4_Effector), (obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$CD8_Memory), 
                 (obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$CD4_Memory), (obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$CD8_Effector_Memory),
                 (obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$CD8_Central_Memory), (obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$CD8_Resident_Memory),
                 (obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$CD8_Exhausted), (obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$Activated),
                 (obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$Proliferating))

maxpltval <- 0.7
minpltval <- -0.4

FeaturePlot(obj_T_fltr_rclstr_sbclstr_fnl_TCR, features = c("NKT", "CD8_Naive", "CD8_Effector", "CD8_Memory", 
                                                            "CD8_Effector_Memory", "CD8_Central_Memory",
                                                            "CD8_Resident_Memory", "CD8_Exhausted",
                                                            "CD4_Naive", "CD4_Effector", "CD4_Memory",
                                                            "Activated", "Proliferating"), 
            reduction = "integrated_sct_umap", ncol = 4#, cols = c("white", "blue")
) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),
                         breaks=c(minpltval,maxpltval), limits=c(minpltval,maxpltval)) 



#------------------------------------------------------------------------------#
### Save object ----
#------------------------------------------------------------------------------#
saveRDS(obj_T_fltr_rclstr_sbclstr_fnl_TCR, 
        "/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/CSF_T_cell_obj_2023.rds")

#------------------------------------------------------------------------------#
## Full object ----
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### Add response information full object ----
#------------------------------------------------------------------------------#
integrated_obj@meta.data$UPN_Cycle <- paste(integrated_obj@meta.data$UPN, integrated_obj@meta.data$Cycle, sep = "_")
levels(as.factor(integrated_obj@meta.data$UPN_Cycle))

integrated_obj@meta.data$ResponseNon <- NA # or any other initialization value
integrated_obj@meta.data$ResponseNon[integrated_obj@meta.data$UPN_Cycle == "514_5"] <- "Non-Response" 
integrated_obj@meta.data$ResponseNon[integrated_obj@meta.data$UPN_Cycle == "514_8"] <- "Response"
integrated_obj@meta.data$ResponseNon[integrated_obj@meta.data$UPN_Cycle == "515_4"] <- "Response"
integrated_obj@meta.data$ResponseNon[integrated_obj@meta.data$UPN_Cycle == "574_8"] <- "Non-Response" 
integrated_obj@meta.data$ResponseNon[integrated_obj@meta.data$UPN_Cycle == "514_11"] <- "Non-Response" 
integrated_obj@meta.data$ResponseNon[integrated_obj@meta.data$UPN_Cycle == "515_8"] <- "Non-Response" 
integrated_obj@meta.data$ResponseNon[integrated_obj@meta.data$UPN_Cycle == "574_3"] <- "Non-Response" 

integrated_obj@meta.data$Responses <- NA # or any other initialization value
integrated_obj@meta.data$Responses[integrated_obj@meta.data$UPN_Cycle == "514_5"] <- "Stable" 
integrated_obj@meta.data$Responses[integrated_obj@meta.data$UPN_Cycle == "514_8"] <- "Response"
integrated_obj@meta.data$Responses[integrated_obj@meta.data$UPN_Cycle == "515_4"] <- "Response"
integrated_obj@meta.data$Responses[integrated_obj@meta.data$UPN_Cycle == "574_8"] <- "Stable"
integrated_obj@meta.data$Responses[integrated_obj@meta.data$UPN_Cycle == "514_11"] <- "Progression"
integrated_obj@meta.data$Responses[integrated_obj@meta.data$UPN_Cycle == "515_8"] <- "Progression"
integrated_obj@meta.data$Responses[integrated_obj@meta.data$UPN_Cycle == "574_3"] <- "Pseudoprogression" 

table(integrated_obj@meta.data$ResponseNon)
table(integrated_obj@meta.data$Responses)


#------------------------------------------------------------------------------#
### Add cell type annotations to full object ----
#------------------------------------------------------------------------------#

#### Sub-cluster ----
# Sub-cluster 1 and 12
Idents(integrated_obj) <- "integrated_sct_snn_res.0.2"
integrated_obj_subcluster_1 <- FindSubCluster(integrated_obj, "1", resolution = 0.07, graph.name = "integrated_sct_snn")

Idents(integrated_obj_subcluster_1) <- "sub.cluster"
integrated_obj_subcluster_1_12 <- FindSubCluster(integrated_obj_subcluster_1, "12", resolution = 0.05, graph.name = "integrated_sct_snn")


## Get cell IDs for non-T cells clusters that were in potential T cell object ----
# 1. Read in object with all
# T_NK_integrated_obj from above

# NK clusters: 4, 11
Idents(T_NK_integrated_obj) <- "integrated_sct_snn_res.0.2"
nk_clusters <- c("4", "11")
NK_integrated_obj <- subset(T_NK_integrated_obj, idents = nk_clusters)
NK_integrated_obj@meta.data$CT <- "NK"

# Pull out metadata
meta_NK_integrated_obj <- NK_integrated_obj@meta.data
colnames(meta_NK_integrated_obj)


# 2. Read in object with NK filtered out 
# T_only_obj from above
integrated_object_T_only <- T_only_obj

# LYZ cluster: 8
Idents(integrated_object_T_only) <- "integrated_sct_snn_res.0.2"
lyz_clusters <- c("8")
lyz_integrated_obj <- subset(integrated_object_T_only, idents = lyz_clusters)

lyz_integrated_obj@meta.data$CT <- "Monocytes"

# Pull out metadata
meta_lyz_integrated_obj <- lyz_integrated_obj@meta.data
colnames(meta_lyz_integrated_obj)


# 3. Read in object with NK and LYZ cluster filtered out 
s_object_subcluster_2_0_0_7 <- obj_T_fltr_rclstr_sbclstr_fnl_TCR

# CD4: "8", "0_3", "2_1", "0_2", "7_1", "0_0", "0_1", "4_0"
# CD8: "9", "1_1", "1_0", "1_2", "1_3", "3_0", "4_1", "6", "2_0_1", "3_2", "5", "2_0_2", "3_1", "2_0_0", "7_0"
# Potentially Tregs
s_object_subcluster_2_0_0_7@meta.data$CT <- s_object_subcluster_2_0_0_7@meta.data$sub.cluster

s_object_subcluster_2_0_0_7@meta.data$CT[
  s_object_subcluster_2_0_0_7@meta.data$sub.cluster %in% c("0_3", "2_1", 
                                                           "0_2", "7_1", "0_0", 
                                                           "0_1", "4_0")] <- "CD4+ T"
s_object_subcluster_2_0_0_7@meta.data$CT[
  s_object_subcluster_2_0_0_7@meta.data$sub.cluster %in% c("8")] <- "Treg"

s_object_subcluster_2_0_0_7@meta.data$CT[
  s_object_subcluster_2_0_0_7@meta.data$sub.cluster %in% c("9", "1_1", "1_0", 
                                                           "1_2", "1_3", "3_0", 
                                                           "4_1", "6", "2_0_1",
                                                           "3_2", "5", "2_0_2",
                                                           "3_1", "2_0_0", 
                                                           "7_0")] <- "CD8+ T"
table(s_object_subcluster_2_0_0_7@meta.data$CT)

meta_s_object_subcluster_2_0_0_7 <- s_object_subcluster_2_0_0_7@meta.data
colnames(meta_s_object_subcluster_2_0_0_7)

# rbind NK metadata, LYZ metadata, and T cell metadata - make sure all cols are 
# the same
cols_to_keep_tmp <- intersect(colnames(meta_NK_integrated_obj), colnames(meta_s_object_subcluster_2_0_0_7))
cols_to_keep <- intersect(cols_to_keep_tmp, colnames(meta_lyz_integrated_obj))

meta_NK_integrated_obj <- meta_NK_integrated_obj %>% dplyr::select(CT)
meta_lyz_integrated_obj <- meta_lyz_integrated_obj %>% dplyr::select(CT)
meta_s_object_subcluster_2_0_0_7 <- meta_s_object_subcluster_2_0_0_7 %>% dplyr::select(CT)

meta_all <- rbind(meta_NK_integrated_obj, meta_lyz_integrated_obj, meta_s_object_subcluster_2_0_0_7)

# Merge this metadata with integrated_obj_subcluster_1_12 metadata by row name which is by = 0
meta_integrated_obj <- integrated_obj_subcluster_1_12@meta.data

meta_integrated_obj_merge <- merge(meta_integrated_obj, meta_all, by = 0, all = T)

table(meta_integrated_obj_merge$CT)

####  Add in cell types ---
meta_integrated_obj_merge$CT[
  meta_integrated_obj_merge$sub.cluster %in% c("1_0", "9")] <- "Monocytes"
meta_integrated_obj_merge$CT[
  meta_integrated_obj_merge$sub.cluster %in% c("1_1")] <- "Macrophages"
meta_integrated_obj_merge$CT[
  meta_integrated_obj_merge$sub.cluster %in% c("14")] <- "pDC"
meta_integrated_obj_merge$CT[
  meta_integrated_obj_merge$sub.cluster %in% c("6")] <- "cDC2" 
meta_integrated_obj_merge$CT[
  meta_integrated_obj_merge$sub.cluster %in% c("16")] <- "cDC1"
meta_integrated_obj_merge$CT[
  meta_integrated_obj_merge$sub.cluster %in% c("15")] <- "mDC"
meta_integrated_obj_merge$CT[
  meta_integrated_obj_merge$sub.cluster %in% c("12_0")] <- "B cells"
meta_integrated_obj_merge$CT[
  meta_integrated_obj_merge$sub.cluster %in% c("12_1")] <- "Plasma"

table(meta_integrated_obj_merge$CT)

# full join 
tmp <- integrated_obj_subcluster_1_12@meta.data
tmp_join <- full_join(tmp, meta_integrated_obj_merge)
integrated_obj_subcluster_1_12@meta.data$CT <- tmp_join$CT

table(integrated_obj_subcluster_1_12@meta.data$CT)


## Remove monocytes that cluster with T cells and re-UMAP ----

# Get cell IDs of cells to keep
lyz_cells_to_remove <- rownames(lyz_integrated_obj@meta.data)
cells_to_keep <-  colnames(integrated_obj_subcluster_1_12)[!(colnames(integrated_obj_subcluster_1_12) %in% lyz_cells_to_remove)]

# Make cell name column
integrated_obj_subcluster_1_12@meta.data$CellName <- rownames(integrated_obj_subcluster_1_12@meta.data)
# Subset object
integrated_obj_subcluster_1_12_filter <- subset(integrated_obj_subcluster_1_12, 
                                                subset = (CellName %in% cells_to_keep) 
)
# Check to make sure the cells were removed
length(lyz_cells_to_remove)
table(integrated_obj_subcluster_1_12@meta.data$CT) - table(integrated_obj_subcluster_1_12_filter@meta.data$CT)


# Re-UMAP
npcs <- min(get_pcs(integrated_obj_subcluster_1_12_filter, reduction_name = "integrated_sct_pca"))
npcs
integrated_obj_subcluster_1_12_filter <- RunUMAP(integrated_obj_subcluster_1_12_filter,
                                                 reduction = "integrated_sct_pca",
                                                 reduction.name = "integrated_sct_umap",
                                                 dims = 1:npcs,
                                                 return.model = TRUE)


#------------------------------------------------------------------------------#
### Add CAR-T information to full object ----
#------------------------------------------------------------------------------#
DefaultAssay(integrated_obj_subcluster_1_12_filter) <- "RNA"

integrated_obj_subcluster_1_12_filter@meta.data$IL13OPCounts <- integrated_obj_subcluster_1_12_filter@assays$RNA@counts["IL13OP",]
integrated_obj_subcluster_1_12_filter@meta.data$CART <- ifelse(
  (integrated_obj_subcluster_1_12_filter@meta.data$IL13OPCounts >= 3 & (integrated_obj_subcluster_1_12_filter@meta.data$CT == "CD8+ T" | integrated_obj_subcluster_1_12_filter@meta.data$CT == "CD4+ T" | integrated_obj_subcluster_1_12_filter@meta.data$CT == "Treg")), "Positive", "Negative")
unique(integrated_obj_subcluster_1_12_filter@meta.data$CART)
table(integrated_obj_subcluster_1_12_filter@meta.data$CART)


#------------------------------------------------------------------------------#
### Save object ----
#------------------------------------------------------------------------------#
saveRDS(integrated_obj_subcluster_1_12_filter, 
        "/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/CSF_all_cells_obj_2023.rds")


#==============================================================================#
# Author : Angela M. Oill, aoill@tgen.org
# Date: 2023/06/09
# Project: Pediatric CAR-T 
# Description: Preprocessing and integration of PBMC samples
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
sample_metadata_multiplexed <- sample_metadata %>% filter(Multiplexed == "Yes")

sample_metadata_not_multiplexed <- sample_metadata %>% filter(Multiplexed == "No")


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

# Keep only PBMC samples
pbmc_seurat_list_demultiplexed_singlets <- list()
for (i in 1:length(seurat_list_demultiplexed_singlets)) {
  print(i)
  sub_i <- subset(seurat_list_demultiplexed_singlets[[i]], subset = Sample_Type == "PBMC")
  pbmc_seurat_list_demultiplexed_singlets <- c(pbmc_seurat_list_demultiplexed_singlets, sub_i)
}
names(pbmc_seurat_list_demultiplexed_singlets) <- c("Batch37_filtered", "Batch38_filtered", "Batch39_filtered", "Batch40_filtered")


pbmc_seurat_list <- pbmc_seurat_list_demultiplexed_singlets


#==============================================================================#
# Filter ----
#==============================================================================#
pbmc_seurat_list_filtered <- filter_manual(pbmc_seurat_list, 
                                           pt_mt = 10, 
                                           nFeature = 650,
                                           nCount = 1200)

seurat_list_separated <- list()
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
seurat_merge <- merge(x = pbmc_seurat_list[[1]], y = pbmc_seurat_list[2:length(pbmc_seurat_list)])
seurat_merge@meta.data$UPN_Cycle_Day_Sample_Type_Batch <- paste0(seurat_merge@meta.data$UPN, "_",
                                                                 seurat_merge@meta.data$Cycle,  "_", 
                                                                 seurat_merge@meta.data$Day, "_",
                                                                 seurat_merge@meta.data$Sample_Type, "_",
                                                                 seurat_merge@meta.data$Batch)
table(seurat_merge@meta.data$UPN_Cycle_Day_Sample_Type_Batch)

seurat_merge <- merge(x = pbmc_seurat_list_filtered[[1]], y = pbmc_seurat_list_filtered[2:length(pbmc_seurat_list_filtered)])
seurat_merge@meta.data$UPN_Cycle_Day_Sample_Type_Batch <- paste0(seurat_merge@meta.data$UPN, "_",
                                                                 seurat_merge@meta.data$Cycle,  "_", 
                                                                 seurat_merge@meta.data$Day, "_",
                                                                 seurat_merge@meta.data$Sample_Type, "_",
                                                                 seurat_merge@meta.data$Batch)
table(seurat_merge@meta.data$UPN_Cycle_Day_Sample_Type_Batch)


#==============================================================================#
# Run SoupX ----
#==============================================================================#
pbmc_samples <- names(seurat_list_separated)

#-----------------------------------------#
## Run SoupX on demultiplexed samples ----
#-----------------------------------------#
# Extract product samples from these above objects
sample_list <- pbmc_samples

pbmc_seurat_list_separated_SoupX <- sapply(sample_list,  function(i){
  print(i)
  # Read in count and droplet data
  d10x_toc <- Read10X(paste0("/scratch/aoill/projects/CAR-T/00_final/soupX/demultiplexed_", i))
  
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
    batch_path <- sample_metadata_multiplexed %>% filter(Batch == 40) %>% 
      pull(CellRanger_path) %>% unique()
    d10x_tod <- Read10X(paste0(batch_path, "/outs/raw_feature_bc_matrix/"))
  }
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
  
})
names(pbmc_seurat_list_separated_SoupX) <- sample_list

# Add metadata
batch_ID = "FID_GEXFB"
cellRanger_path = "CellRanger_path"
for (i in names(pbmc_seurat_list_separated_SoupX)) {
  upn_id <- str_split(i, "_")[[1]][1]
  sample_type_id <- str_split(i, "_")[[1]][2]
  cycle_type_id <- str_split(i, "_")[[1]][3]
  day_type_id <- str_split(i, "_")[[1]][4]
  batch_type_id <- str_split(i, "_")[[1]][5]
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
    sample_seurat_list[[i]] <- SCTransform(sample_seurat_list[[i]], 
                                           vst.flavor = "v2", verbose = F)
    DefaultAssay(sample_seurat_list[[i]]) <- "SCT"
    sample_seurat_list[[i]] <- CellCycleScoring(sample_seurat_list[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident = F) 
    DefaultAssay(sample_seurat_list[[i]]) <- rna_assay
  }
  
  return(sample_seurat_list)
}

pbmc_seurat_list_separated_SoupX <- add_cell_cycle_score_2(pbmc_seurat_list_separated_SoupX, rna_assay = "SoupX_RNA")


#-----------------------------------------#
## Re-normalize with cell cycle scores ----
#-----------------------------------------#
for (i in 1:length(pbmc_seurat_list_separated_SoupX)) {
  DefaultAssay(pbmc_seurat_list_separated_SoupX[[i]]) = "SoupX_RNA"
  pbmc_seurat_list_separated_SoupX[[i]] = SCTransform(pbmc_seurat_list_separated_SoupX[[i]], 
                                                      method = "glmGamPoi", 
                                                      vars.to.regress = c("S.Score", "G2M.Score"),
                                                      vst.flavor = "v2",
                                                      verbose = T)
}


#==============================================================================#
# Remove RB and MT genes ----
#==============================================================================#
# Remove ribosomal and mt genes
pbmc_seurat_list_separated_SoupX_noRBSMT <- list()
message("Removing ribosomal and mitochondrial genes from Seurat object")
for (i in 1:length(pbmc_seurat_list_separated_SoupX)) {
  RBMTgenes <- grep(pattern = "^RP[SL]|^MRP[SL]|^MT-", x = rownames(pbmc_seurat_list_separated_SoupX[[i]]@assays$RNA@data), value = TRUE, invert = TRUE)
  pbmc_seurat_list_separated_SoupX_noRBSMT[[i]] = subset(pbmc_seurat_list_separated_SoupX[[i]], features = RBMTgenes)
  
}
names(pbmc_seurat_list_separated_SoupX_noRBSMT) <- names(pbmc_seurat_list_separated_SoupX)


#==============================================================================#
# Integration ----
#==============================================================================#

#------------------------------------------------------------------------------#
## Add some more metadata ----
#------------------------------------------------------------------------------#
for (i in 1:length(pbmc_seurat_list_separated_SoupX_noRBSMT)) {
  pbmc_seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$UPN_Cycle_Day_Sample_Type_Batch <- paste0(pbmc_seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$UPN, "_",
                                                                                                    pbmc_seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$Cycle,  "_", 
                                                                                                    pbmc_seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$Day, "_",
                                                                                                    pbmc_seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$Sample_Type, "_",
                                                                                                    pbmc_seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$Batch)
  
  pbmc_seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$Cycle_Day <- paste0("Cycle", pbmc_seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$Cycle,  "_Day", 
                                                                              pbmc_seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$Day)
}



#------------------------------------------------------------------------------#
## Integrate PBMCs ----
#------------------------------------------------------------------------------#
pbmc_seurat_list <- pbmc_seurat_list_separated_SoupX_noRBSMT

pbmc_seurat_list <- lapply(pbmc_seurat_list, function(xx){
  # Rerunning SCTransform
  DefaultAssay(xx) <- "SoupX_RNA"
  xx <- SCTransform(xx, 
                    method = "glmGamPoi",
                    vars.to.regress = c("S.Score", "G2M.Score"), 
                    vst.flavor = "v2",
                    verbose = T)
  
})

pbmc_integrated <- basic_sct_rpca_integration(pbmc_seurat_list, npcs=50, k_weight=100, numfeatures = 2500)


#==============================================================================#
# T cell integration ----
#==============================================================================#
#------------------------------------------------------------------------------#
## Subset potential T cell clusters and re-integrate ----
#------------------------------------------------------------------------------#
# potential T cell clusters
potential_T_clusters <- c("1", "2", "3", "5", "6", "7", "13", "14")

# we are using resolution 0.2 
Idents(pbmc_integrated) <- "integrated_sct_snn_res.0.2"
# get number of cells in each cluster
table(Idents(pbmc_integrated))
# Get cell names from the clusters that I think are T cells
integrated_obj_potential_T_clusters <- subset(pbmc_integrated, idents = potential_T_clusters)
table(integrated_obj_potential_T_clusters@meta.data$integrated_sct_snn_res.0.2)

# Drop levels
integrated_obj_potential_T_clusters@meta.data[["integrated_sct_snn_res.0.2"]] <- droplevels(integrated_obj_potential_T_clusters@meta.data[["integrated_sct_snn_res.0.2"]])

# Split object by FID_GEXFB, since this is how they were split originally 
integrated_obj_potential_T_clusters_batches_list <- SplitObject(integrated_obj_potential_T_clusters, split.by = "FID_GEXFB")

# Re-run SCTransform
integrated_obj_potential_T_clusters_batches_list <- lapply(integrated_obj_potential_T_clusters_batches_list, function(xx){
  # Rerunning SCTransform
  DefaultAssay(xx) <- "SoupX_RNA"
  xx <- SCTransform(xx, 
                    method = "glmGamPoi",
                    vars.to.regress = c("S.Score", "G2M.Score"), 
                    vst.flavor = "v2",
                    verbose = T)
  
})


# Integrate
potential_T_clusters_integrated_obj <- basic_sct_rpca_integration(integrated_obj_potential_T_clusters_batches_list, npcs=50, k_weight=100, numfeatures = 2500)

#------------------------------------------------------------------------------#
### Subset T cells and re-cluster ----
#------------------------------------------------------------------------------#
Idents(potential_T_clusters_integrated_obj) <- "integrated_sct_snn_res.0.5"
T_only_clusters <- c("2", "3", "4", "6", "7", "8", "9", "12", "13")
T_only_obj <- subset(potential_T_clusters_integrated_obj, idents = T_only_clusters)


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


#==============================================================================#
# Add metadata ----
#==============================================================================#

#------------------------------------------------------------------------------#
## T cell object ----
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### Add response information T cell object ----
#------------------------------------------------------------------------------#
integrated_object_T_only@meta.data$UPN_Cycle <- paste(integrated_object_T_only@meta.data$UPN, integrated_object_T_only@meta.data$Cycle, sep = "_")
levels(as.factor(integrated_object_T_only@meta.data$UPN_Cycle))

integrated_object_T_only@meta.data$ResponseNon <- NA # or any other initialization value
integrated_object_T_only@meta.data$ResponseNon[integrated_object_T_only@meta.data$UPN_Cycle == "514_5"] <- "Non-Response" 
integrated_object_T_only@meta.data$ResponseNon[integrated_object_T_only@meta.data$UPN_Cycle == "514_8"] <- "Response"
integrated_object_T_only@meta.data$ResponseNon[integrated_object_T_only@meta.data$UPN_Cycle == "515_4"] <- "Response"
integrated_object_T_only@meta.data$ResponseNon[integrated_object_T_only@meta.data$UPN_Cycle == "574_8"] <- "Non-Response" 
integrated_object_T_only@meta.data$ResponseNon[integrated_object_T_only@meta.data$UPN_Cycle == "514_11"] <- "Non-Response" 
integrated_object_T_only@meta.data$ResponseNon[integrated_object_T_only@meta.data$UPN_Cycle == "515_8"] <- "Non-Response" 
integrated_object_T_only@meta.data$ResponseNon[integrated_object_T_only@meta.data$UPN_Cycle == "574_3"] <- "Non-Response" 
integrated_object_T_only@meta.data$ResponseNon[integrated_object_T_only@meta.data$UPN_Cycle == "625_1"] <- "Non-Response" 
integrated_object_T_only@meta.data$ResponseNon[integrated_object_T_only@meta.data$UPN_Cycle == "625_3"] <- "Non-Response" 
integrated_object_T_only@meta.data$ResponseNon[integrated_object_T_only@meta.data$UPN_Cycle == "625_4"] <- "Non-Response" 
integrated_object_T_only@meta.data$ResponseNon[integrated_object_T_only@meta.data$UPN_Cycle == "626_1"] <- "Non-Response" 
integrated_object_T_only@meta.data$ResponseNon[integrated_object_T_only@meta.data$UPN_Cycle == "626_4"] <- "Response" 
integrated_object_T_only@meta.data$ResponseNon[integrated_object_T_only@meta.data$UPN_Cycle == "626_7"] <- "Non-Response" 

integrated_object_T_only@meta.data$Responses <- NA # or any other initialization value
integrated_object_T_only@meta.data$Responses[integrated_object_T_only@meta.data$UPN_Cycle == "514_5"] <- "Stable" 
integrated_object_T_only@meta.data$Responses[integrated_object_T_only@meta.data$UPN_Cycle == "514_8"] <- "Response"
integrated_object_T_only@meta.data$Responses[integrated_object_T_only@meta.data$UPN_Cycle == "515_4"] <- "Response"
integrated_object_T_only@meta.data$Responses[integrated_object_T_only@meta.data$UPN_Cycle == "574_8"] <- "Stable"
integrated_object_T_only@meta.data$Responses[integrated_object_T_only@meta.data$UPN_Cycle == "514_11"] <- "Progression"
integrated_object_T_only@meta.data$Responses[integrated_object_T_only@meta.data$UPN_Cycle == "515_8"] <- "Progression"
integrated_object_T_only@meta.data$Responses[integrated_object_T_only@meta.data$UPN_Cycle == "574_3"] <- "Pseudoprogression" 
integrated_object_T_only@meta.data$Responses[integrated_object_T_only@meta.data$UPN_Cycle == "625_1"] <- "Baseline" 
integrated_object_T_only@meta.data$Responses[integrated_object_T_only@meta.data$UPN_Cycle == "625_3"] <- "Progression" 
integrated_object_T_only@meta.data$Responses[integrated_object_T_only@meta.data$UPN_Cycle == "625_4"] <- "Progression" 
integrated_object_T_only@meta.data$Responses[integrated_object_T_only@meta.data$UPN_Cycle == "626_1"] <- "Baseline" 
integrated_object_T_only@meta.data$Responses[integrated_object_T_only@meta.data$UPN_Cycle == "626_4"] <- "Response" 
integrated_object_T_only@meta.data$Responses[integrated_object_T_only@meta.data$UPN_Cycle == "626_7"] <- "Progression" 

table(integrated_object_T_only@meta.data$ResponseNon)
table(integrated_object_T_only@meta.data$Responses)


#------------------------------------------------------------------------------#
### Add lymphodepletion information T cell object ----
#------------------------------------------------------------------------------#
integrated_object_T_only@meta.data$Lymphodepletion <- NA # or any other initialization value
integrated_object_T_only@meta.data$Lymphodepletion[integrated_object_T_only@meta.data$UPN == "514"] <- "Non-Lymphodepleted"
integrated_object_T_only@meta.data$Lymphodepletion[integrated_object_T_only@meta.data$UPN == "515"] <- "Non-Lymphodepleted"
integrated_object_T_only@meta.data$Lymphodepletion[integrated_object_T_only@meta.data$UPN == "574"] <- "Non-Lymphodepleted" 
integrated_object_T_only@meta.data$Lymphodepletion[integrated_object_T_only@meta.data$UPN == "625"] <- "Lymphodepleted" 
integrated_object_T_only@meta.data$Lymphodepletion[integrated_object_T_only@meta.data$UPN == "626"] <- "Lymphodepleted" 

table(integrated_object_T_only@meta.data$Lymphodepletion)


#------------------------------------------------------------------------------#
### Add T cell state annotations to T cell object ----
#------------------------------------------------------------------------------#
#### Sub-cluster ----

# Sub-cluster 1
Idents(integrated_object_T_only) <- "integrated_sct_snn_res.0.2"
integrated_object_T_only_subcluster_1 <- FindSubCluster(integrated_object_T_only, "1", resolution = 0.2, graph.name = "integrated_sct_snn")

# Sub-cluster 2
Idents(integrated_object_T_only_subcluster_1) <- "sub.cluster"
integrated_object_T_only_subcluster_1_2 <- FindSubCluster(integrated_object_T_only_subcluster_1, "2", resolution = 0.2, graph.name = "integrated_sct_snn")


# Sub-cluster 3
Idents(integrated_object_T_only_subcluster_1_2) <- "sub.cluster"
integrated_object_T_only_subcluster_1_2_3 <- FindSubCluster(integrated_object_T_only_subcluster_1_2, "3", resolution = 0.1, graph.name = "integrated_sct_snn")

# Sub-cluster 0
# Cluster 0_2 has plasma marker expression (CD79A and MZB1). So remove this
# cluster and at the end re-cluster
Idents(integrated_object_T_only_subcluster_1_2_3) <- "sub.cluster"
integrated_object_T_only_subcluster_1_2_3_0 <- FindSubCluster(integrated_object_T_only_subcluster_1_2_3, "0", resolution = 0.1, graph.name = "integrated_sct_snn")

# Remove plasma cluster (0_2)
Idents(integrated_object_T_only_subcluster_1_2_3_0) <- "sub.cluster"
integrated_object_T_only_subcluster_1_2_3_0_sub <- subset(integrated_object_T_only_subcluster_1_2_3_0, subset = sub.cluster != "0_2")

# Re-UMAP
npcs <- min(get_pcs(integrated_object_T_only_subcluster_1_2_3_0_sub, reduction_name = "integrated_sct_pca"))
npcs

obj_T_fltr_rclstr_sbclstr_fnl <- RunUMAP(integrated_object_T_only_subcluster_1_2_3_0_sub,
                                              reduction = "integrated_sct_pca",
                                              reduction.name = "integrated_sct_umap",
                                              dims = 1:npcs,
                                              return.model = TRUE)


#DimPlot(obj_T_fltr_rclstr_sbclstr_fnl, group.by = "sub.cluster", reduction = "integrated_sct_umap", label = T) 


#### Add cell state annotations ----
obj_T_fltr_rclstr_sbclstr_fnl@meta.data$ct_T_states <- obj_T_fltr_rclstr_sbclstr_fnl@meta.data$sub.cluster

obj_T_fltr_rclstr_sbclstr_fnl@meta.data$ct_T_states[
  obj_T_fltr_rclstr_sbclstr_fnl@meta.data$sub.cluster %in% c("1_0", "3_0", "4")] <- "CD8 Effector"
obj_T_fltr_rclstr_sbclstr_fnl@meta.data$ct_T_states[
  obj_T_fltr_rclstr_sbclstr_fnl@meta.data$sub.cluster %in% c("1_1", "5")] <- "CD8 Effector Memory"
obj_T_fltr_rclstr_sbclstr_fnl@meta.data$ct_T_states[
  obj_T_fltr_rclstr_sbclstr_fnl@meta.data$sub.cluster %in% c("0_0", "0_1")] <- "CD4 Memory"
obj_T_fltr_rclstr_sbclstr_fnl@meta.data$ct_T_states[
  obj_T_fltr_rclstr_sbclstr_fnl@meta.data$sub.cluster %in% c("2_1")] <- "CD4 Naive"
obj_T_fltr_rclstr_sbclstr_fnl@meta.data$ct_T_states[
  obj_T_fltr_rclstr_sbclstr_fnl@meta.data$sub.cluster %in% c("2_0", "3_1")] <- "CD8 Naive"
obj_T_fltr_rclstr_sbclstr_fnl@meta.data$ct_T_states[
  obj_T_fltr_rclstr_sbclstr_fnl@meta.data$sub.cluster %in% c("6")] <- "CD8 Activated"
obj_T_fltr_rclstr_sbclstr_fnl@meta.data$ct_T_states[
  obj_T_fltr_rclstr_sbclstr_fnl@meta.data$sub.cluster %in% c("7")] <- "CD8 Proliferating"

#DimPlot(obj_T_fltr_rclstr_sbclstr_fnl, group.by = "ct_T_states", label = T)


#------------------------------------------------------------------------------#
### Add CAR-T info to T cell object ----
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


##### Process the contig data ----

# Split by batch
seurat_obj <- obj_T_fltr_rclstr_sbclstr_fnl

seurat_list <- SplitObject(seurat_obj, split.by = "Batch")

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

# Get the contig list for the multiplexed batches
contig_list_37 <- createHTOContigList(F05037, seurat_list_cells_renamed[["37"]], group.by = c("TCRseq_ID", "UPN", "Cycle", "Day", "Sample_Type"))
contig_list_38 <- createHTOContigList(F05041, seurat_list_cells_renamed[["38"]], group.by = c("TCRseq_ID", "UPN", "Cycle", "Day", "Sample_Type"))
contig_list_39 <- createHTOContigList(F05625, seurat_list_cells_renamed[["39"]], group.by = c("TCRseq_ID", "UPN", "Cycle", "Day", "Sample_Type"))
contig_list_40 <- createHTOContigList(F05626, seurat_list_cells_renamed[["40"]], group.by = c("TCRseq_ID", "UPN", "Cycle", "Day", "Sample_Type"))

# Merge all contigs and contig lists into one contig list
contig_list_pbmc <- append(contig_list_37, contig_list_38)
contig_list_pbmc <- append(contig_list_pbmc, contig_list_39)
contig_list_pbmc <- append(contig_list_pbmc, contig_list_40)

names(contig_list_pbmc)
# [1] "F05037.514.8.1.PBMC"  "F05037.514.5.1.PBMC"  "F05037.574.3.1.PBMC"  "F05037.515.8.1.PBMC"  "F05041.514.11.1.PBMC" "F05041.515. 4.2.PBMC"
# [7] "F05041.574. 8.1.PBMC" "F05625.626.1.0.PBMC"  "F05625.625.3.1.PBMC"  "F05625.626.4.1.PBMC"  "F05626.626.7.1.PBMC"  "F05626.625.4.1.PBMC" 
# [13] "F05626.625.1.0.PBMC" 

names(contig_list_pbmc) <- c("F05037_514_8_1_PBMC", "F05037_514_5_1_PBMC", 
                             "F05037_574_3_1_PBMC", "F05037_515_8_1_PBMC", 
                             "F05041_514_11_1_PBMC", "F05041_515_4_2_PBMC",
                             "F05041_574_8_1_PBMC", "F05625_626_1_0_PBMC", 
                             "F05625_625_3_1_PBMC", "F05625_626_4_1_PBMC", 
                             "F05626_626_7_1_PBMC", "F05626_625_4_1_PBMC", 
                             "F05626_625_1_0_PBMC")

##### Combine the contigs ----

# Will want to use the parameter removeNA set to true to remove any cell barcode 
# with an NA value in at least one of the chains
combined <- combineTCR(contig_list_pbmc, 
                       samples =  c("F05037_514_8_1_PBMC", "F05037_514_5_1_PBMC", 
                                    "F05037_574_3_1_PBMC", "F05037_515_8_1_PBMC", 
                                    "F05041_514_11_1_PBMC", "F05041_515_4_2_PBMC",
                                    "F05041_574_8_1_PBMC", "F05625_626_1_0_PBMC", 
                                    "F05625_625_3_1_PBMC", "F05625_626_4_1_PBMC", 
                                    "F05626_626_7_1_PBMC", "F05626_625_4_1_PBMC", 
                                    "F05626_625_1_0_PBMC"),
                       cells ="T-AB", filterMulti = FALSE, removeNA = TRUE
)


##### Add TCR info to the Seurat objects ----

### Add a new barcode column in the TCR list to match the Seurat objects
# I need to replace the barcode column in combined which matching barcode to 
# seurat object 
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


names(combined2) <- c("F05037_514_8_1_PBMC", "F05037_514_5_1_PBMC", 
                      "F05037_574_3_1_PBMC", "F05037_515_8_1_PBMC", 
                      "F05041_514_11_1_PBMC", "F05041_515_4_2_PBMC",
                      "F05041_574_8_1_PBMC", "F05625_626_1_0_PBMC", 
                      "F05625_625_3_1_PBMC", "F05625_626_4_1_PBMC", 
                      "F05626_626_7_1_PBMC", "F05626_625_4_1_PBMC", 
                      "F05626_625_1_0_PBMC")



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
saveRDS(combined2, "/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/PBMC_TCR_clonotype_scReperroire.rds")


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
# Re-do this analysis with the set of genes from Colt
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

FeaturePlot(obj_T_fltr_rclstr_sbclstr_fnl_TCR, features = c(#"NKT", 
  "CD8_Naive", "CD8_Effector", "CD8_Memory", 
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
        "/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/PBMC_T_cell_obj_2023.rds")


#------------------------------------------------------------------------------#
## Full object ----
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### Add response information full object ----
#------------------------------------------------------------------------------#
pbmc_integrated@meta.data$UPN_Cycle <- paste(pbmc_integrated@meta.data$UPN, pbmc_integrated@meta.data$Cycle, sep = "_")
levels(as.factor(pbmc_integrated@meta.data$UPN_Cycle))

pbmc_integrated@meta.data$ResponseNon <- NA # or any other initialization value
pbmc_integrated@meta.data$ResponseNon[pbmc_integrated@meta.data$UPN_Cycle == "514_5"] <- "Non-Response" # Stable, I originally had this as response but after talking to Leo we decided to change this to stable/non-response
pbmc_integrated@meta.data$ResponseNon[pbmc_integrated@meta.data$UPN_Cycle == "514_8"] <- "Response"
pbmc_integrated@meta.data$ResponseNon[pbmc_integrated@meta.data$UPN_Cycle == "515_4"] <- "Response"
pbmc_integrated@meta.data$ResponseNon[pbmc_integrated@meta.data$UPN_Cycle == "574_8"] <- "Non-Response" # Stable
pbmc_integrated@meta.data$ResponseNon[pbmc_integrated@meta.data$UPN_Cycle == "514_11"] <- "Non-Response" # Progression
pbmc_integrated@meta.data$ResponseNon[pbmc_integrated@meta.data$UPN_Cycle == "515_8"] <- "Non-Response" # Progression
pbmc_integrated@meta.data$ResponseNon[pbmc_integrated@meta.data$UPN_Cycle == "574_3"] <- "Non-Response" # Pseudoprogression
pbmc_integrated@meta.data$ResponseNon[pbmc_integrated@meta.data$UPN_Cycle == "625_1"] <- "Non-Response" 
pbmc_integrated@meta.data$ResponseNon[pbmc_integrated@meta.data$UPN_Cycle == "625_3"] <- "Non-Response" 
pbmc_integrated@meta.data$ResponseNon[pbmc_integrated@meta.data$UPN_Cycle == "625_4"] <- "Non-Response" 
pbmc_integrated@meta.data$ResponseNon[pbmc_integrated@meta.data$UPN_Cycle == "626_1"] <- "Non-Response" 
pbmc_integrated@meta.data$ResponseNon[pbmc_integrated@meta.data$UPN_Cycle == "626_4"] <- "Response" 
pbmc_integrated@meta.data$ResponseNon[pbmc_integrated@meta.data$UPN_Cycle == "626_7"] <- "Non-Response" 

pbmc_integrated@meta.data$Responses <- NA # or any other initialization value
pbmc_integrated@meta.data$Responses[pbmc_integrated@meta.data$UPN_Cycle == "514_5"] <- "Stable" 
pbmc_integrated@meta.data$Responses[pbmc_integrated@meta.data$UPN_Cycle == "514_8"] <- "Response"
pbmc_integrated@meta.data$Responses[pbmc_integrated@meta.data$UPN_Cycle == "515_4"] <- "Response"
pbmc_integrated@meta.data$Responses[pbmc_integrated@meta.data$UPN_Cycle == "574_8"] <- "Stable"
pbmc_integrated@meta.data$Responses[pbmc_integrated@meta.data$UPN_Cycle == "514_11"] <- "Progression"
pbmc_integrated@meta.data$Responses[pbmc_integrated@meta.data$UPN_Cycle == "515_8"] <- "Progression"
pbmc_integrated@meta.data$Responses[pbmc_integrated@meta.data$UPN_Cycle == "574_3"] <- "Pseudoprogression" # Do we change to stable?
pbmc_integrated@meta.data$Responses[pbmc_integrated@meta.data$UPN_Cycle == "625_1"] <- "Baseline" 
pbmc_integrated@meta.data$Responses[pbmc_integrated@meta.data$UPN_Cycle == "625_3"] <- "Progression" 
pbmc_integrated@meta.data$Responses[pbmc_integrated@meta.data$UPN_Cycle == "625_4"] <- "Progression" 
pbmc_integrated@meta.data$Responses[pbmc_integrated@meta.data$UPN_Cycle == "626_1"] <- "Baseline" 
pbmc_integrated@meta.data$Responses[pbmc_integrated@meta.data$UPN_Cycle == "626_4"] <- "Response" 
pbmc_integrated@meta.data$Responses[pbmc_integrated@meta.data$UPN_Cycle == "626_7"] <- "Progression" 

table(pbmc_integrated@meta.data$ResponseNon)
table(pbmc_integrated@meta.data$Responses)


#------------------------------------------------------------------------------#
### Add lymphodepletion information full object ----
#------------------------------------------------------------------------------#
pbmc_integrated@meta.data$Lymphodepletion <- NA # or any other initialization value
pbmc_integrated@meta.data$Lymphodepletion[pbmc_integrated@meta.data$UPN == "514"] <- "Non-Lymphodepleted"
pbmc_integrated@meta.data$Lymphodepletion[pbmc_integrated@meta.data$UPN == "515"] <- "Non-Lymphodepleted"
pbmc_integrated@meta.data$Lymphodepletion[pbmc_integrated@meta.data$UPN == "574"] <- "Non-Lymphodepleted" 
pbmc_integrated@meta.data$Lymphodepletion[pbmc_integrated@meta.data$UPN == "625"] <- "Lymphodepleted" 
pbmc_integrated@meta.data$Lymphodepletion[pbmc_integrated@meta.data$UPN == "626"] <- "Lymphodepleted" 

table(pbmc_integrated@meta.data$Lymphodepletion)


#------------------------------------------------------------------------------#
### Add cell type annotations to full object ----
#------------------------------------------------------------------------------#
#### Sub-cluster ----
# Sub-cluster 4
Idents(pbmc_integrated) <- "integrated_sct_snn_res.0.2"
pbmc_integrated_sbclstr_4 <- FindSubCluster(pbmc_integrated, "4", resolution = 0.05, graph.name = "integrated_sct_snn")
DimPlot(pbmc_integrated_sbclstr_4, group.by = "sub.cluster", label = T)

# Sub-cluster 12
Idents(pbmc_integrated_sbclstr_4) <- "sub.cluster"
pbmc_integrated_sbclstr_4_12 <- FindSubCluster(pbmc_integrated_sbclstr_4, "12", resolution = 0.01, graph.name = "integrated_sct_snn")
DimPlot(pbmc_integrated_sbclstr_4_12, group.by = "sub.cluster", label = T)


# Add a label in meta called CT
# Pull out the following clusters from potential T object: 5, 0, 1, 10, 11, 14
Idents(potential_T_clusters_integrated_obj) <- "integrated_sct_snn_res.0.5"
nk_clusters <- c("1", "10", "11")
NK_integrated_obj <- subset(potential_T_clusters_integrated_obj, idents = nk_clusters)
NK_integrated_obj@meta.data$CT <- "NK"

mono_clusters <- c("5")
mono_extra_integrated_obj <- subset(potential_T_clusters_integrated_obj, idents = mono_clusters)
mono_extra_integrated_obj@meta.data$CT <- "CD14+ Monocytes"

proliferating_clusters <- c("0")
proliferating_extra_integrated_obj <- subset(potential_T_clusters_integrated_obj, idents = proliferating_clusters)
proliferating_extra_integrated_obj@meta.data$CT <- "Proliferating"

plasma_clusters <- c("14")
plasma_extra_integrated_obj <- subset(potential_T_clusters_integrated_obj, idents = plasma_clusters)
plasma_extra_integrated_obj@meta.data$CT <- "Plasma"

# CTs from final T object:
# CD4+ T: 0_0, 0_1, 2_1
# CD8+ T: 1_0, 1_1, 2_0, 4, 5, 6, 7, 3_0, 3_1
integrated_object_T_only_subcluster_1_2_3_0_sub@meta.data$CT <- integrated_object_T_only_subcluster_1_2_3_0_sub@meta.data$sub.cluster

integrated_object_T_only_subcluster_1_2_3_0_sub@meta.data$CT[
  integrated_object_T_only_subcluster_1_2_3_0_sub@meta.data$sub.cluster %in% c("0_0", "0_1", "2_1")] <- "CD4+ T"
integrated_object_T_only_subcluster_1_2_3_0_sub@meta.data$CT[
  integrated_object_T_only_subcluster_1_2_3_0_sub@meta.data$sub.cluster %in% c("1_0", "1_1", "2_0", "3_0", "3_1", "4", "5", "6", "7")] <- "CD8+ T"


table(integrated_object_T_only_subcluster_1_2_3_0_sub@meta.data$CT)


# Pull out CT from all metadata objects ensuring rownames are maintained 
# from all except full object
meta_NK_integrated_obj <- NK_integrated_obj@meta.data
colnames(meta_NK_integrated_obj)

meta_mono_extra_integrated_obj <- mono_extra_integrated_obj@meta.data
colnames(meta_mono_extra_integrated_obj)

meta_proliferating_extra_integrated_obj <- proliferating_extra_integrated_obj@meta.data
colnames(meta_proliferating_extra_integrated_obj)

meta_plasma_extra_integrated_obj <- plasma_extra_integrated_obj@meta.data
colnames(meta_plasma_extra_integrated_obj)

# T cell object
meta_integrated_object_T_only_subcluster_1_2_3_0_sub <- integrated_object_T_only_subcluster_1_2_3_0_sub@meta.data
colnames(meta_integrated_object_T_only_subcluster_1_2_3_0_sub)


# Merge these subsetted object metadata with meta from the full object
# and add in the rest of the cluster labels
cols_to_keep_tmp <- intersect(colnames(meta_NK_integrated_obj), colnames(meta_mono_extra_integrated_obj))
cols_to_keep_tmp2 <- intersect(cols_to_keep_tmp, colnames(meta_proliferating_extra_integrated_obj))
cols_to_keep_tmp3 <- intersect(cols_to_keep_tmp2, colnames(meta_plasma_extra_integrated_obj))
cols_to_keep_tmp4 <- intersect(cols_to_keep_tmp3, colnames(meta_integrated_object_T_only_subcluster_1_2_3_0_sub))

meta_NK_integrated_obj <- meta_NK_integrated_obj %>% dplyr::select(CT)
meta_mono_extra_integrated_obj <- meta_mono_extra_integrated_obj %>% dplyr::select(CT)
meta_proliferating_extra_integrated_obj <- meta_proliferating_extra_integrated_obj %>% dplyr::select(CT)
meta_plasma_extra_integrated_obj <- meta_plasma_extra_integrated_obj %>% dplyr::select(CT)
meta_integrated_object_T_only_subcluster_1_2_3_0_sub <- meta_integrated_object_T_only_subcluster_1_2_3_0_sub %>% dplyr::select(CT)

meta_all <- rbind(meta_NK_integrated_obj, meta_mono_extra_integrated_obj, meta_proliferating_extra_integrated_obj,
                  meta_plasma_extra_integrated_obj, meta_integrated_object_T_only_subcluster_1_2_3_0_sub)

# Merge this metadata with pbmc_integrated_sbclstr_4_12 metadata by row name which is by = 0
meta_integrated_obj <- pbmc_integrated_sbclstr_4_12@meta.data
nrow(meta_integrated_obj)
nrow(meta_all)

meta_integrated_obj_merge <- merge(meta_integrated_obj, meta_all, by = 0, all = T)
nrow(meta_integrated_obj)

table(meta_integrated_obj_merge$CT)


####  Add in cell types ---
# Add in the other cell types from full object
# Full object CTs:
# CD14+ Monocytes: 0, 9, 4_0
# FCGR3A+ Monocytes: 10
# B: 8, 11
# Plasma: 12_0
# DC: 4_1
# pDC: 12_1
#pbmc_integrated_sbclstr_4_12@meta.data$CT <- pbmc_integrated_sbclstr_4_12@meta.data$sub.cluster

meta_integrated_obj_merge$CT[
  meta_integrated_obj_merge$sub.cluster %in% c("0", "9", "4_0")] <- "CD14+ Monocytes"
meta_integrated_obj_merge$CT[
  meta_integrated_obj_merge$sub.cluster %in% c("10")] <- "FCGR3A+ Monocytes"
meta_integrated_obj_merge$CT[
  meta_integrated_obj_merge$sub.cluster %in% c("8", "11")] <- "B cells"
meta_integrated_obj_merge$CT[
  meta_integrated_obj_merge$sub.cluster %in% c("12_0")] <- "Plasma"
meta_integrated_obj_merge$CT[
  meta_integrated_obj_merge$sub.cluster %in% c("4_1")] <- "cDC"
meta_integrated_obj_merge$CT[
  meta_integrated_obj_merge$sub.cluster %in% c("12_1")] <- "pDC"

table(meta_integrated_obj_merge$CT)

# Merge this metadata back into full object
tmp <- pbmc_integrated_sbclstr_4_12@meta.data
tmp_join <- full_join(tmp, meta_integrated_obj_merge)
pbmc_integrated_sbclstr_4_12@meta.data$CT <- tmp_join$CT

table(pbmc_integrated_sbclstr_4_12@meta.data$CT)

DimPlot(pbmc_integrated_sbclstr_4_12, group.by = "sub.cluster", 
        reduction = "integrated_sct_umap",
        label = T) 
DimPlot(pbmc_integrated_sbclstr_4_12, group.by = "CT", label = T
)


#------------------------------------------------------------------------------#
### Add CAR-T info to full object ----
#------------------------------------------------------------------------------#
DefaultAssay(pbmc_integrated_sbclstr_4_12) <- "RNA"

pbmc_integrated_sbclstr_4_12@meta.data$IL13OPCounts <- pbmc_integrated_sbclstr_4_12@assays$RNA@counts["IL13OP",]
pbmc_integrated_sbclstr_4_12@meta.data$CART <- ifelse(
  (pbmc_integrated_sbclstr_4_12@meta.data$IL13OPCounts >= 3 & (pbmc_integrated_sbclstr_4_12@meta.data$CT == "CD8+ T" | pbmc_integrated_sbclstr_4_12@meta.data$CT == "CD4+ T" | pbmc_integrated_sbclstr_4_12@meta.data$CT == "Treg" )), "Positive", "Negative")

table(pbmc_integrated_sbclstr_4_12@meta.data$CART)


#------------------------------------------------------------------------------#
### Save object ----
#------------------------------------------------------------------------------#
saveRDS(pbmc_integrated_sbclstr_4_12, 
        "/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/PBMC_all_cells_obj_2023.rds")


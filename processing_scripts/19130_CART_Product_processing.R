#==============================================================================#
# Author : Angela M. Oill, aoill@tgen.org
# Date: 2023/06/09
# Project: Pediatric CAR-T 
# Description: Preprocessing and integration of product samples
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

# Keep only product samples
product_seurat_list_demultiplexed_singlets <- list()
for (i in 1:length(seurat_list_demultiplexed_singlets)) {
  print(i)
  sub_i <- subset(seurat_list_demultiplexed_singlets[[i]], subset = Sample_Type == "Product")
  product_seurat_list_demultiplexed_singlets <- c(product_seurat_list_demultiplexed_singlets, sub_i)
}
names(product_seurat_list_demultiplexed_singlets) <- c("Batch37_filtered", "Batch38_filtered", "Batch39_filtered")


#==============================================================================#
# Read in product sample that wasn't multiplexed (626) ----
#==============================================================================#
sample_metadata_not_multiplexed_626 <- sample_metadata_not_multiplexed %>% filter(FID_GEXFB == "F05565")
seurat_list_626 <- prep_seurat_list(sample_metadata_not_multiplexed_626, batch_ID = "FID_GEXFB", cellRanger_path = "CellRanger_path", cell_ID_prefix = "Cell_Prefix", run_soupX = FALSE)

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
  DropletUtils::write10xCounts(path = paste0("/scratch/aoill/projects/CAR-T/00_final/soupX/demultiplexed_",names(seurat_list_separated[i])), 
                               x = obj.sub[["RNA"]]@data, 
                               barcodes = colnames(obj.sub[["RNA"]]@data), # cell names
                               gene.id = rownames(obj.sub[["RNA"]]@data), # Gene identifiers, one per row of X
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
    sample_seurat_list[[i]] <- SCTransform(sample_seurat_list[[i]], 
                                           vst.flavor = "v2", verbose = F)
    DefaultAssay(sample_seurat_list[[i]]) <- "SCT"
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
# Remove RB and MT genes ----
#==============================================================================#
# Remove ribosomal and mt genes
product_seurat_list_separated_SoupX_noRBSMT <- list()
message("Removing ribosomal and mitochondrial genes from Seurat object")
for (i in 1:length(product_seurat_list_separated_SoupX)) {
  RBMTgenes <- grep(pattern = "^RP[SL]|^MRP[SL]|^MT-", x = rownames(product_seurat_list_separated_SoupX[[i]]@assays$RNA@data), value = TRUE, invert = TRUE)
  product_seurat_list_separated_SoupX_noRBSMT[[i]] = subset(product_seurat_list_separated_SoupX[[i]], features = RBMTgenes)
  
}
names(product_seurat_list_separated_SoupX_noRBSMT) <- names(product_seurat_list_separated_SoupX)


#==============================================================================#
# Integration ----
#==============================================================================#

#------------------------------------------------------------------------------#
## Add some more metadata ----
#------------------------------------------------------------------------------#
for (i in 1:length(product_seurat_list_separated_SoupX_noRBSMT)) {
  product_seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$UPN_Cycle_Day_Sample_Type_Batch <- paste0(product_seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$UPN, "_",
                                                                                                       product_seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$Cycle,  "_", 
                                                                                                       product_seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$Day, "_",
                                                                                                       product_seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$Sample_Type, "_",
                                                                                                       product_seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$Batch)
  
  product_seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$Cycle_Day <- paste0("Cycle", product_seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$Cycle,  "_Day", 
                                                                                 product_seurat_list_separated_SoupX_noRBSMT[[i]]@meta.data$Day)
}



#------------------------------------------------------------------------------#
## Integrate Product ----
#------------------------------------------------------------------------------#
product_seurat_list <- product_seurat_list_separated_SoupX_noRBSMT

product_seurat_list <- lapply(product_seurat_list, function(xx){
  # Rerunning SCTransform
  DefaultAssay(xx) <- "SoupX_RNA"
  xx <- SCTransform(xx, 
                    method = "glmGamPoi",
                    vars.to.regress = c("S.Score", "G2M.Score"), 
                    vst.flavor = "v2",
                    verbose = T)
  
})

product_integrated_NEW <- basic_sct_rpca_integration(product_seurat_list, npcs=50, k_weight=100, numfeatures = 2500)


# Remove cluster 8 and re-cluster
# This is an LYZ cluster
Idents(product_integrated_NEW) <- "integrated_sct_snn_res.0.5"
Idents(product_integrated_NEW)
product_integrated_NEW_subset <- subset(product_integrated_NEW, 
                                        idents = c("0", "1", "2", "3", "4", "5",
                                                   "6", "7", "9"))
table(product_integrated_NEW_subset@meta.data$integrated_sct_snn_res.0.5)

# Re-cluster
npcs <- min(get_pcs(product_integrated_NEW_subset, reduction_name = "integrated_sct_pca"))
npcs
integrated_object <- RunUMAP(product_integrated_NEW_subset,
                             reduction = "integrated_sct_pca",
                             reduction.name = "integrated_sct_umap",
                             dims = 1:npcs,
                             return.model = TRUE)
integrated_object <- FindNeighbors(integrated_object,
                                   reduction = "integrated_sct_pca",
                                   dims = 1:npcs,
                                   graph.name = c("integrated_sct_nn",
                                                  "integrated_sct_snn"))
integrated_object <- FindClusters(integrated_object,
                                  resolution = c(0.1,0.2,0.3,0.5,0.8,1),
                                  graph.name = "integrated_sct_snn")


#==============================================================================#
# Add metadata (IN PROGRESS) ----
#==============================================================================#

#------------------------------------------------------------------------------#
### Add T cell state annotations to T cell object ----
#------------------------------------------------------------------------------#
#### Sub-cluster ----
# sub-cluster 1
Idents(integrated_object) <- "integrated_sct_snn_res.1"
integrated_object_sbclst_1 <- FindSubCluster(integrated_object, "1", resolution = 0.2, graph.name = "integrated_sct_snn")

# sub-cluster 6
Idents(integrated_object_sbclst_1) <- "sub.cluster"
integrated_object_sbclst_1_6 <- FindSubCluster(integrated_object_sbclst_1, "6", resolution = 0.25, graph.name = "integrated_sct_snn")

# sub-cluster 7
Idents(integrated_object_sbclst_1_6) <- "sub.cluster"
integrated_object_sbclst_1_6_7 <- FindSubCluster(integrated_object_sbclst_1_6, "7", resolution = 0.3, graph.name = "integrated_sct_snn")

# sub-cluster 11
Idents(integrated_object_sbclst_1_6_7) <- "sub.cluster"
integrated_object_sbclst_1_6_7_11 <- FindSubCluster(integrated_object_sbclst_1_6_7, "11", resolution = 0.3, graph.name = "integrated_sct_snn")

# sub-cluster 12
Idents(integrated_object_sbclst_1_6_7_11) <- "sub.cluster"
integrated_object_sbclst_1_6_7_11_12 <- FindSubCluster(integrated_object_sbclst_1_6_7_11, "12", resolution = 0.2, graph.name = "integrated_sct_snn")

#### Add cell state annotations ----
Idents(integrated_object_sbclst_1_6_7_11_12) <- "sub.cluster"

integrated_object_sbclst_1_6_7_11_12@meta.data$ct_T_states <- integrated_object_sbclst_1_6_7_11_12@meta.data$sub.cluster

integrated_object_sbclst_1_6_7_11_12@meta.data$ct_T_states[
  integrated_object_sbclst_1_6_7_11_12@meta.data$sub.cluster %in% c("1_0", "1_1", "11_0", "11_2", "14", "15", 
                                                                    "16", "5", "6_1", "6_2", "7_0")] <- "CD8 Proliferating"
integrated_object_sbclst_1_6_7_11_12@meta.data$ct_T_states[
  integrated_object_sbclst_1_6_7_11_12@meta.data$sub.cluster %in% c("1_2", "11_1", "7_1", "7_2")] <- "CD4 Proliferating"
integrated_object_sbclst_1_6_7_11_12@meta.data$ct_T_states[
  integrated_object_sbclst_1_6_7_11_12@meta.data$sub.cluster %in% c("3", "10")] <- "CD8 Effector"
integrated_object_sbclst_1_6_7_11_12@meta.data$ct_T_states[
  integrated_object_sbclst_1_6_7_11_12@meta.data$sub.cluster %in% c("8", "12_1")] <- "CD8 Exhausted"
integrated_object_sbclst_1_6_7_11_12@meta.data$ct_T_states[
  integrated_object_sbclst_1_6_7_11_12@meta.data$ct_T_states %in% c("0")] <- "CD8 Central Memory (Activated)"
integrated_object_sbclst_1_6_7_11_12@meta.data$ct_T_states[
  integrated_object_sbclst_1_6_7_11_12@meta.data$ct_T_states %in% c("4", "12_0")] <- "CD4 Memory (LTB+)"
integrated_object_sbclst_1_6_7_11_12@meta.data$ct_T_states[
  integrated_object_sbclst_1_6_7_11_12@meta.data$sub.cluster %in% c("2", "6_0")] <- "CD4 Memory"
integrated_object_sbclst_1_6_7_11_12@meta.data$ct_T_states[
  integrated_object_sbclst_1_6_7_11_12@meta.data$sub.cluster %in% c("9")] <- "CD8 T (OX40+)" #TNFRSF4
integrated_object_sbclst_1_6_7_11_12@meta.data$ct_T_states[
  integrated_object_sbclst_1_6_7_11_12@meta.data$sub.cluster %in% c("13")] <- "NKT" 


table(integrated_object_sbclst_1_6_7_11_12@meta.data$ct_T_states)


#------------------------------------------------------------------------------#
### Add CART information to T cell object ----
#------------------------------------------------------------------------------#
# Add IL13OP info to metadata
DefaultAssay(integrated_object_sbclst_1_6_7_11_12) <- "RNA"

integrated_object_sbclst_1_6_7_11_12@meta.data$IL13OPCounts <- integrated_object_sbclst_1_6_7_11_12@assays$RNA@counts["IL13OP",]
integrated_object_sbclst_1_6_7_11_12@meta.data$CART <- ifelse(integrated_object_sbclst_1_6_7_11_12@meta.data$IL13OPCounts >= 3, "Positive", "Negative")
table(integrated_object_sbclst_1_6_7_11_12@meta.data$CART)

#------------------------------------------------------------------------------#
### Add TCR data to T cell object ----
#------------------------------------------------------------------------------#
#### Add TCR data determined from scRepertoire ----

##### Load the contig data ----
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


##### Process the contig data ----

product_obj <- integrated_object

product_list <- SplitObject(product_obj, split.by = "Batch")

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

##### Combine the contigs ----

# Will want to use the parameter removeNA set to true to remove any cell barcode 
# with an NA value in at least one of the chains
combined <- combineTCR(contig_list_product, 
                       samples =  c("F05037_515_Product", "F05041_574_Product", 
                                    "F05041_514_Product", "F05625_625_Product", 
                                    "F05629_626_Product"),
                       cells ="T-AB", filterMulti = FALSE, removeNA = TRUE
)


##### Add TCR info to the Seurat objects ----

### Add a new barcod column in the TCR list to match the Seurat objects  ----
# I need to replace the barcode column in combined which matching barcode to 
# seurat object 
# Path to google sheet with metadata for your samples
metapath <- "https://docs.google.com/spreadsheets/d/1fmEuOFXm893mfS2T1zpEsMrN7MpFc_Aym24ZODbw2bA/edit#gid=0"
sheet_name <- "Sheet1"
sample_metadata <- load_metadata(metadata_path = metapath, data_type = "googlesheet", sheet_name = sheet_name)

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

obj_T_fltr_rclstr_sbclstr_fnl_TCR <- combineExpression(combined2, integrated_object_sbclst_1_6_7_11_12,
                                                       #group.by = "UPN_Cycle_Day_Sample_Type",
                                                       group.by = "sample",
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
saveRDS(combined2, "/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/product_TCR_clonotype_scRepertoire.rds")



#### Get normalized TCR frequencies and categories ----
# TCRs with a frequency of 1 will be categorized as non-expanded 
# Distribution of TCR frequencies normalized will only include TCRs with a 
# frequency > 1

# calculate Frequency_norm
tmp_meta <- obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data

# Add a T cell number count for each UPN patient and cycle
tmp_meta <- obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data %>%
  add_count(UPN, name = "UPN_n")

# Get normalized TCR frequency
tmp_meta$Frequency_norm <- tmp_meta$Frequency/tmp_meta$UPN_n

obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$UPN_n <- tmp_meta$UPN_n
obj_T_fltr_rclstr_sbclstr_fnl_TCR@meta.data$Frequency_norm <- tmp_meta$Frequency_norm


tmp_meta_sub <- tmp_meta %>%
  dplyr::select(UPN, Sample_Type, CART, CTstrict, Frequency,
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
        "/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/product_T_cell_obj_2023.rds")

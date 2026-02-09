#==============================================================================#
# Author(s) : Angela M. Oill (aoill@tgen.org) 
#             Annika Vannan (avannan@tgen.org)
#             Heini Natri (hnatri@tgen.org)
#             Linh Bui (lbui@tgen.org)
#
# Date: 2023/02/13
#
# Project: scRNA processing pipeline
#
# Description: Module 1 - Pre-processing and QC
#==============================================================================#



#==============================================================================#
# Load Libraries ----
#==============================================================================#
suppressMessages(library(Seurat))
suppressMessages(library(googlesheets4))
suppressMessages(library(SoupX))
suppressMessages(library(scater))
suppressMessages(library(DoubletFinder))
#suppressMessages(library(scDblFinder))
suppressMessages(library(plyr))



#==============================================================================#
# Functions ----
#==============================================================================#


#---------------------------#
## Load metadata into R ----
#---------------------------#
# Description: The purpose of this function is to read metadata associated with
#              single cell data into R as a data frame.
# Input(s):
#         metadata_path: This is the path to the users metadata. Can be a link 
#                        if it is a google sheet or a path if it is a tsv/csv.
#                        This path should be called in the users master script.
#         data_type: This is the format of the users metadata. Accepted formats:
#                    "tsv", "csv", "googlesheet".
#         sheet_name: This is the specific sheet from the google sheet where the
#                     metadata is stored. 
# Output(s): 
#         metadata: Data frame containing metadata associated with the users
#                   single cell data.
load_metadata <- function(metadata_path, data_type = "googlesheet", sheet_name = "Sheet1"){
  if (data_type == "tsv"){
    metadata <- read.table(metadata_path, header = T, sep = "\t")
    return(metadata)
  } else if (data_type == "csv"){
    metadata <- read.table(metadata_path, header = T, sep = ",")
    return(metadata)
  } else if (data_type == "googlesheet"){
    metadata <- gs4_get(metadata_path)
    sheet_names(metadata)
    metadata <- read_sheet(metadata, sheet = sheet_name)
    return(metadata)
  }
}


#------------------------------------#
## Prepare list of Seurat objects ----
#------------------------------------#
# Description: This function will load 10x data by user specified batch, 
#              calculate percent mt and ribosomal genes, and add a prefix to
#              cell IDs. Optionally, this function can also run SoupX and add 
#              the corrected counts to each object in the list.
# Input(s):
#          metadata: Data frame containing metadata associated with the 
#                    users single cell data.
#          batch_ID: Column name in "metadata" that contains batch ID. This 
#                    corresponds to each sequence sample to make a Seurat object 
#                    for.
#          cellRanger_path: Column name in "metadata" that contains the path to 
#                           the Cellranger data. Do not include 
#                           "/outs/filtered_feature_bc_matrix/" or 
#                           "/outs/raw_feature_bc_matrix/" in the path name.
#          cell_ID_prefix: Column name in "metadata" with a prefix to add to the
#                          cell IDs. Make sure these are unique to each 10x data
#                          that will be in the Seurat list.
#          run_soupX: Optional. If specify TRUE, SoupX will be run. Default is 
#                     FALSE.
# Output(s):
#          sample_seurat_list: This is the pre-processed list of Seurat objects.
prep_seurat_list <- function(metadata, batch_ID = "Batch_ID", cellRanger_path = "CellRanger_path", cell_ID_prefix = "Batch_ID", run_soupX = FALSE){
  sample_list <- unique(metadata[[batch_ID]])
  
  prefix_list <- metadata[[cell_ID_prefix]]
  if (length(prefix_list) != length(unique(prefix_list))) {
    message("WARNING: Cell ID prefixes are not unique across samples!!")
  }
  
  if (run_soupX == FALSE){
    sample_seurat_list <- lapply(sample_list, function(i){
      message(i)
      sample_10x_data <- Read10X(paste0(metadata[which(metadata[[batch_ID]]==i),][[cellRanger_path]], "/outs/filtered_feature_bc_matrix/"))
      seurat_object <- CreateSeuratObject(counts = sample_10x_data)
      seurat_object <- PercentageFeatureSet(seurat_object, pattern = "^MT-", col.name = "percent.mt_RNA", assay = "RNA")
      seurat_object <- PercentageFeatureSet(seurat_object, pattern = "^RP[SL]|^MRP[SL]", col.name = "percent.ribo_RNA", assay = "RNA")
      # Rename cell IDs, adding a prefix specified by the user
      seurat_object <- RenameCells(seurat_object,
                                   new.names = paste0(metadata[which(metadata[[batch_ID]]==i),][[cell_ID_prefix]], 
                                                      "_", colnames(seurat_object)))
      
    })
    names(sample_seurat_list) <- sample_list
    
    # add metadata
    for (i in names(sample_seurat_list)) {
      meta.data <- metadata[which(metadata[[batch_ID]]==i),]
      # Add sample metadata
      for (m in colnames(metadata)){
        if (m == cellRanger_path) {
          next
        } else {
          sample_seurat_list[[i]]@meta.data[[m]] = meta.data[[m]]
          }
      }
    }
    
    return(sample_seurat_list) 
    
  } else if (run_soupX == TRUE){
    sample_seurat_list <- sapply(sample_list,  function(i){
      # Read in count and droplet data
      d10x_toc <- Read10X(paste0(metadata[which(metadata[[batch_ID]]==i),][[cellRanger_path]], "/outs/filtered_feature_bc_matrix/"))
      d10x_tod <- Read10X(paste0(metadata[which(metadata[[batch_ID]]==i),][[cellRanger_path]], "/outs/raw_feature_bc_matrix/"))
      
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
      sc <- autoEstCont(sc)
      out <- adjustCounts(sc)
      
      # Create Seurat object using corrected data
      d10x_seu <- CreateSeuratObject(out, assay = "SoupX_RNA")
      d10x_seu[["RNA"]] <- toc_seu@assays[["RNA"]]
      d10x_seu <- PercentageFeatureSet(d10x_seu, pattern = "^MT-", col.name = "percent.mt_RNA", assay = "RNA")
      d10x_seu <- PercentageFeatureSet(d10x_seu, pattern = "^RP[SL]|^MRP[SL]", col.name = "percent.ribo_RNA", assay = "RNA")
      d10x_seu <- PercentageFeatureSet(d10x_seu, pattern = "^MT-", col.name = "percent.mt_SoupX_RNA", assay = "SoupX_RNA")
      d10x_seu <- PercentageFeatureSet(d10x_seu, pattern = "^RP[SL]|^MRP[SL]", col.name = "percent.ribo_SoupX_RNA", assay = "SoupX_RNA")
      
      # Rename cell IDs, adding a prefix specified by the user
      d10x_seu <- RenameCells(d10x_seu,
                                   new.names = paste0(metadata[which(metadata[[batch_ID]]==i),][[cell_ID_prefix]], 
                                                      "_", colnames(d10x_seu)))
    })
    names(sample_seurat_list) <- sample_list
    
    # add metadata
    for (i in names(sample_seurat_list)) {
      meta.data <- metadata[which(metadata[[batch_ID]]==i),]
      # Add sample metadata
      for (m in colnames(metadata)){
        if (m == cellRanger_path) {
          next
        } else {
          sample_seurat_list[[i]]@meta.data[[m]] = meta.data[[m]]
        }
      }
    }
    
    return(sample_seurat_list)
  }
}


#------------------------------------------------------------#
## Prepare list of Seurat objects for multiplexed samples ----
#------------------------------------------------------------#
# Description: This function will load 10x multiplexed data by user specified 
#              batch, calculate percent mt and ribosomal genes, add a prefix to
#              cell IDs, run cell hashing, and add metadata to the list of 
#              Seurat objects. 
# Input(s):
#          metadata: Data frame containing metadata associated with the 
#                    users single cell data.
#          batch_ID: Column name in "metadata" that contains batch ID. This 
#                    corresponds to each sequence sample to make a Seurat object 
#                    for.
#          cellRanger_path: Column name in "metadata" that contains the path to 
#                           the Cellranger data. Do not include 
#                           "/outs/filtered_feature_bc_matrix/" or 
#                           "/outs/raw_feature_bc_matrix/" in the path name.
#          cell_ID_prefix: Column name in "metadata" with a prefix to add to the
#                          cell IDs. Make sure these are unique to each 10x data
#                          that will be in the Seurat list.
#          CellHashing_Ab: Column name in "metadata" that has the cell hashing 
#                          antibody names.
# Output(s):
#          sample_seurat_list: This is the pre-processed list of Seurat objects.
prep_seurat_list_multiplexed <- function(metadata, batch_ID = "Batch", 
                                         cellRanger_path = "CellRanger_path", 
                                         cell_ID_prefix = "Batch_ID", 
                                         CellHashing_Ab = "CellHashing_Ab", 
                                         p_quantile = 0.99){
  sample_list <- unique(metadata[[batch_ID]])
  
  message("Reading in 10X files")
  sample_seurat <- lapply(sample_list, function(i){
    #message(i)
    sample_10x_data <- Read10X(unique(paste0(metadata[which(metadata[[batch_ID]]==i),][[cellRanger_path]], "/outs/filtered_feature_bc_matrix/")))
    
    # Fix antibody names
    rownames(sample_10x_data$`Antibody Capture`) = strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\2', rownames(sample_10x_data$`Antibody Capture`)), ' ')
    rownames(sample_10x_data$`Antibody Capture`) = str_remove(rownames(sample_10x_data$`Antibody Capture`), "mouse_")
    rownames(sample_10x_data$`Antibody Capture`) = str_remove(rownames(sample_10x_data$`Antibody Capture`), "rat_")
    rownames(sample_10x_data$`Antibody Capture`) = str_remove(rownames(sample_10x_data$`Antibody Capture`), "human_")
    
    sample_10x_data
  })
  names(sample_seurat) <- sample_list
  
  hash_antibodies = c("TotalSeqC0251_Hashtag1", 
                      "TotalSeqC0252_Hashtag2",  
                      "TotalSeqC0253_Hashtag3", 
                      "TotalSeqC0254_Hashtag4",
                      "TotalSeqC0256_Hashtag6",
                      "TotalSeqC0257_Hashtag7", 
                      "TotalSeqC0258_Hashtag8",
                      "TotalSeqC0259_Hashtag9")
  
  # Create a batch list 
  # NOTE: no CITEseq, but the protein assay was created
  message("Creating a list of Seurat objects, pre cell-hashing")
  sample_seurat_list = list()
  for (i in 1:length(sample_seurat)) {
    #message(names(sample_seurat)[i])
    # Set up the Seurat object
    # Splitting out the gene expression and the antibody matrices
    GEX = sample_seurat[[i]][[1]]
    AB = sample_seurat[[i]][[2]]
    # Creating the object with the expression matrix
    sample_seurat_list[[i]] = CreateSeuratObject(counts = GEX)
    sample_seurat_list[[i]] = PercentageFeatureSet(sample_seurat_list[[i]], pattern = "^MT-", col.name = "percent.mt_RNA", assay = "RNA")
    sample_seurat_list[[i]] = PercentageFeatureSet(sample_seurat_list[[i]], pattern = "^RP[SL]|^MRP[SL]", col.name = "percent.ribo_RNA", assay = "RNA")
    
    hash = AB[rownames(AB) %in% hash_antibodies,]
    citeseq = AB[!(rownames(AB) %in% hash_antibodies),]
    
    # Creating the protein and hash assays
    sample_seurat_list[[i]][["Protein"]] = CreateAssayObject(counts = citeseq)
    sample_seurat_list[[i]][["Hash"]] = CreateAssayObject(counts = hash)
    
    # Rename cell IDs, adding a prefix specified by the user
    sample_seurat_list[[i]] <- RenameCells(sample_seurat_list[[i]],
                                           new.names = paste0(unique(metadata[which(metadata[[batch_ID]]==names(sample_seurat)[i]),][[cell_ID_prefix]]), 
                                                              "_", colnames(sample_seurat_list[[i]])))
  }
  names(sample_seurat_list) <- sample_list
  
  # Normalizing and run cell hashing
  # Metadata will also be added in this loop
  message("Starting cell hashing")
  for (i in 1:length(sample_seurat_list)) {
    # Keeping antibodies that were used in the panel
    batch_num <- names(sample_seurat)[[i]]
    #message(batch_num)
    
    meta.data <- metadata[which(metadata[[batch_ID]]==names(sample_seurat)[[i]]),]
    meta.data[[CellHashing_Ab]] <- gsub("_", "-", meta.data[[CellHashing_Ab]])
    #print(meta.data)
    
    #sample_seurat_list[[i]][["Hash"]] = CreateAssayObject(counts = sample_seurat_list[[i]]@assays$Hash[which(rownames(sample_seurat_list[[i]]@assays$Hash) %in% meta.data[[CellHashing_Ab]]),])
    sample_seurat_list[[i]][["Hash"]] = CreateAssayObject(counts = sample_seurat_list[[i]]@assays$Hash$counts[which(rownames(sample_seurat_list[[i]]@assays$Hash) %in% meta.data[[CellHashing_Ab]]),])
    
    # Normalizing and scaling
    sample_seurat_list[[i]] = NormalizeData(sample_seurat_list[[i]], assay = "Protein", normalization.method = "CLR")
    sample_seurat_list[[i]] = ScaleData(sample_seurat_list[[i]], assay = "Protein")
    sample_seurat_list[[i]] = NormalizeData(sample_seurat_list[[i]], assay = "Hash", normalization.method = "CLR")
    sample_seurat_list[[i]] = ScaleData(sample_seurat_list[[i]], assay = "Hash")
    
    # Assign single cells back to their sample origins
    # Demultiplexing based on the hashing antibodies
    # Singlets will be kept based on "hash classification global"
    sample_seurat_list[[i]] = HTODemux(sample_seurat_list[[i]], assay = "Hash", positive.quantile = p_quantile, verbose = F)
    
    message("Cell hashing complete")
    
    # Add sample metadata
    message("Adding metadata")
    for (m in colnames(metadata)){
      if (m == cellRanger_path) {
        next
      } else {
        sample_seurat_list[[i]]@meta.data[[m]] = plyr::mapvalues(x = sample_seurat_list[[i]]@meta.data$hash.ID,
                                                                 from = meta.data$CellHashing_Ab,
                                                                 to = as.character(meta.data[[m]]))
      }
    }
  }
  
  
  return(sample_seurat_list) 
}


#--------------------------------------------------------------#
## Filter out doublets and negative cells post CellHashing ----
#--------------------------------------------------------------#
# Description: This function will filter out doublets and negative cells from
#              each Seurat object in a list where de-multiplexing using 
#              CellHashing was performed. Since cells were filtered certain
#              metadata columns might loose factors, resulting in empty factors
#              in the metadata, thus unused levels in metadata are also dropped 
#              in this function.
# Input(s):
#          sample_seurat_list: list of Seurat objects where de-multiplexing using
#                              CellHashing was performed.
# Output(s):
#          sample_seurat_list_filtered: List of Seurat objects with only 
#                                       singlets.
filter_doublets_cellhashing <- function(sample_seurat_list){
  sample_seurat_list_filtered = list()
  for (i in 1:length(sample_seurat_list)) { 
    # Keeping singlets only
    message("Keeping singlets only")
    Idents(sample_seurat_list[[i]]) = sample_seurat_list[[i]]$Hash_classification.global
    sample_seurat_list_filtered[[i]] = subset(sample_seurat_list[[i]], idents = "Singlet")
    
    # Drop the unused levels in metadata
    message("Drop unused levels in metadata")
    for (m in base::colnames(sample_seurat_list_filtered[[i]]@meta.data)){
      if (base::is.factor(sample_seurat_list_filtered[[i]]@meta.data[[m]])){
        sample_seurat_list_filtered[[i]]@meta.data[[m]] <- base::droplevels(sample_seurat_list_filtered[[i]]@meta.data[[m]])  # need to drop levels of the removed values
      }
    }
  }
  
  batch_filt_names = paste0(names(sample_seurat_list), "_filtered")
  names(sample_seurat_list_filtered) = batch_filt_names
  
  
  return(sample_seurat_list_filtered)
}


#---------------------#
## Plot QC metrics ----
#---------------------#
# Description:
# 
# Input(s):
# 
# Output(s):
#
plot_QC <- function(sample_seurat_list){
  # First merge the list of Seurat objects
  seurat_merge <- merge(x = sample_seurat_list[[1]], y = sample_seurat_list[2:length(sample_seurat_list)])
  
  # Make a violin plots
  p1 <- VlnPlot(seurat_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt_RNA", "percent.ribo_RNA"), ncol = 2, pt.size = 0)
  print(p1)
  
  # Make density plots
  smoothScatter(seurat_merge@meta.data$percent.mt, seurat_merge@meta.data$nFeature_RNA, ylab = "nFeature_RNA", xlab = "% MT")

  smoothScatter(seurat_merge@meta.data$percent.mt, seurat_merge@meta.data$nCount_RNA, ylab = "nCount_RNA", xlab = "% MT")
}


#-----------------------------------#
## Annotated smooth scatter plot ----
#-----------------------------------#
plot_smoothscatter_annot <- function(sample_seurat_list){
  # First merge the list of Seurat objects
  seurat_merge <- merge(x = sample_seurat_list[[1]], y = sample_seurat_list[2:length(sample_seurat_list)])
  
  # Calculate medians and mads for log2-transformed values
  median_log_nCount <- median(log2(seurat_merge$nCount_RNA))
  median_log_nFeat <- median(log2(seurat_merge$nFeature_RNA))
  median_log_perc_mt <- median(log2(seurat_merge$percent.mt_RNA))
  mad_log_nCount <- mad(log2(seurat_merge$nCount_RNA))
  mad_log_nFeat <- mad(log2(seurat_merge$nFeature_RNA))
  #mad_log_perc_mt <- mad(log2(seurat_merge$percent.mt_RNA))
  mad_log_perc_mt <- log2(10)
  
  # Mads line markings
  lines_log_nCount <- c(median_log_nCount + (mad_log_nCount)*c(seq(2, 5, 1), seq(-2, -5, -1)))
  lines_unlog_nCount <- 2^(lines_log_nCount)
  lines_log_nFeat <- c(median_log_nFeat + (mad_log_nFeat)*c(seq(2, 5, 1), seq(-2, -5, -1)))
  lines_unlog_nFeat <- 2^(lines_log_nFeat)
  #lines_unlog_perc_mt <- 2^(mad_log_perc_mt + (mad_log_perc_mt)*c(seq(2, 5, 1), seq(-2, -5, -1)))
  lines_unlog_perc_mt <- 10
  line_colors = c(rep(c("red", "blue", "green", "orange"), 2),
                  rep(c("red", "blue", "green", "orange"), 2))
  
  # nCount vs. nFeature
  smoothScatter(log2(seurat_merge$nCount_RNA), log2(seurat_merge$nCount_RNA),
                xlab = "log2(nCount_RNA)", ylab = "log2(nFeature_RNA)")
  abline(v = lines_log_nCount, h = lines_log_nFeat,       
         lty = "dashed", lwd = 1.25, col = line_colors)
  smoothScatter(seurat_merge$nCount_RNA, seurat_merge$nCount_RNA,
                xlab = "nCount_RNA", ylab = "nFeature_RNA")
  abline(v = lines_unlog_nCount, h = lines_unlog_nFeat,     
         lty = "dashed", lwd = 1.25, col = line_colors)
  
  # nCount vs. percent.mt_RNA
  smoothScatter(seurat_merge$percent.mt_RNA, log2(seurat_merge$nCount_RNA),
                xlab = "% MT", ylab = "log2(nCount_RNA)")
  abline(v = lines_unlog_perc_mt, h = lines_log_nCount, 
         lty = "dashed", lwd = 1.25, col = line_colors)
  smoothScatter(seurat_merge$percent.mt_RNA, seurat_merge$nCount_RNA,
                xlab = "% MT", ylab = "nCount_RNA")
  abline(v = lines_unlog_perc_mt, h = lines_unlog_nCount,
         lty = "dashed", lwd = 1.25, col = line_colors)
  
  # nFeature vs. percent.mt_RNA
  smoothScatter(seurat_merge$percent.mt_RNA, log2(seurat_merge$nFeature_RNA),
                xlab = "% MT", ylab = "log2(nFeature_RNA)")
  abline(v = lines_unlog_perc_mt, h = lines_log_nFeat,     
         lty = "dashed", lwd = 1.25, col = line_colors)
  smoothScatter(seurat_merge$percent.mt_RNA, seurat_merge$nFeature_RNA,
                xlab = "% MT", ylab = "nFeature_RNA")
  abline(v = lines_unlog_perc_mt, h = lines_unlog_nFeat, 
         lty = "dashed", lwd = 1.25, col = line_colors)
  
  }



#-------------------------#
## Filtering with MADs ----
#-------------------------#
# Description: Filter cells based on a user specified median absolute deviations 
#              (MADs) for nFeature and nCount, and hard filters cells based on 
#              a user specified percent MT. After filtering, mitochondrial and 
#              ribosomal genes are removed.
# Input(s):
#         sample_seurat_list: List of Seurat objects to filter
#         batch_ID: Column name in "metadata" that contains batch ID. This 
#                   corresponds to each sequence sample to make a Seurat object 
#                   for.
#         pt_mt: Percent mitochondrial filter. Cells will be removed with more 
#                than the set pt_mt.
#         nFeature: MADs threshold for nFeature. If set, will set upper and lower 
#                   the same.
#         nFeature_upper: Upper MADs threshold for nFeature. 
#         nFeature_lower: Lower MADs threshold for nFeature.
#         nCount: MADs threshold for nCount If set, will set upper and lower 
#                 the same.
#         nCount_upper: Upper MADs threshold for nCount
#         nCount_lower: Lower MADs threshold for nCount
# Output(s):
#         filter_seurat_list_noRBSMT: Filtered list of Seurat objects.
filter_mads_mt_rb_genes <- function(sample_seurat_list, batch_ID = "Batch_ID", pt_mt = 20, nFeature = 5, nFeature_upper = nFeature, nFeature_lower = nFeature, nCount = 5, nCount_upper = nCount, nCount_lower = nCount){
  # merge and create a single cell object
  seurat_merge <- merge(x = sample_seurat_list[[1]], y = sample_seurat_list[2:length(sample_seurat_list)])
  pre_sce <- as.SingleCellExperiment(seurat_merge, assay = "RNA")
  
  # get list of cells to filter
  lenient_mads_outs <- colnames(pre_sce)[which(
    isOutlier(pre_sce$nCount_RNA, nmads = nCount_upper, type = "higher", log = T) |
      isOutlier(pre_sce$nCount_RNA, nmads = nCount_lower, type = "lower", log = T) |
      isOutlier(pre_sce$nFeature_RNA, nmads = nFeature_upper, type = "higher", log = T) |
      isOutlier(pre_sce$nFeature_RNA, nmads = nFeature_lower, type = "lower", log = T) |
      pre_sce$percent.mt_RNA > pt_mt)]
  
  # Remove outliers from merged Seurat object
  filter_seurat_merge <- subset(seurat_merge, cells = lenient_mads_outs, invert = TRUE)
  
  # Create split 
  filter_seurat_list <- SplitObject(filter_seurat_merge, split.by = batch_ID)
  
  # Remove ribosomal and mt genes
  filter_seurat_list_noRBSMT = list()
  
  sample_list <- names(sample_seurat_list)
  for (i in 1:length(sample_list)) {
    RBMTgenes <- grep(pattern = "^RP[SL]|^MRP[SL]|^MT-", x = rownames(filter_seurat_list[[i]]@assays$RNA@data), value = TRUE, invert = TRUE)
    filter_seurat_list_noRBSMT[[i]] = subset(filter_seurat_list[[i]], features = RBMTgenes)
    
  }
  names(filter_seurat_list_noRBSMT) <- sample_list

  return(filter_seurat_list_noRBSMT)
}


#-----------------------#
## Manual filtering ----
#-----------------------#
# NOTE: This will only do a lower filter for nFeature_RNA and nCount_RNA and an
#       upper filter for percent MT.
filter_manual_mt_rb_genes <- function(sample_seurat_list, pt_mt = 20, nFeature = 350, nCount = 450){
  filter_seurat_list <- list()
  for (i in 1:length(sample_seurat_list)) { 
    message("Filtering cells from sample ", i)
    filter_seurat_list[[i]] <- subset(sample_seurat_list[[i]], subset = nFeature_RNA > nFeature & nCount_RNA > nCount & percent.mt_RNA < pt_mt)
  }
  names(filter_seurat_list) <- names(sample_seurat_list)
  
  # Remove ribosomal and mt genes
  filter_seurat_list_noRBSMT <- list()
  message("Removing ribosomal and mitochondrial genes from Seurat object")
  for (i in 1:length(sample_seurat_list)) {
    RBMTgenes <- grep(pattern = "^RP[SL]|^MRP[SL]|^MT-", x = rownames(filter_seurat_list[[i]]@assays$RNA@data), value = TRUE, invert = TRUE)
    filter_seurat_list_noRBSMT[[i]] = subset(filter_seurat_list[[i]], features = RBMTgenes)
    
  }
  names(filter_seurat_list_noRBSMT) <- names(sample_seurat_list)
  
  return(filter_seurat_list_noRBSMT)
}

filter_manual <- function(sample_seurat_list, pt_mt = 20, nFeature = 350, nCount = 450){
  filter_seurat_list <- list()
  for (i in 1:length(sample_seurat_list)) { 
    message("Filtering cells from sample ", i)
    filter_seurat_list[[i]] <- subset(sample_seurat_list[[i]], subset = nFeature_RNA > nFeature & nCount_RNA > nCount & percent.mt_RNA < pt_mt)
  }
  names(filter_seurat_list) <- names(sample_seurat_list)

  return(filter_seurat_list)
}



#--------------------------#
## Add cell cycle score ----
#--------------------------#
# Description: Calculates cell cycle phase scores based on canonical markers 
#              using Seurat package.
# Input(s):    
#           sample_seurat_list: A list of Seurat objects to calculate cell cycle
#                               score on.
# Output(s): This function outputs a Seurat list with cell cycle scores added to
#            the metadata.
add_cell_cycle_score <- function(sample_seurat_list){
  # Add cell cycle score
  s.genes <- cc.genes.updated.2019$s.genes
  g2m.genes <- cc.genes.updated.2019$g2m.genes
  
  for (i in 1:length(sample_seurat_list)){
    message(names(sample_seurat_list)[i])
    DefaultAssay(sample_seurat_list[[i]]) <- "RNA"
    sample_seurat_list[[i]] <- SCTransform(sample_seurat_list[[i]], 
                                           vst.flavor = "v2", verbose = F)
    DefaultAssay(sample_seurat_list[[i]]) <- "SCT"
    sample_seurat_list[[i]] <- CellCycleScoring(sample_seurat_list[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident = F) 
    DefaultAssay(sample_seurat_list[[i]]) <- "RNA"
  }
  
  return(sample_seurat_list)
}


#--------------------#
## DoubletFinder ----
#--------------------#
# Description: Runs DoubletFinder, adding doublet information to each Seurat
#              object in the list.
# Input(s):
#         sample_seurat_list: A list of Seurat objects to run DoubletFinder on.
#         norm_method: Normalization method. Options are: "SCT" or "norm_data".
#                      "SCT" uses SCTransform() version 2 and "norm_data" uses 
#                      NormalizeData().
#         manual_dbr: Doublet rate. Default is set to 0.05 (5%).
# Output(s): 
#         sample_seurat_list: A list of Seurat objects with scDblFinder results.
run_doubletfinder <- function(sample_seurat_list, norm_method = "SCT", manual_dbr = 0.05) {
  message("Doublet detection with DoubletFinder")
  sapply(names(sample_seurat_list), function(XX) {
    if (norm_method == "SCT") {
      message("Normalization method: SCTransform V2")
      message("DoubletFinder for sample ", XX)
      message("Normalization, PCA, UMAP, and clustering")
      DefaultAssay(sample_seurat_list[[XX]]) <- "RNA"
      sample_seurat_list[[XX]] <- SCTransform(sample_seurat_list[[XX]], vst.flavor = "v2", 
                                              verbose = F) %>%
        RunPCA(verbose = F) %>%
        RunUMAP(dims = 1:10, verbose = F) %>%
        FindNeighbors(dims = 1:10, verbose = F) %>%
        FindClusters(resolution = 1, verbose = F)
      
      message("Find parameters")
      sweep.list_test <- paramSweep_v3(sample_seurat_list[[XX]], PCs = 1:10, sct = T) 
      sweep.stats_test <- summarizeSweep(sweep.list_test, GT = F)
      bcmvn <- find.pK(sweep.stats_test)
      
      message("Estimate homotypic doublet proportion")
      homo.prop <- modelHomotypic(sample_seurat_list[[XX]]@meta.data$seurat_clusters)
      
      message("Select nExp_poi assuming 5% doublet rate (from 10X manual)")
      nExp_poi <- round(manual_dbr*ncol(sample_seurat_list[[XX]]))
      nExp_poi.adj <- round(nExp_poi*(1-homo.prop))
      pk <- as.numeric(as.vector(bcmvn$pK)[which.max(bcmvn$BCmetric)])
      
      message("Run DoubletFinder")
      sample_seurat_list[[XX]] <- doubletFinder_v3(sample_seurat_list[[XX]], PCs = 1:10, pN = 0.25, 
                                                   pK = pk, nExp = nExp_poi, sct = TRUE, 
                                                   reuse.pANN = FALSE)
      message("Run DoubletFinder again w/ pANN")
      pANN_col_pos <- ncol(sample_seurat_list[[XX]]@meta.data) - 1
      sample_seurat_list[[XX]] <- doubletFinder_v3(
        sample_seurat_list[[XX]], PCs = 1:10, pN = 0.25, pK = pk, nExp = nExp_poi.adj, 
        sct = TRUE, reuse.pANN = colnames(sample_seurat_list[[XX]]@meta.data)[pANN_col_pos])
      
    } else if (norm_method == "norm_data") {
      message("Normalization method: NormalizeData")
      message("Normalization, PCA, UMAP, and clustering")
      DefaultAssay(sample_seurat_list[[XX]]) <- "RNA"
      sample_seurat_list[[XX]] <- NormalizeData(sample_seurat_list[[XX]]) %>%
        ScaleData() %>%
        FindVariableFeatures() %>%
        RunPCA(verbose = F) %>%
        RunUMAP(dims = 1:10, verbose = F) %>%
        FindNeighbors(dims = 1:10, verbose = F) %>%
        FindClusters(resolution = 1, verbose = F)
      
      message("Find parameters")
      sweep.list_test <- paramSweep_v3(sample_seurat_list[[XX]], PCs = 1:10, sct = F) 
      sweep.stats_test <- summarizeSweep(sweep.list_test, GT = F)
      bcmvn <- find.pK(sweep.stats_test)
      
      message("Estimate homotypic doublet proportion")
      homo.prop <- modelHomotypic(sample_seurat_list[[XX]]@meta.data$seurat_clusters)
      
      message("Select nExp_poi assuming 5% doublet rate (from 10X manual)")
      nExp_poi <- round(manual_dbr*ncol(sample_seurat_list[[XX]]))
      nExp_poi.adj <- round(nExp_poi*(1-homo.prop))
      pk <- as.numeric(as.vector(bcmvn$pK)[which.max(bcmvn$BCmetric)])
      
      message("Run DoubletFinder")
      sample_seurat_list[[XX]] <- doubletFinder_v3(sample_seurat_list[[XX]], PCs = 1:10, pN = 0.25, 
                                                   pK = pk, nExp = nExp_poi, sct = FALSE, 
                                                   reuse.pANN = FALSE) 
      message("Run DoubletFinder again w/ pANN")
      pANN_col_pos <- ncol(sample_seurat_list[[XX]]@meta.data) - 1
      sample_seurat_list[[XX]] <- doubletFinder_v3(
        sample_seurat_list[[XX]], PCs = 1:10, pN = 0.25, pK = pk, nExp = nExp_poi.adj, 
        sct = FALSE, reuse.pANN = colnames(sample_seurat_list[[XX]]@meta.data)[pANN_col_pos])
    }
    
    # Rename columns
    names(sample_seurat_list[[XX]]@meta.data)[pANN_col_pos] <- "pANN_col"
    names(sample_seurat_list[[XX]]@meta.data)[pANN_col_pos+1] <- "DF.class_col"
    names(sample_seurat_list[[XX]]@meta.data)[pANN_col_pos+2] <- "doublet_finder"
    
    message("Finished doublet detection")
    return(sample_seurat_list[[XX]])
    
  })
}


#------------------#
## scDblFinder ----
#------------------#
# Description: Runs scDblFinder, adding doublet information to each Seurat object
#              in the list.
# Input(s):
#         sample_seurat_list: A list of Seurat objects to run scDblFinder on.
#         manual_dbr: Doublet rate. Default is set to 0.05 (5%). If you want to 
#                     scDblFinder to auto-detect the doublet rate, set manual_dbr
#                     equal to NULL.
#         batch_col: Column name in "metadata" that contains batch ID. This 
#                    corresponds to each sequence sample to make a Seurat object 
#                    for.
# Output(s):
#         sample_seurat_list: A list of Seurat objects with scDblFinder results.
run_scdblfinder <- function(sample_seurat_list, manual_dbr = 0.05, batch_col = "Batch_ID"){
  # Merge Seurat list
  sample_seurat_obj <- merge(sample_seurat_list[[1]],
                             y = sample_seurat_list[2:length(sample_seurat_list)])
  
  # Create SingleCellExperiment object
  rawRNA_counts <- GetAssayData(sample_seurat_obj, assay = "RNA", slot = "counts")
  scdbl_sce <- SingleCellExperiment(assays = list(counts = rawRNA_counts),
                                    colData = sample_seurat_obj@meta.data)
  
  message("Doublet detection with scDblFinder")
  # Choose to set dbr to NULL (automatic) or 0.05
  if (is.null(manual_dbr)) {
    message("Automatically setting doublet rate")
    # Run scDblFinder
    scdbl_sce <- scDblFinder(scdbl_sce, samples = batch_col,
                             includePCs = 1:10) 
    # Extract names of doublet cells
    scdbl_cells <- rownames(scdbl_sce@colData[which(
      scdbl_sce@colData$scDblFinder.class == "doublet"), ])
    # Create new metadata column in Seurat object
    sample_seurat_obj$scdbl_auto <- ""
    sample_seurat_obj$scdbl_auto[colnames(sample_seurat_obj) %in% scdbl_cells] <- "Doublet"
    sample_seurat_obj$scdbl_auto[sample_seurat_obj$scdbl_auto != "Doublet"] <- "Singlet"
  } else {
    message(paste0("Using manually doublet rate (", manual_dbr, ")"))
    # Run scDblFinder
    scdbl_sce <- scDblFinder(scdbl_sce, samples = batch_col, 
                             includePCs = 1:10, dbr = manual_dbr) 
    # Extract names of doublet cells
    scdbl_cells <- rownames(scdbl_sce@colData[which(
      scdbl_sce@colData$scDblFinder.class == "doublet"), ]) 
    # Create new metadata column in Seurat object
    sample_seurat_obj$scdbl_manual <- ""
    sample_seurat_obj$scdbl_manual[colnames(sample_seurat_obj) %in% scdbl_cells] <- "Doublet"
    sample_seurat_obj$scdbl_manual[sample_seurat_obj$scdbl_manual != "Doublet"] <- "Singlet"
  }
  
  sample_seurat_list <- SplitObject(sample_seurat_obj, split.by = "Batch_ID")
  
  return(sample_seurat_list)
}

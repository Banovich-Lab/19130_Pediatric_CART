#==============================================================================#
# Author : Angela M. Oill, aoill@tgen.org
# Date: 2025/12/12
# Project: Pediatric CAR-T 
# Description: Integration of all CSF and PBMC samples
#==============================================================================#


#==============================================================================#
# Load Libraries, helper modules and functions ----
#==============================================================================#
#library(Matrix)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(dittoSeq)
library(scProportionTest)
library(scales)

source("/home/aoill/projects/SingleCellBestPractices/scripts/preprocessing_qc_module.R")
source("/home/aoill/projects/SingleCellBestPractices/scripts/helper_functions_module.R")


#==============================================================================#
# Set variables ----
#==============================================================================#
set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)


#==============================================================================#
# Load in PBMC and CSF objects ----
#==============================================================================#
## All - Lymphoid + Myeloid ----
csf_obj_lym <- readRDS("/scratch/aoill/projects/CAR-T/00_new_2025/geo_2026/csf_lymphoid_seurat_obj.rds")
pbmc_obj_lym <- readRDS("/scratch/aoill/projects/CAR-T/00_new_2025/geo_2026/pbmc_lymphoid_seurat_obj.rds")
pbmc_obj_lym <- JoinLayers(pbmc_obj_lym)

# add in myeloid
csf_obj_mye <- readRDS("/scratch/aoill/projects/CAR-T/00_new_2025/geo_2026/csf_myeloid_seurat_obj.rds")
csf_obj_mye <- JoinLayers(csf_obj_mye)
csf_obj_mye@meta.data$ct_final_all <- csf_obj_mye@meta.data$ct_final

pbmc_obj_mye <- readRDS("/scratch/aoill/projects/CAR-T/00_new_2025/geo_2026/pbmc_myeloid_seurat_obj.rds")
pbmc_obj_mye <- JoinLayers(pbmc_obj_mye)
pbmc_obj_mye@meta.data$ct_final_all <- pbmc_obj_mye@meta.data$ct_final


# Split object by how it will be integrated (old way of splitting)
csf_obj_lym_list <- SplitObject(csf_obj_lym, split.by = "UPN_Cycle_Day_Sample_Type_Batch")
csf_obj_mye_list <- SplitObject(csf_obj_mye, split.by = "UPN_Cycle_Day_Sample_Type_Batch")

pbmc_obj_lym_list <- SplitObject(pbmc_obj_lym, split.by = "UPN_Cycle_Day_Sample_Type_Batch")
pbmc_obj_mye_list <- SplitObject(pbmc_obj_mye, split.by = "UPN_Cycle_Day_Sample_Type_Batch")


# Add as one list of objects
csf_pbmc_obj_list <- c(csf_obj_lym_list, 
                           csf_obj_mye_list,
                           pbmc_obj_lym_list,
                           pbmc_obj_mye_list)
csf_pbmc_obj_merge <- merge(x = csf_pbmc_obj_list[[1]], y = csf_pbmc_obj_list[2:length(csf_pbmc_obj_list)])

# Join and resplit by UPN_Cycle_Day_Sample_Type_Batch
DefaultAssay(csf_pbmc_obj_merge) <- "SoupX_RNA"
csf_pbmc_obj_integrated_no_TRBV <- JoinLayers(csf_pbmc_obj_merge)

# Split layers by RUN
csf_pbmc_obj_integrated_no_TRBV[["SoupX_RNA"]] <- split(csf_pbmc_obj_integrated_no_TRBV[["SoupX_RNA"]], 
                                                        f = csf_pbmc_obj_integrated_no_TRBV$UPN_Cycle_Day_Sample_Type_Batch)
csf_pbmc_obj_integrated_no_TRBV
table(csf_pbmc_obj_integrated_no_TRBV$UPN_Cycle_Day_Sample_Type_Batch) # to get min number of cells in a sample
min(table(csf_pbmc_obj_integrated_no_TRBV$UPN_Cycle_Day_Sample_Type_Batch)) # to get min number of cells in a sample
#[1] 94

# SCTransform is performed per sample, since they are split into layers
csf_pbmc_obj_integrated_no_TRBV <- SCTransform(csf_pbmc_obj_integrated_no_TRBV, assay = "SoupX_RNA", 
                                               vst.flavor = "v2",
                                               vars.to.regress = c("S.Score", "G2M.Score"),
                                               return.only.var.genes = FALSE,
                                               variable.features.n = 2500)

# Get variable features minus TRBV genes
orig_var_feat_tcells <- csf_pbmc_obj_integrated_no_TRBV@assays[["SCT"]]@var.features
filt_var_feat_tcells <- grep("^TRBV", orig_var_feat_tcells, invert = TRUE, value = TRUE)

csf_pbmc_obj_integrated_no_TRBV <- RunPCA(csf_pbmc_obj_integrated_no_TRBV, features = filt_var_feat_tcells)
npcs <- min(get_pcs(csf_pbmc_obj_integrated_no_TRBV))
npcs # 14 (2500)
Sys.time()

# Integrate the layers using the SCT values
csf_pbmc_obj_integrated_no_TRBV <- IntegrateLayers(object = csf_pbmc_obj_integrated_no_TRBV, 
                                                   method = RPCAIntegration,
                                                   assay = "SCT", # either specify here or run default assay to SCT
                                                   orig.reduction = "pca", 
                                                   new.reduction = "integrated.rpca",
                                                   verbose = F,
                                                   normalization.method = "SCT", 
                                                   dims = 1:npcs,
                                                   features = filt_var_feat_tcells,
                                                   k.weight = 50
)

csf_pbmc_obj_integrated_no_TRBV <- FindNeighbors(csf_pbmc_obj_integrated_no_TRBV, dims = 1:npcs, verbose = FALSE, reduction = "integrated.rpca")
csf_pbmc_obj_integrated_no_TRBV <- RunUMAP(csf_pbmc_obj_integrated_no_TRBV, dims = 1:npcs, verbose = FALSE, 
                                           reduction = "integrated.rpca", reduction.name = "umap.integrated.rpca")
csf_pbmc_obj_integrated_no_TRBV <- FindClusters(csf_pbmc_obj_integrated_no_TRBV, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), verbose = FALSE)

DimPlot(csf_pbmc_obj_integrated_no_TRBV, reduction = "umap.integrated.rpca", group.by = "SCT_snn_res.0.2",
        raster= FALSE, label = T, #cols = randomcoloR::distinctColorPalette(20)
) + NoLegend()


# Add in cell type annotations 
csf_obj_lym_cell_types <- csf_obj_lym@meta.data %>% 
  dplyr::select(cellID, ct_final_all)

pbmc_obj_lym_cell_types <- pbmc_obj_lym@meta.data %>% 
  dplyr::select(cellID, ct_final_all)

all_cell_types <- rbind(csf_obj_lym_cell_types, pbmc_obj_lym_cell_types)


csf_pbmc_obj_integrated_no_TRBV_meta <- csf_pbmc_obj_integrated_no_TRBV@meta.data %>%
  dplyr::select(cellID)

csf_pbmc_obj_integrated_no_TRBV_meta_l_join <- left_join(
  csf_pbmc_obj_integrated_no_TRBV_meta, 
  all_cell_types)


table(csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all)
table(all_cell_types$ct_final_all)

ct_colors <- c(
  "T (Memory)" = "#88CCEE", "T (GAPDH+ Glycolysis)" = "#332288", 
  "T (Effector-memory)" = "#44AA99", "T (Undifferentiated)" = "#117733",
  "T (Cytotoxic)" = "#AA4499", 
  "T (Heat-shock)" = "#CC6677", "T (Proliferating)" = "#DDCC77", 
  "Treg" = "#882255",
  "T (Effector)" = "#75871B", "T (Naïve)" = "#DAB2FF", #"#E88C23"
  "T (Central memory)" = "#FF637E", #"#A25933",
  "T (MAIT)" = "#3A4330",
  
  "T cells" = "#AA4499", 
  
  
  "NK" = "#0000ff", 
  "NK (CD56-bright, CD16+)" = "#0000ff",  
  "NK (CD56-dim, CD16+)" = "#6767D6", "B" = "#994F00", "pDC" = "#d55e00",
  
  "Classical Monocyte" = "#bcd095", "Mo/Mac" = "#6dce52", "Non-classical Monocyte" = "#697f3f",
  "S100 high, HLA class II low Monocyte" = "#64e9b3",
  "Intermediate Monocyte" = "#81e858",
  "Monocyte" = "#bcd095",
  
  "Macrophage" = "#64A373", 
  "cDC1" = "#927754", "cDC2" = "#e3bb87", 
  "migDC" = "#de9e35"
)


DimPlot(csf_pbmc_obj_integrated_no_TRBV, 
        reduction = "umap.integrated.rpca", group.by = "ct_final_all",
        raster= FALSE, label = T
) +
  scale_color_manual(values = ct_colors) +
  NoLegend()


DimPlot(csf_pbmc_obj_integrated_no_TRBV, 
        reduction = "umap.integrated.rpca", group.by = "ct_final_all",
        raster= FALSE, label = T,
        split.by = "Sample_Type"
) +
  scale_color_manual(values = ct_colors) +
  NoLegend()

# Change colors to simplify annotations (so T cells merged and NK merged)
unique(csf_pbmc_obj_integrated_no_TRBV$ct_final_all)

csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all_simplified <- csf_pbmc_obj_integrated_no_TRBV$ct_final_all # or any other initialization value
csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all_simplified[csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all == "T (Memory)"] <- "T cells" 
csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all_simplified[csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all == "T (Undifferentiated)"] <- "T cells" 
csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all_simplified[csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all == "T (Effector-memory)"] <- "T cells" 
csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all_simplified[csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all == "T (Cytotoxic)"] <- "T cells" 
csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all_simplified[csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all == "Treg"] <- "T cells" 
csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all_simplified[csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all == "T (GAPDH+ Glycolysis)"] <- "T cells" 
csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all_simplified[csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all == "T (Heat-shock)"] <- "T cells" 
csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all_simplified[csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all == "T (Proliferating)"] <- "T cells" 
csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all_simplified[csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all == "T (Effector)"] <- "T cells" 
csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all_simplified[csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all == "T (Naïve)"] <- "T cells" 
csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all_simplified[csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all == "T (MAIT)"] <- "T cells" 
csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all_simplified[csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all == "T (Central memory)"] <- "T cells" 

csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all_simplified[csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all == "NK (CD56-dim, CD16+)"] <- "NK" 
csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all_simplified[csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all == "NK (CD56-bright, CD16+)"] <- "NK" 


csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all_simplified[csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all == "Classical Monocyte"] <- "Monocyte" 
csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all_simplified[csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all == "Non-classical Monocyte"] <- "Monocyte" 
csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all_simplified[csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all == "S100 high, HLA class II low Monocyte"] <- "Monocyte" 
csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all_simplified[csf_pbmc_obj_integrated_no_TRBV@meta.data$ct_final_all == "Intermediate Monocyte"] <- "Monocyte" 


unique(csf_pbmc_obj_integrated_no_TRBV$ct_final_all_simplified)


#-----------------------------------------------------------------------------#
## SAVE ----
#-----------------------------------------------------------------------------#
saveRDS(csf_pbmc_obj_integrated_no_TRBV, "/scratch/aoill/projects/CAR-T/00_new_2025/geo_2026/csf_and_pbmc_integrated_together_seurat_obj.rds")


#==============================================================================#
# Run linear regression ----
#==============================================================================#
# Keep only common cell types (ones found in both CSF and PBMC)
obj <- csf_pbmc_obj_integrated_no_TRBV

DefaultAssay(obj) <- "SoupX_RNA"
obj <- JoinLayers(obj)

all_CT_results <- list() # this is where the lm test results will go. It will be
# a list of list for the results for each cell type (each cell type is in its 
# own list)
levels(as.factor(obj@meta.data$ct_final_all_simplified))
#[1] "B"          "cDC1"       "cDC2"       "Macrophage" "migDC"      "Mo/Mac"     "Monocyte"   "NK"         "pDC"        "T cells"   

cts_to_analyze <- levels(as.factor(obj@meta.data$ct_final_all_simplified))
# Only keep cell types found in both CSF and PBMC


#ct = "T cells"
#ct = "Mo/Mac"
# "pDCs", "Plasma" do not have enough samples (2 each)
for (ct in cts_to_analyze) {
  #for (ct in c("B", "Monocyte", "NK", "T cells")) {
  print(ct)
  
  # I think I need to separate out the PBMCs and CSFs into two different object/matricies
  # because when I do the test it will be between these two sample types
  # check number of cells for these
  ncells_ct_pbmcs <- obj@meta.data %>% filter(Sample_Type == "PBMC" & 
                                                ct_final_all_simplified == ct & 
                                                UPN %in% c("514", "515", "574", "689", "692", "705", "716")
  ) %>%
    nrow()
  ncells_ct_csf <- obj@meta.data %>% filter(Sample_Type == "CSF" & 
                                              ct_final_all_simplified == ct & 
                                              UPN %in% c("514", "515", "574", "689", "692", "705", "716")
  ) %>%
    nrow()
  
  if (ncells_ct_pbmcs == 0 | ncells_ct_csf == 0) {
    print(paste(ct, " not observed in either PBMCs or CSF", sep = ""))
  } else {
    obj_CT_pbmcs <- subset(obj, subset = (Sample_Type == "PBMC" & 
                                            ct_final_all_simplified == ct & 
                                            UPN %in% c("514", "515", "574", "689", "692", "705", "716"))) # remove 625 and 626
    obj_CT_csf <- subset(obj, subset = (Sample_Type == "CSF" & ct_final_all_simplified == ct & 
                                          UPN %in% c("514", "515", "574", "689", "692", "705", "716")))
    
    
    # new
    # Extract the normalized expression data (data is log-transformed counts)
    #obj_CT_pbmcs_expression_int <- as.data.frame(as.matrix(obj_CT_pbmcs@assays$SoupX_RNA@layers$data))
    obj_CT_pbmcs_expression_int <- as.matrix(GetAssayData(object = obj_CT_pbmcs, assay = "SoupX_RNA", slot = "data")) # rows are genes, columns are cells
    obj_CT_pbmcs_expression <- t(obj_CT_pbmcs_expression_int) # columns are genes, rows are cells
    obj_CT_pbmcs_expression <- as.data.frame(obj_CT_pbmcs_expression)
    # Add sample information to the data
    obj_CT_pbmcs_expression$UPN_Cycle <- obj_CT_pbmcs@meta.data$UPN_Cycle
    # cells are rows, cols genes
    
     
    # new
    # Extract the normalized expression data (data is log-transformed counts)
    #obj_CT_csf_expression_int <- as.data.frame(as.matrix(obj_CT_pbmcs@assays$SoupX_RNA@layers$data))
    obj_CT_csf_expression_int <- as.matrix(GetAssayData(object = obj_CT_csf, assay = "SoupX_RNA", slot = "data"))
    obj_CT_csf_expression <- t(obj_CT_csf_expression_int)
    obj_CT_csf_expression <- as.data.frame(obj_CT_csf_expression)
    # Add sample information to the data
    obj_CT_csf_expression$UPN_Cycle <- obj_CT_csf@meta.data$UPN_Cycle
    
    
    # new, at least 10 cells per UPN_Cycle for a given cell type
    keep_pbmcs <- which(table(obj_CT_pbmcs_expression$UPN_Cycle)>10)
    keep_csf <- which(table(obj_CT_csf_expression$UPN_Cycle)>10)
    
    keep <- intersect(names(keep_pbmcs), names(keep_csf))
    # Here count the number of samples if less than 3 skip (print warning), else 
    # continue
    if (length(keep) < 3){
      print(paste("Not enough sample for this cell type. ", length(keep), " total samples for ", ct, ".", sep = ""))
    } else {
      
      
      # NEW
      #obj_CT_pbmcs_expression$UPN_Cycle
      #rownames(obj_CT_pbmcs_expression)
      # filter to keep samples
      # transform rows and cols
      # make sure rownames are retined 
      
      
      obj_CT_pbmcs_expression_filtered <- obj_CT_pbmcs_expression %>% filter(UPN_Cycle %in% keep)
      obj_CT_csf_expression_filtered <- obj_CT_csf_expression %>% filter(UPN_Cycle %in% keep)
      
      
      obj_CT_pbmcs_expression_filtered_mean <- obj_CT_pbmcs_expression_filtered %>%
        dplyr::group_by(UPN_Cycle) %>%
        dplyr::summarise(across(everything(), mean, na.rm = TRUE))
      obj_CT_pbmcs_expression_filtered_mean <- as.data.frame(obj_CT_pbmcs_expression_filtered_mean)
      rownames(obj_CT_pbmcs_expression_filtered_mean) <- obj_CT_pbmcs_expression_filtered_mean$UPN_Cycle
      obj_CT_pbmcs_expression_filtered_mean$UPN_Cycle <- NULL
      
      
      obj_CT_csf_expression_filtered_mean <- obj_CT_csf_expression_filtered %>%
        dplyr::group_by(UPN_Cycle) %>%
        dplyr::summarise(across(everything(), mean, na.rm = TRUE))
      obj_CT_csf_expression_filtered_mean <- as.data.frame(obj_CT_csf_expression_filtered_mean)
      rownames(obj_CT_csf_expression_filtered_mean) <- obj_CT_csf_expression_filtered_mean$UPN_Cycle      
      obj_CT_csf_expression_filtered_mean$UPN_Cycle <- NULL
      
      # I think rows need to be genes so need to transform data
      obj_CT_pbmcs_expression_filtered_mean_t <- t(obj_CT_pbmcs_expression_filtered_mean)
      obj_CT_csf_expression_filtered_mean_t <- t(obj_CT_csf_expression_filtered_mean)
      
      
      
      # We want to analyze genes that are actually expressed
      # Keep genes where at least half of the samples have expression >0
      # remove genes with at least half the number of samples with 0 expression
      smpl_num_half <- length(keep)/2
      
      obj_CT_pbmcs_expression_filtered_mean_t_filtered <- obj_CT_pbmcs_expression_filtered_mean_t[rowSums(obj_CT_pbmcs_expression_filtered_mean_t > 0) >= smpl_num_half, ]
      obj_CT_csf_expression_filtered_mean_t_filtered <- obj_CT_csf_expression_filtered_mean_t[rowSums(obj_CT_csf_expression_filtered_mean_t > 0) >= smpl_num_half, ]
      
      # only analyze the data between the two
      ngenes <- length(intersect(rownames(obj_CT_pbmcs_expression_filtered_mean_t_filtered), rownames(obj_CT_csf_expression_filtered_mean_t_filtered)))
      print(paste(ct, ". Total number of genes that pass low expression filter: ", ngenes, sep = ""))
      
      
      # Sort columns in ascending alphabetical order. Want both PBMCs and CSF data 
      # frames to be in the same order so the test is comparing each matched UPN/Cycle
      # to eachother
      new_order_pbmcs <- sort(colnames(obj_CT_pbmcs_expression_filtered_mean_t_filtered))
      new_order_csf <- sort(colnames(obj_CT_csf_expression_filtered_mean_t_filtered))
      
      obj_CT_pbmcs_expression_filtered_mean_t_filtered_order <- obj_CT_pbmcs_expression_filtered_mean_t_filtered[, new_order_pbmcs]
      obj_CT_pbmcs_expression_filtered_mean_t_filtered_order <- as.data.frame(obj_CT_pbmcs_expression_filtered_mean_t_filtered_order)
      #head(obj_CT_pbmcs_expression_filtered_mean_t_filtered_order)
      
      obj_CT_csf_expression_filtered_mean_t_filtered_order <- obj_CT_csf_expression_filtered_mean_t_filtered[, new_order_csf]
      obj_CT_csf_expression_filtered_mean_t_filtered_order <- as.data.frame(obj_CT_csf_expression_filtered_mean_t_filtered_order)
      #head(obj_CT_csf_expression_filtered_mean_t_filtered_order)
      
      # Now that the gene expression matrices are in the desired format, I can now
      # run each test per gene and output the results for all genes
      
      # Get list of genes that overlap between pbmcs and csf (should be the same but
      # running for good measure)
      # make a merged gene list
      genes_pbmcs <- rownames(obj_CT_pbmcs_expression_filtered_mean_t_filtered_order)
      genes_csf <- rownames(obj_CT_csf_expression_filtered_mean_t_filtered_order)
      
      genes_intersect <- (intersect(genes_pbmcs, genes_csf))
      
      
      # Add gene names to data frames so genes can be extracted
      obj_CT_pbmcs_expression_filtered_mean_t_filtered_order$genes <- rownames(obj_CT_pbmcs_expression_filtered_mean_t_filtered_order)
      obj_CT_csf_expression_filtered_mean_t_filtered_order$genes <- rownames(obj_CT_csf_expression_filtered_mean_t_filtered_order)
      
      
      # Loop through each gene and run test
      # establish data frame
      lm_results <- c()
      
      for (gene in 1:length(genes_intersect)) {
        #for (gene in 1:400) {
        genei <- genes_intersect[gene]
        #print(gene)
        #print(genei)
        
        # get expression data for the gene
        t_obj_CT_pbmcs_expression_sum_df_order_gene <- (obj_CT_pbmcs_expression_filtered_mean_t_filtered_order %>% filter(genes == genei))
        t_obj_CT_pbmcs_expression_sum_df_order_gene$genes <- NULL
        t_obj_CT_pbmcs_expression_sum_df_order_gene <- t(t_obj_CT_pbmcs_expression_sum_df_order_gene)
        
        t_obj_CT_csf_expression_sum_df_order_gene <- (obj_CT_csf_expression_filtered_mean_t_filtered_order %>% filter(genes == genei))
        t_obj_CT_csf_expression_sum_df_order_gene$genes <- NULL
        t_obj_CT_csf_expression_sum_df_order_gene <- t(t_obj_CT_csf_expression_sum_df_order_gene)
        
        # prep data frames
        colnames(t_obj_CT_pbmcs_expression_sum_df_order_gene) <- c("PBMCs")
        colnames(t_obj_CT_csf_expression_sum_df_order_gene) <- c("CSF")
        merge_LM_dat <- merge(t_obj_CT_pbmcs_expression_sum_df_order_gene, t_obj_CT_csf_expression_sum_df_order_gene, 
                              by = "row.names")
        
        n_sample_pbmcs <- merge_LM_dat %>% filter(PBMCs > 0) %>% nrow()
        n_sample_csf <- merge_LM_dat %>% filter(CSF > 0) %>% nrow()
        
        # run lm
        merge_LM_result <- lm(CSF~PBMCs, data = merge_LM_dat)
        
        cor_test_test <- cor.test(merge_LM_dat$CSF, merge_LM_dat$PBMCs)
        
        # extract results
        #pval_result <- summary(merge_LM_result)$coefficients[2, 4] # p value
        pval_result <- tryCatch(
          #try to do this
          {
            #some expression
            summary(merge_LM_result)$coefficients[2, 4]
          },
          #if an error occurs, tell me the error
          error=function(e) {
            message('An Error Occurred')
            #print(e)
            return(NA)
          },
          #if a warning occurs, tell me the warning
          warning=function(w) {
            message('A Warning Occurred')
            print(w)
            return(summary(merge_LM_result)$coefficients[2, 4])
          }
        )
        
        cor_result <- cor_test_test$estimate # correlation coefficient
        
        # add results to df
        # add number of samples used in lm
        n_samples <- nrow(t_obj_CT_pbmcs_expression_sum_df_order_gene)
        row_to_add <- cbind(ct, genei, n_samples, n_sample_pbmcs, n_sample_csf, cor_result, pval_result)
        lm_results <- rbind(lm_results, row_to_add)
      }
      
      
      lm_results_df <- as.data.frame(lm_results)
      
      lm_results_df_clean <- na.omit(lm_results_df) 
      
      # Remove AC genes and LINCs
      lm_results_df_clean_gene_filter <- lm_results_df_clean %>% 
        filter(grepl("\\.", genei) == FALSE, grepl("LINC", genei) == FALSE)
      
      # Correct for multiple testing
      lm_results_df_clean_gene_filter$padj <- p.adjust(lm_results_df_clean_gene_filter$pval_result, method = "fdr")
      
      # Add results to a list of list
      all_CT_results <- append(all_CT_results, list(lm_results_df_clean_gene_filter))
    }
    
  }
  
}

names(all_CT_results) <- c("B", "cDC2", "Monocyte", "NK", "pDC", "T cells")


## SAVE ## 
saveRDS(all_CT_results, "/scratch/aoill/projects/CAR-T/00_new_2025/rds_files/All_CSF_PBMC_correlations_2026_01_13_exp_fix_02.rds")
#all_CT_results <- readRDS("/scratch/aoill/projects/CAR-T/00_new_2025/rds_files/All_CSF_PBMC_correlations_2026_01_13_exp_fix_02.rds")

all_CT_results_one_df <- rbind(all_CT_results[["B"]], all_CT_results[["cDC2"]],
                               all_CT_results[["Monocyte"]], all_CT_results[["NK"]],
                               all_CT_results[["pDC"]], all_CT_results[["T cells"]])
all_CT_results_one_df$padj2 <- p.adjust(all_CT_results_one_df$pval_result, method = "fdr")


# Get # sig and # sig after FDR
# B cells
all_CT_results[["B"]] %>% 
  nrow()
# 9543
all_CT_results[["B"]] %>% 
  filter(pval_result < 0.05) %>%
  nrow()
# 343

all_CT_results[["B"]] %>% 
  filter(padj < 0.05) %>%
  nrow()
# 35

# cDCs
all_CT_results[["cDC2"]] %>% 
  nrow()
# 10828
all_CT_results[["cDC2"]] %>% 
  filter(pval_result < 0.05) %>%
  nrow()
# 1513
all_CT_results[["cDC2"]] %>% 
  filter(padj < 0.05) %>%
  nrow()
# 149

# Monocytes
all_CT_results[["Monocyte"]] %>% 
  nrow()
# 13123
all_CT_results[["Monocyte"]] %>% 
  filter(pval_result < 0.05) %>%
  nrow()
# 2842
all_CT_results[["Monocyte"]] %>% 
  filter(padj < 0.05) %>%
  nrow()
# 984
#all_CT_results_one_df %>% 
#  filter(padj2 < 0.05) %>%
#  nrow()

# NK
all_CT_results[["NK"]] %>% 
  nrow()
# 11368
all_CT_results[["NK"]] %>% 
  filter(pval_result < 0.05) %>%
  nrow()
# 2357
all_CT_results[["NK"]] %>% 
  filter(padj < 0.05) %>%
  nrow()
# 685


# pDC
all_CT_results[["pDC"]] %>% 
  nrow()
# 9561
all_CT_results[["pDC"]] %>% 
  filter(pval_result < 0.05) %>%
  nrow()
# 479
all_CT_results[["pDC"]] %>% 
  filter(padj < 0.05) %>%
  nrow()
# 0 

# T cells
all_CT_results[["T cells"]] %>% 
  nrow()
# 13164
all_CT_results[["T cells"]] %>% 
  filter(pval_result < 0.05) %>%
  nrow()
# 4119
all_CT_results[["T cells"]] %>% 
  filter(padj < 0.05) %>%
  nrow()
# 2408

# Make output results for table s7 ----
all_CT_results_one_df <- rbind(
  all_CT_results[["B"]],
  all_CT_results[["cDC2"]],
  all_CT_results[["Monocyte"]],
  all_CT_results[["NK"]],
  all_CT_results[["pDC"]],
  all_CT_results[["T cells"]]
)
rownames(all_CT_results_one_df) <- NULL
colnames(all_CT_results_one_df) <- c("Cell Type", "Gene", "Sample Number", "PBMC_smpls_expressed", "CSF_smpls_expressed",
                                     "R2", "p-val", "FDR")

write.csv(all_CT_results_one_df, "/home/aoill/tables/csf_pbmc_correlations_2026_01_13_exp_fix_02.csv",
          quote = F, row.names = F)


# get number of sig
all_CT_results_one_df %>% nrow()
# 67587
all_CT_results_one_df %>% filter(FDR <=0.05) %>% nrow()
# 4261
# 4261/67587 = 6.3%


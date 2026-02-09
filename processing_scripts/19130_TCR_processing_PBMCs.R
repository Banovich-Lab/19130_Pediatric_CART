#==============================================================================#
# PBMC ----
#==============================================================================#

#------------------------------------------------------------------------------#
## Add TCR data to T cell object ----
#------------------------------------------------------------------------------#

### Add TCR data determined from scRepertoire ----

#### Load the contig data ----

# TCR info for each sample

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

## NEW SAMPLES (batches 7-10, 1)
# Batch 1
F07949 = read.csv("/tgen_labs/banovich/SingleCell/CellRanger/8_0_1/Projects/PEDCAR/VDJ/PEDCAR_0001_1_CS_Whole_C1_X5TCR_F07949_22YCMWLT3/outs/filtered_contig_annotations.csv")
# Batches 7-10
F07955 = read.csv("/tgen_labs/banovich/SingleCell/CellRanger/8_0_1/Projects/PEDCAR/VDJ/PEDCAR_0013_1_PB_Whole_C1_X5TCR_F07955_22YCMWLT3/outs/filtered_contig_annotations.csv")
F07956 = read.csv("/tgen_labs/banovich/SingleCell/CellRanger/8_0_1/Projects/PEDCAR/VDJ/PEDCAR_0014_1_PB_Whole_C1_X5TCR_F07956_CONCAT/outs/filtered_contig_annotations.csv")
F07957 = read.csv("/tgen_labs/banovich/SingleCell/CellRanger/8_0_1/Projects/PEDCAR/VDJ/PEDCAR_0015_1_PB_Whole_C1_X5TCR_F07957_22YCMWLT3/outs/filtered_contig_annotations.csv")
F07958 = read.csv("/tgen_labs/banovich/SingleCell/CellRanger/8_0_1/Projects/PEDCAR/VDJ/PEDCAR_0016_1_PB_Whole_C1_X5TCR_F07958_22YCMWLT3/outs/filtered_contig_annotations.csv")


#### Process the contig data ----
seurat_obj <- readRDS("/scratch/aoill/projects/CAR-T/00_new_2025/rds_files/pbmc_DoubletFinder_integrated_2025_09_24.rds")

metapath <- "/home/aoill/projects/CAR-T/00_new_2025/19130_all_batches_metadata_new_edited.csv"
sheet_name <- "Sheet1"
sample_metadata <- load_metadata(metadata_path = metapath, data_type = "csv", sheet_name = sheet_name)

seurat_obj_meta <- seurat_obj@meta.data
seurat_obj_meta$TCRseq_CellRanger_path <- NULL
seurat_obj_meta$TCRseq_ID_long <- NULL
seurat_obj_meta$TCRseq_ID <- NULL
seurat_obj_meta_join <- left_join(seurat_obj_meta, sample_metadata)

# add back to metadata in object
seurat_obj@meta.data$TCRseq_CellRanger_path <- seurat_obj_meta_join$TCRseq_CellRanger_path
seurat_obj@meta.data$TCRseq_ID_long <- seurat_obj_meta_join$TCRseq_ID_long
seurat_obj@meta.data$TCRseq_ID <- seurat_obj_meta_join$TCRseq_ID


seurat_list <- SplitObject(seurat_obj, split.by = "Cell_Prefix")

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
# [1] "Batch37" "Batch38" "Batch39" "Batch40" "Batch7"  "Batch8"  "Batch9"  "Batch10" "Batch1" 

# Merge all contigs and contig lists into one contig list
# Batched data
# Get the contig list for the multiplexed batches - original batches
contig_list_37 <- createHTOContigList(F05037, seurat_list_cells_renamed[["Batch37"]], group.by = c("TCRseq_ID", "UPN", "Cycle", "Day", "Sample_Type"))
contig_list_38 <- createHTOContigList(F05041, seurat_list_cells_renamed[["Batch38"]], group.by = c("TCRseq_ID", "UPN", "Cycle", "Day", "Sample_Type"))
contig_list_39 <- createHTOContigList(F05625, seurat_list_cells_renamed[["Batch39"]], group.by = c("TCRseq_ID", "UPN", "Cycle", "Day", "Sample_Type"))
contig_list_40 <- createHTOContigList(F05626, seurat_list_cells_renamed[["Batch40"]], group.by = c("TCRseq_ID", "UPN", "Cycle", "Day", "Sample_Type"))

# Get the contig list for the multiplexed batches - new batches
contig_list_7 <- createHTOContigList(F07955, seurat_list_cells_renamed[["Batch7"]], group.by = c("TCRseq_ID", "UPN", "Cycle", "Day", "Sample_Type"))
contig_list_8 <- createHTOContigList(F07956, seurat_list_cells_renamed[["Batch8"]], group.by = c("TCRseq_ID", "UPN", "Cycle", "Day", "Sample_Type"))
contig_list_9 <- createHTOContigList(F07957, seurat_list_cells_renamed[["Batch9"]], group.by = c("TCRseq_ID", "UPN", "Cycle", "Day", "Sample_Type"))
contig_list_10 <- createHTOContigList(F07958, seurat_list_cells_renamed[["Batch10"]], group.by = c("TCRseq_ID", "UPN", "Cycle", "Day", "Sample_Type"))
contig_list_1 <- createHTOContigList(F07949, seurat_list_cells_renamed[["Batch1"]], group.by = c("TCRseq_ID", "UPN", "Cycle", "Day", "Sample_Type"))


# Merge all contigs and contig lists into one contig list
contig_list_pbmc <- append(contig_list_37, contig_list_38)
contig_list_pbmc <- append(contig_list_pbmc, contig_list_39)
contig_list_pbmc <- append(contig_list_pbmc, contig_list_40)
contig_list_pbmc <- append(contig_list_pbmc, contig_list_7)
contig_list_pbmc <- append(contig_list_pbmc, contig_list_8)
contig_list_pbmc <- append(contig_list_pbmc, contig_list_9)
contig_list_pbmc <- append(contig_list_pbmc, contig_list_10)
contig_list_pbmc <- append(contig_list_pbmc, contig_list_1)

names(contig_list_pbmc)
#[1] "F05037.514.8.1.PBMC"  "F05037.514.5.1.PBMC"  "F05037.574.3.1.PBMC"  "F05037.515.8.1.PBMC"  "F05041.514.11.1.PBMC" "F05041.515. 4.2.PBMC"
#[7] "F05041.574. 8.1.PBMC" "F05625.626.1.0.PBMC"  "F05625.625.3.1.PBMC"  "F05625.626.4.1.PBMC"  "F05626.626.7.1.PBMC"  "F05626.625.4.1.PBMC" 
#[13] "F05626.625.1.0.PBMC"  "F07955.514.2.1.PBMC"  "F07955.705.5.1.PBMC"  "F07955.689.1.1.PBMC"  "F07955.716.8.1.PBMC"  "F07956.574.1.1.PBMC" 
#[19] "F07956.692.2.1.PBMC"  "F07956.716.4.1.PBMC"  "F07957.705.8.1.PBMC"  "F07957.689.4.1.PBMC"  "F07957.716.2.1.PBMC"  "F07958.689.12.1.PBMC"
#[25] "F07958.692. 4.1.PBMC" "F07958.514.13.1.PBMC" "F07949.574. 3.1.PBMC" "F07949.689.17.1.PBMC" "F07949.705. 1.1.PBMC"



names(contig_list_pbmc) <- c("F05037_514_8_1_PBMC",  "F05037_514_5_1_PBMC",  
                            "F05037_574_3_1_PBMC",  "F05037_515_8_1_PBMC",  
                            "F05041_514_11_1_PBMC", "F05041_515_4_2_PBMC", 
                            "F05041_574_8_1_PBMC", "F05625_626_1_0_PBMC",  
                            "F05625_625_3_1_PBMC", "F05625_626_4_1_PBMC",  
                            "F05626_626_7_1_PBMC",  "F05626_625_4_1_PBMC", 
                            "F05626_625_1_0_PBMC",  "F07955_514_2_1_PBMC",  
                            "F07955_705_5_1_PBMC",  "F07955_689_1_1_PBMC",  
                            "F07955_716_8_1_PBMC", "F07956_574_1_1_PBMC", 
                            "F07956_692_2_1_PBMC",  "F07956_716_4_1_PBMC",  
                            "F07957_705_8_1_PBMC",  "F07957_689_4_1_PBMC",  
                            "F07957_716_2_1_PBMC",  "F07958_689_12_1_PBMC", 
                            "F07958_692_4_1_PBMC", "F07958_514_13_1_PBMC", 
                            "F07949_574_3_1_PBMC", "F07949_689_17_1_PBMC", 
                            "F07949_705_1_1_PBMC"
)


#### Combine the contigs ----
# Will want to use the parameter removeNA set to true to remove any cell barcode 
# with an NA value in at least one of the chains
# https://github.com/BorchLab/scRepertoire/issues/373
combined <- combineTCR(contig_list_pbmc, 
                       samples =  c("F05037_514_8_1_PBMC",  "F05037_514_5_1_PBMC",  
                                    "F05037_574_3_1_PBMC",  "F05037_515_8_1_PBMC",  
                                    "F05041_514_11_1_PBMC", "F05041_515_4_2_PBMC", 
                                    "F05041_574_8_1_PBMC", "F05625_626_1_0_PBMC",  
                                    "F05625_625_3_1_PBMC", "F05625_626_4_1_PBMC",  
                                    "F05626_626_7_1_PBMC",  "F05626_625_4_1_PBMC", 
                                    "F05626_625_1_0_PBMC",  "F07955_514_2_1_PBMC",  
                                    "F07955_705_5_1_PBMC",  "F07955_689_1_1_PBMC",  
                                    "F07955_716_8_1_PBMC", "F07956_574_1_1_PBMC", 
                                    "F07956_692_2_1_PBMC",  "F07956_716_4_1_PBMC",  
                                    "F07957_705_8_1_PBMC",  "F07957_689_4_1_PBMC",  
                                    "F07957_716_2_1_PBMC",  "F07958_689_12_1_PBMC", 
                                    "F07958_692_4_1_PBMC", "F07958_514_13_1_PBMC", 
                                    "F07949_574_3_1_PBMC", "F07949_689_17_1_PBMC", 
                                    "F07949_705_1_1_PBMC"),
                       # cells ="T-AB", # depreciated argument 
                       filterMulti = FALSE, removeNA = TRUE
)

View(combined)


#### Add TCR info to the Seurat objects ----

### Add a new barcode column in the TCR list to match the Seurat objects
# Replace the barcode column in combined which matching barcode to seurat object 
# Path to metadata for your samples
metapath <- "/home/aoill/projects/CAR-T/00_new_2025/19130_all_batches_metadata_new_edited.csv"
sheet_name <- "Sheet1"
sample_metadata <- load_metadata(metadata_path = metapath, data_type = "csv", sheet_name = sheet_name)

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

names(combined2) <- c("F05037_514_8_1_PBMC",  "F05037_514_5_1_PBMC",  
                      "F05037_574_3_1_PBMC",  "F05037_515_8_1_PBMC",  
                      "F05041_514_11_1_PBMC", "F05041_515_4_2_PBMC", 
                      "F05041_574_8_1_PBMC", "F05625_626_1_0_PBMC",  
                      "F05625_625_3_1_PBMC", "F05625_626_4_1_PBMC",  
                      "F05626_626_7_1_PBMC",  "F05626_625_4_1_PBMC", 
                      "F05626_625_1_0_PBMC",  "F07955_514_2_1_PBMC",  
                      "F07955_705_5_1_PBMC",  "F07955_689_1_1_PBMC",  
                      "F07955_716_8_1_PBMC", "F07956_574_1_1_PBMC", 
                      "F07956_692_2_1_PBMC",  "F07956_716_4_1_PBMC",  
                      "F07957_705_8_1_PBMC",  "F07957_689_4_1_PBMC",  
                      "F07957_716_2_1_PBMC",  "F07958_689_12_1_PBMC", 
                      "F07958_692_4_1_PBMC", "F07958_514_13_1_PBMC", 
                      "F07949_574_3_1_PBMC", "F07949_689_17_1_PBMC", 
                      "F07949_705_1_1_PBMC"
)


## SAVE ##
saveRDS(combined2, "/scratch/aoill/projects/CAR-T/00_new_2025/geo_2026/pbmc_scRepertiore_output.rds")


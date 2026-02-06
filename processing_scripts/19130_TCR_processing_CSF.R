#==============================================================================#
# Author : Angela M. Oill, aoill@tgen.org
# Date: 2025/05/15
# Project: Pediatric CAR-T 
# Description: TCR processing
#==============================================================================#


#==============================================================================#
# Load Libraries, helper modules and functions ----
#==============================================================================#
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(scRepertoire)
library(RColorBrewer)
library(stringr)
source("/home/aoill/projects/SingleCellBestPractices/scripts/preprocessing_qc_module.R")
source("/home/aoill/projects/SingleCellBestPractices/scripts/helper_functions_module.R")


#==============================================================================#
# CSF ----
#==============================================================================#

#------------------------------------------------------------------------------#
## Add TCR data to T cell object ----
#------------------------------------------------------------------------------#

### Add TCR data determined from scRepertoire ----

#### Load the contig data ----

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

## NEW SAMPLES (batches 2-6, 11)
# Batch 11
F07959 = read.csv("/tgen_labs/banovich/SingleCell/CellRanger/8_0_1/Projects/PEDCAR/VDJ/PEDCAR_0017_1_PB_Whole_C1_X5TCR_F07959_22YCMWLT3/outs/filtered_contig_annotations.csv")
# Batches 2-6
F07950 = read.csv("/tgen_labs/banovich/SingleCell/CellRanger/8_0_1/Projects/PEDCAR/VDJ/PEDCAR_0008_1_CS_Whole_C1_X5TCR_F07950_22YCMWLT3/outs/filtered_contig_annotations.csv")
F07951 = read.csv("/tgen_labs/banovich/SingleCell/CellRanger/8_0_1/Projects/PEDCAR/VDJ/PEDCAR_0009_1_CS_Whole_C1_X5TCR_F07951_22YCMWLT3/outs/filtered_contig_annotations.csv")
F07952 = read.csv("/tgen_labs/banovich/SingleCell/CellRanger/8_0_1/Projects/PEDCAR/VDJ/PEDCAR_0010_1_CS_Whole_C1_X5TCR_F07952_22YCMWLT3/outs/filtered_contig_annotations.csv")
F07953 = read.csv("/tgen_labs/banovich/SingleCell/CellRanger/8_0_1/Projects/PEDCAR/VDJ/PEDCAR_0011_1_CS_Whole_C1_X5TCR_F07953_22YCMWLT3/outs/filtered_contig_annotations.csv")
F07954 = read.csv("/tgen_labs/banovich/SingleCell/CellRanger/8_0_1/Projects/PEDCAR/VDJ/PEDCAR_0012_1_CS_Whole_C1_X5TCR_F07954_22YCMWLT3/outs/filtered_contig_annotations.csv")


#### Process the contig data ----
# do lymphoid csf object for now
seurat_obj <- readRDS("/scratch/aoill/projects/CAR-T/00_new_2025/rds_files/csf_merge_rpca_subcluster_all_2025_09_05.rds")
metapath <- "/home/aoill/projects/CAR-T/00_new_2025/19130_all_batches_metadata_new_edited.csv"
sheet_name <- "Sheet1"
sample_metadata <- load_metadata(metadata_path = metapath, data_type = "csv", sheet_name = sheet_name)

seurat_obj_meta <- seurat_obj@meta.data
# remove any cols with TCR info because I hadn't added that info to the metadata
# originally 
seurat_obj_meta$TCRseq_CellRanger_path <- NULL
seurat_obj_meta$TCRseq_ID_long <- NULL
seurat_obj_meta$TCRseq_ID <- NULL
seurat_obj_meta_join <- left_join(seurat_obj_meta, sample_metadata)

# add back to metadata in object
seurat_obj@meta.data$TCRseq_CellRanger_path <- seurat_obj_meta_join$TCRseq_CellRanger_path
seurat_obj@meta.data$TCRseq_ID_long <- seurat_obj_meta_join$TCRseq_ID_long
seurat_obj@meta.data$TCRseq_ID <- seurat_obj_meta_join$TCRseq_ID


DimPlot(seurat_obj, reduction = "umap.integrated.rpca", 
        group.by = "sub.cluster",
        #split.by = "Batch",
        raster=FALSE,
        label = T
        #, ncol = 5
        ) + NoLegend()

DimPlot(seurat_obj, reduction = "umap.integrated.rpca", 
        group.by = "sub.cluster",
        split.by = "UPN_Cycle_Day_Sample_Type_Batch",
        raster=FALSE,
        label = T, ncol = 7
) + NoLegend()

seurat_list <- SplitObject(seurat_obj, split.by = "Cell_Prefix")
View(seurat_list
     )
length(seurat_list)
names(seurat_list)

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
names(seurat_list_cells_renamed) <- names(seurat_list) # these are the scRNA names not tcr names

# Merge all contigs and contig lists into one contig list
# not multiplexed samples
contig_list_csf <- list(F05038, F05039, F05040, F05042, 
                        F05043, F05044, F05045, 
                         
                        F05047, F05046, F05048#, 
                        #F07959, F07950, F07951, F07952, 
                        #F07953, F07954
                        )


names(contig_list_csf)
names(contig_list_csf) <- c("F05038_514_11_1_CSF", "F05039_515_4_2_CSF", "F05040_574_8_1_CSF", "F05042_514_8_1_CSF", 
                            "F05043_515_8_1_CSF", "F05044_574_3_1_CSF", "F05045_514_5_1_CSF", "F05047_515_4_2_CSF", 
                            "F05046_514_11_1_CSF", "F05048_574_8_1_CSF"
                            )

length(contig_list_csf)

# Batched data
# Get the contig list for the multiplexed batches
contig_list_11 <- createHTOContigList(F07959, seurat_list_cells_renamed[["Batch11"]], group.by = c("TCRseq_ID", "UPN", "Cycle", "Day", "Sample_Type"))
contig_list_2 <- createHTOContigList(F07950, seurat_list_cells_renamed[["Batch2"]], group.by = c("TCRseq_ID", "UPN", "Cycle", "Day", "Sample_Type"))
contig_list_3 <- createHTOContigList(F07951, seurat_list_cells_renamed[["Batch3"]], group.by = c("TCRseq_ID", "UPN", "Cycle", "Day", "Sample_Type"))
contig_list_4 <- createHTOContigList(F07952, seurat_list_cells_renamed[["Batch4"]], group.by = c("TCRseq_ID", "UPN", "Cycle", "Day", "Sample_Type"))
contig_list_5 <- createHTOContigList(F07953, seurat_list_cells_renamed[["Batch5"]], group.by = c("TCRseq_ID", "UPN", "Cycle", "Day", "Sample_Type"))
contig_list_6 <- createHTOContigList(F07954, seurat_list_cells_renamed[["Batch6"]], group.by = c("TCRseq_ID", "UPN", "Cycle", "Day", "Sample_Type"))

sum(length(names(contig_list_11)),
    length(names(contig_list_2)),
    length(names(contig_list_3)),
    length(names(contig_list_4)),
    length(names(contig_list_5)),
    length(names(contig_list_6))
)
# 16 
# 10 from the non multiplexed 

# Merge all contigs and contig lists into one contig list
contig_list_csf <- append(contig_list_csf, contig_list_11)
contig_list_csf <- append(contig_list_csf, contig_list_2)
contig_list_csf <- append(contig_list_csf, contig_list_3)
contig_list_csf <- append(contig_list_csf, contig_list_4)
contig_list_csf <- append(contig_list_csf, contig_list_5)
contig_list_csf <- append(contig_list_csf, contig_list_6)


names(contig_list_csf)
# NEW
#[1] "F05038_514_11_1_CSF" "F05039_515_4_2_CSF"  "F05040_574_8_1_CSF"  "F05042_514_8_1_CSF"  "F05043_515_8_1_CSF"  "F05044_574_3_1_CSF" 
#[7] "F05045_514_5_1_CSF"  "F05047_515_4_2_CSF"  "F05046_514_11_1_CSF" "F05048_574_8_1_CSF"  "F07959.514. 2.1.CSF" "F07959.689.17.1.CSF"
#[13] "F07959.705. 8.1.CSF" "F07950.574. 1.1.CSF" "F07950.716. 4.1.CSF" "F07950.689.12.1.CSF" "F07951.692.2.1.CSF"  "F07951.689.8.1.CSF" 
#[19] "F07952.514.13.1.CSF" "F07952.716. 8.1.CSF" "F07952.689. 4.1.CSF" "F07953.689.1.1.CSF"  "F07953.692.4.1.CSF"  "F07953.574.3.1.CSF" 
#[25] "F07954.705.5.1.CSF"  "F07954.716.2.1.CSF" 


# NEW
names(contig_list_csf) <- c("F05038_514_11_1_CSF",
  "F05039_515_4_2_CSF",
  "F05040_574_8_1_CSF",
  "F05042_514_8_1_CSF",
  "F05043_515_8_1_CSF",
  "F05044_574_3_1_CSF",
  "F05045_514_5_1_CSF",
  "F05047_515_4_2_CSF",
  "F05046_514_11_1_CSF",
  "F05048_574_8_1_CSF",
  "F07959_514_2_1_CSF",
  "F07959_689_17_1_CSF", 
  "F07959_705_8_1_CSF",
  "F07950_574_1_1_CSF",
  "F07950_716_4_1_CSF",
  "F07950_689_12_1_CSF",
  "F07951_692_2_1_CSF",
  "F07951_689_8_1_CSF", 
  "F07952_514_13_1_CSF",
  "F07952_716_8_1_CSF",
  "F07952_689_4_1_CSF",
  "F07953_689_1_1_CSF",
  "F07953_692_4_1_CSF",
  "F07953_574_3_1_CSF", 
  "F07954_705_5_1_CSF",
  "F07954_716_2_1_CSF" 
)

#### Combine the contigs ----
# Will want to use the parameter removeNA set to true to remove any cell barcode 
# with an NA value in at least one of the chains
# https://github.com/BorchLab/scRepertoire/issues/373
# NEW
combined <- combineTCR(contig_list_csf, 
                       samples =  c("F05038_514_11_1_CSF",
                                    "F05039_515_4_2_CSF",
                                    "F05040_574_8_1_CSF",
                                    "F05042_514_8_1_CSF",
                                    "F05043_515_8_1_CSF",
                                    "F05044_574_3_1_CSF",
                                    "F05045_514_5_1_CSF",
                                    "F05047_515_4_2_CSF",
                                    "F05046_514_11_1_CSF",
                                    "F05048_574_8_1_CSF",
                                    "F07959_514_2_1_CSF",
                                    "F07959_689_17_1_CSF", 
                                    "F07959_705_8_1_CSF",
                                    "F07950_574_1_1_CSF",
                                    "F07950_716_4_1_CSF",
                                    "F07950_689_12_1_CSF",
                                    "F07951_692_2_1_CSF",
                                    "F07951_689_8_1_CSF", 
                                    "F07952_514_13_1_CSF",
                                    "F07952_716_8_1_CSF",
                                    "F07952_689_4_1_CSF",
                                    "F07953_689_1_1_CSF",
                                    "F07953_692_4_1_CSF",
                                    "F07953_574_3_1_CSF", 
                                    "F07954_705_5_1_CSF",
                                    "F07954_716_2_1_CSF" 
                                    ),
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


# New
names(combined2) <- c("F05038_514_11_1_CSF",
                      "F05039_515_4_2_CSF",
                      "F05040_574_8_1_CSF",
                      "F05042_514_8_1_CSF",
                      "F05043_515_8_1_CSF",
                      "F05044_574_3_1_CSF",
                      "F05045_514_5_1_CSF",
                      "F05047_515_4_2_CSF",
                      "F05046_514_11_1_CSF",
                      "F05048_574_8_1_CSF",
                      "F07959_514_2_1_CSF",
                      "F07959_689_17_1_CSF", 
                      "F07959_705_8_1_CSF",
                      "F07950_574_1_1_CSF",
                      "F07950_716_4_1_CSF",
                      "F07950_689_12_1_CSF",
                      "F07951_692_2_1_CSF",
                      "F07951_689_8_1_CSF", 
                      "F07952_514_13_1_CSF",
                      "F07952_716_8_1_CSF",
                      "F07952_689_4_1_CSF",
                      "F07953_689_1_1_CSF",
                      "F07953_692_4_1_CSF",
                      "F07953_574_3_1_CSF", 
                      "F07954_705_5_1_CSF",
                      "F07954_716_2_1_CSF"
)
saveRDS(combined2, "/scratch/aoill/projects/CAR-T/00_new_2025/geo_2026/csf_scRepertoire_output.rds")


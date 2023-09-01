#==============================================================================#
# Author : Angela M. Oill, aoill@tgen.org
# Date: 2023/06/13
# Project: Pediatric CAR-T 
# Description: Figure S05
#==============================================================================#

#==============================================================================#
# Load libraries ----
#==============================================================================#
library(ggplot2)
library(Seurat)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(dittoSeq)
library(scProportionTest)
library(viridis)
library(ComplexHeatmap)


#==============================================================================#
# Read in integrated object ----
#==============================================================================#
pbmc_obj_T <- readRDS("/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/PBMC_T_cell_obj_2023.rds")


#==============================================================================#
# Set up cluster colors ----
#==============================================================================#
t_cell_colors <- c("CD4 Naive" = "#FEEA9A", "CD4 Effector" = "#FECD6A",
                   "CD4 Memory" = "#FDA245", "CD4 Memory (LTB+)" = "#FC6931", 
                   "CD4 Activated" = "#EA2820",
                   "CD4 Proliferating" = "#C20324", "Treg" = "#800026", 
                   
                   "CD8 Naive" = "#D5E3EF", "CD8 Effector" = "#A9C3DE",
                   "CD8 Effector Memory" = "#8C96C6", "CD8 Memory" = "#8A5DAA", 
                   "CD8 Resident Memory-like" = "#831F87",  "CD8 Central Memory (Activated)" = "#4D004B",
                   "CD8 Activated" = "#D3EDCC", "CD8 Proliferating" = "#98D493", 
                   "CD8 Exhausted" = "#4BAF61", 
                   "CD8 T (OX40+)" = "#147E3A", "NKT" = "#00441B"
                   
                   
)

# I don't want all of the labels to come up on the legend so use the code
# below to subset the ct_color character vector
# some of the colors are for cell states only in pbmcs or product
tmp_CTs <- as.character(levels(as.factor(pbmc_obj_T@meta.data$ct_T_states)))
tmp_color <- setNames(names(t_cell_colors), t_cell_colors) # swap names and values so I can subset
tmp_color_sub <- tmp_color[tmp_color %in% c(tmp_CTs)] # subset ct_color to only keep colors in PBMC UMAP
ct_colors_sub <- setNames(names(tmp_color_sub), tmp_color_sub) # re-swap names and values so that values are the color names 


#==============================================================================#
# Dotplot Heatmap ----
#==============================================================================#
# Curated set of genes
all_markers <- c("CD8A", "CD4", "FOXP3", "IL2RA", "CTLA4", 
                 "NKG7", "HOPX", "NCAM1", "CST7", "LEF1", "CCR7", "SATB1", "KLF2", 
                 "SELL", "IL7R", "RORA", "CXCR1", "GZMB", "GZMH", "PRF1", "KLRG1", 
                 "GNLY", "S100A4", "EOMES", "TCF7", "PDCD1", "CD69", "IFNG", 
                 "MKI67", "GZMK", "KLRD1", "CCL5", "ITGAE", "ZNF683", "HAVCR2", 
                 "TOX", "ENTPD1", "LAG3", "CXCL13", "LAYN"
)


pbmc_obj_T_cp <- pbmc_obj_T
DefaultAssay(pbmc_obj_T_cp) <- "SoupX_RNA" 
pbmc_obj_T_cp <- NormalizeData(pbmc_obj_T_cp)
VariableFeatures(pbmc_obj_T_cp) <- rownames(pbmc_obj_T_cp)
pbmc_obj_T_cp <- ScaleData(pbmc_obj_T_cp)

p <- DotPlot(pbmc_obj_T_cp, features = all_markers, 
             group.by = "ct_T_states", scale = TRUE) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
#p


# Reproducing the dotplot using ComplexHeatmap
df <- p$data
head(df)

# Adding log-transformed values
#df$log.avg.exp <- log10(df$avg.exp)

# (Scaled) expression levels
exp_mat <- df %>% 
  dplyr::select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame()

row.names(exp_mat) <- exp_mat$features.plot
exp_mat <- exp_mat[,-1] %>% as.matrix()

# The percentage of cells expressing a feature
percent_mat <- df %>% 
  dplyr::select(-avg.exp, -avg.exp.scaled) %>%  
  pivot_wider(names_from = id, values_from = pct.exp) %>% 
  as.data.frame() 

row.names(percent_mat) <- percent_mat$features.plot  
percent_mat <- percent_mat[,-1] %>% as.matrix()

dim(exp_mat)
dim(percent_mat)

# Get an idea of the ranges of the matrix
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99), na.rm=T)

# relace NA with 0 in matrix
exp_mat[!is.finite(exp_mat)] <- 0
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99))

# Any value that is greater than 2 will be mapped to yellow
col_fun = circlize::colorRamp2(c(-1, 0, 2), viridis(20)[c(1,10, 20)])

# Adding annotation colors
#library(scales) 
#hue_pal()(12) # hex <- hue_pal()(12) # there are 12 clusters
ct_colors_sub_heatmap <- c(list(PBMC = ct_colors_sub))

#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
cluster_anno <- colnames(exp_mat)
column_ha <- HeatmapAnnotation( # HeatmapAnnotation / rowAnnotation
  PBMC = cluster_anno,
  col = ct_colors_sub_heatmap, 
  #col = col, 
  na_col = "grey"
)

layer_fun = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.points(x = x, y = y, 
              gp = gpar(col = col_fun(pindex((exp_mat), i, j))), # t(exp_mat)
              size = pindex((percent_mat), i, j)/100 * unit(4, "mm"), # t(percent_mat)
              pch = 16
  )
}

lgd_list = list(
  Legend( labels = c(0,0.25,0.5,0.75,1), title = "pt", type = "points", 
          pch = 16, size = c(0,0.25,0.5,0.75,1) * unit(4,"mm"),
          legend_gp = gpar(col = "black")))

hp <- ComplexHeatmap::Heatmap((exp_mat), # t(exp_mat)
                              heatmap_legend_param=list(title="Scaled expression"),
                              column_title = "PBMC", 
                              col = col_fun,
                              rect_gp = gpar(type = "none"),
                              layer_fun = layer_fun,
                              row_names_gp = gpar(fontsize = 8),
                              column_names_gp = gpar(fontsize = 8),
                              #row_km = 4,
                              column_km = 4,
                              top_annotation = column_ha, # top_annotation / left_annotation
                              border = "black")


pdf("/home/aoill/plots/Figure_S05.pdf",
    width = 7, height = 7)
draw(hp, annotation_legend_list = lgd_list) 
dev.off()


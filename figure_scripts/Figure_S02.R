#==============================================================================#
# Author : Angela M. Oill, aoill@tgen.org
# Date: 2023/06/09
# Project: Pediatric CAR-T 
# Description: Figure S02
#==============================================================================#

#==============================================================================#
# Load libraries ----
#==============================================================================#
library(ggplot2)
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
# Read in integrated object
pbmc_obj <- readRDS("/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/PBMC_all_cells_obj_2023.rds")

#==============================================================================#
# Set up cluster colors ----
#==============================================================================#
ct_colors <- c("CD4+ T" = "#e09cb7", "CD8+ T" = "#d747ad", "NK" = "#a35776", 
               "Treg" = "#d54170", 
               "Monocytes" = "#bcd095", "CD14+ Monocytes" = "#6dce52", "FCGR3A+ Monocytes" = "#697f3f",
               "Macrophages" = "#c8d048",
               "DC" = "#de9e35", "cDC" = "#de9e35", "cDC1" = "#927754", "cDC2" = "#e3bb87", 
               "mDC" = "#a0712b",
               "pDC" = "#55daa2",  
               "B cells" = "#4a8773", "B" = "#4a8773", 
               "Plasma" = "#89d7c8", 
               #"Proliferating" = "#c7624f"
               "Proliferating" = "#9554b4"
               
)

# I don't want all of the labels to come up on the legend so use the code
# below to subset the ct_color character vector
tmp_CTs <- as.character(levels(as.factor(pbmc_obj@meta.data$CT)))
tmp_color <- setNames(names(ct_colors), ct_colors) # swap names and values so I can subset
tmp_color_sub <- tmp_color[tmp_color %in% c(tmp_CTs)] # subset ct_color to only keep colors in PBMC UMAP
ct_colors_sub <- setNames(names(tmp_color_sub), tmp_color_sub) # re-swap names and values so that values are the color names 


#==============================================================================#
# Dotplot Heatmap ----
#==============================================================================#
# Curated set of genes
all_markers <- c("CD3D", "CD3E", "CD8A", "CD4", "FOXP3", "IL2RA", "CTLA4", 
                 "GNLY", "NKG7", "KLRD1", "LYZ", "MARCO", "FCGR1A", "C1QA", 
                 "CD14", "FCGR3A", "MS4A7", "FCER1A", "CST3", "CLEC9A", "ITGAE", 
                 "ITGAX", "CD141", "XCR1", "CD1C", "CD207", "ITGAM", "NOTCH2", 
                 "SIRPA", "LAMP3", "CD83", "CCR7", "LILRA4", "CLEC4C", "JCHAIN", 
                 "IGHG1", "SDC1", "MS4A1", "CD79A", "CD19",
                 "MKI67"
)


pbmc_obj_cp <- pbmc_obj
DefaultAssay(pbmc_obj_cp) <- "SoupX_RNA" 
pbmc_obj_cp <- NormalizeData(pbmc_obj_cp)
VariableFeatures(pbmc_obj_cp) <- rownames(pbmc_obj_cp)
pbmc_obj_cp <- ScaleData(pbmc_obj_cp)

p <- DotPlot(pbmc_obj_cp, features = all_markers, 
             group.by = "CT", scale = TRUE) + 
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

hp <- Heatmap((exp_mat), # t(exp_mat)
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


pdf("/home/aoill/plots/Figure_S02.pdf",
    width = 7, height = 7)
draw(hp, annotation_legend_list = lgd_list) 
dev.off()



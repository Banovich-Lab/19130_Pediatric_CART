#==============================================================================#
# Author : Angela M. Oill, aoill@tgen.org
# Date: 2023/06/13
# Project: Pediatric CAR-T 
# Description: Figure S11
#==============================================================================#

#==============================================================================#
# Load libraries ----
#==============================================================================#
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(viridis)
library(ComplexHeatmap)


#==============================================================================#
# Read in integrated object ----
#==============================================================================#
tumor_obj <- readRDS("/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/Tumor_all_cells_obj_2023.rds")
tumor_obj_immune <- readRDS("/scratch/aoill/projects/CAR-T/00_final/00_GEO_rds/Tumor_immune_cells_obj_2023.rds")


#==============================================================================#
# All cells UMAP ----
#==============================================================================#
a <- DimPlot(tumor_obj, group.by = "integrated_sct_snn_res.0.2", label = F) +
  ggtitle("Tumor - All Cells") + NoLegend()

pa <- LabelClusters(a, id = "integrated_sct_snn_res.0.2", box = TRUE, size = 3.5, label.size = 0.1, 
                    box.padding = 0.2, seed = 309,
                    alpha = 0.8
)


#==============================================================================#
# All cells immune feature plot ----
#==============================================================================#
DefaultAssay(tumor_obj) <- "SCT"
pb <- FeaturePlot(tumor_obj, 
                  features = c("PTPRC"),
                  reduction = "integrated_sct_umap"
) 


#==============================================================================#
# Immune cells UMAP ----
#==============================================================================#
# Do colors and plotting like in figure 1
ct_colors <- c("T cells" = "#de7aac",
               "Myeloid" = "#5fa54b", #61a148
               "Monocytes" = "#bcd095", 
               "Macrophages" = "#c8d048",
               "cDCs" = "#de9e35", 
               "pDC" = "#55daa2"
               
)


c <- DimPlot(tumor_obj_immune, group.by = "CT", label = F)  +
  scale_color_manual(values = ct_colors) +
  ggtitle("Tumor - Immune Cells") + NoLegend()

pc <- LabelClusters(c, id = "CT", box = TRUE, size = 3.5, label.size = 0.01, 
                    box.padding = 0.7, seed = 309,
                    color = c("black", "white", rep("black", 4)),
                    alpha = 0.8
)


#==============================================================================#
# Immune cells dotplot heatmap ----
#==============================================================================#
# Curated set of genes
all_markers <- c("CD3D", "CD3E", "CD8A", "CD4", "FOXP3", "IL2RA", "CTLA4", 
                 "GNLY", "NKG7", "KLRD1", "LYZ", "MARCO", "FCGR1A", "C1QA", 
                 "CD14", "FCGR3A", "MS4A7", "FCER1A", "CST3", "CLEC9A", "ITGAE", 
                 "ITGAX", "CD141", "XCR1", "CD1C", "CD207", "ITGAM", "NOTCH2", 
                 "SIRPA", "LAMP3", "CD83", "CCR7", "LILRA4", "CLEC4C", "JCHAIN", 
                 "IGHG1", "SDC1", "MS4A1", "CD79A", "CD19"
)

ct_colors_sub <- ct_colors

tumor_obj_immune_cp <- tumor_obj_immune
DefaultAssay(tumor_obj_immune_cp) <- "SoupX_RNA" 
tumor_obj_immune_cp <- NormalizeData(tumor_obj_immune_cp)
VariableFeatures(tumor_obj_immune_cp) <- rownames(tumor_obj_immune_cp)
tumor_obj_immune_cp <- ScaleData(tumor_obj_immune_cp)

p <- DotPlot(tumor_obj_immune_cp, features = all_markers, 
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
ct_colors_sub_heatmap <- c(list(Tumor = ct_colors_sub))

#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
cluster_anno <- colnames(exp_mat)
column_ha <- HeatmapAnnotation( # HeatmapAnnotation / rowAnnotation
  Tumor = cluster_anno,
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
              column_title = "Tumor", 
              col = col_fun,
              rect_gp = gpar(type = "none"),
              layer_fun = layer_fun,
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 8),
              #row_km = 4,
              column_km = 4,
              top_annotation = column_ha, # top_annotation / left_annotation
              border = "black")

#ComplexHeatmap::Heatmap((exp_mat))

pd <- grid.grabExpr(draw(hp, annotation_legend_list = lgd_list,
                         merge_legend = TRUE, annotation_legend_side = "bottom", heatmap_legend_side = "bottom"))

pdf("/home/aoill/plots/Figure_S11.pdf",
    width = 9, height = 10)
ggarrange(ggarrange(pa, pb, pc, 
                    align = "hv", 
                    ncol = 1, 
                    nrow = 3,
                    labels = c("A", "B", "C")), 
          pd,
          ncol = 2,
          nrow = 1,
          widths = c(0.9, 1),
          labels = c("", "D"))
dev.off()

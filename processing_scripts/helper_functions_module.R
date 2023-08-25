# These functions were originally from Heini

library(Seurat)
library(glmGamPoi)

get_pcs <- function(seurat_object, reduction_name="pca") {
  
  # Determine percent of variation associated with each PC
  pct <- seurat_object[[reduction_name]]@stdev / sum(seurat_object[[reduction_name]]@stdev) * 100
  
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  
  # Determine which PC exhibits cumulative percent greater than 90% and % 
  # variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  
  co1
  
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  
  # last point where change of % of variation is more than 0.1%.
  co2
  
  # Minimum of the two calculation
  #pcs <- min(co1, co2)
  
  c(co1, co2)
  
}


basic_sct_rpca_integration <- function(batch_list, npcs=50, k_weight=100, numfeatures = 2500){
  features <- SelectIntegrationFeatures(object.list = batch_list,
                                        nfeatures = numfeatures)
  
  # Running PCA for rPCA integration
  for (i in 1:length(batch_list)) {
    # i <- 1
    message("Running PCA on list item ", i)
    batch_list[[i]] <- RunPCA(batch_list[[i]],
                              npcs = npcs, # changing npcs to avoid an error
                              features = features,
                              approx = F,
                              verbose = F)
  }
  
  message("Running PrepSCTIntegration()")
  batch_list <- PrepSCTIntegration(object.list = batch_list,
                                   anchor.features = features,
                                   verbose = F)
  
  # Identify anchors and integrate the datasets based on RNA
  # Using k=20 neighbors to find anchors (default = 5)?
  # Using rPCA to avoid "problem too large" error with CCA
  message("Running FindIntegrationAnchors()")
  anchors <- FindIntegrationAnchors(object.list = batch_list,
                                    normalization.method = "SCT", 
                                    anchor.features = features,
                                    #reference = seq(1, length(batch_list), by=3),
                                    k.anchor = 20,
                                    dims = 1:30,
                                    reduction = "rpca")
  # Default value for k.weight causes an error here (not enough cells for some 
  # samples?)
  message("Integrating data")
  integrated_object <- IntegrateData(anchorset = anchors,
                                     normalization.method = "SCT", 
                                     new.assay.name = "integrated_sct",
                                     k.weight = k_weight,
                                     verbose = T)
  # No need to run ScaleData if you've used SCT integration
  message("Running dimensionality reduction and clustering on integrated data")
  integrated_object <- RunPCA(integrated_object,
                              reduction.name = "integrated_sct_pca",
                              verbose = F)
  npcs <- min(get_pcs(integrated_object, reduction_name = "integrated_sct_pca"))
  message("Number of PCs to use: ", npcs)
  integrated_object <- RunUMAP(integrated_object,
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
  
  return(integrated_object)
}


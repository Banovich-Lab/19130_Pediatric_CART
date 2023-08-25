# README.md

Single cell analyses for 19130 pediatric CART project

## Scripts

### Preprocessing and integration
All sample types (CSF, PBMCs, product, and tumors) were preprocessed and integrated separately. We also integreated all CSF and PBMC samples together. The below scripts were used to process these data for downstream analyses and visualizations.

`19130_CART_CSF_processing.R`: This script reads in the 10X count matricies for all CSF samples, performs quality control (e.g. filtering and ambient RNA corrections) and other preprocessing steps (e.g. adding cell cycle scores, removing ribosomal and mitochondrial genes), and performs integration. Outputs 3 objects: 1) Seurat integrated object of all CSF cells, 2) Seurat integrated object of CSF T cells, 3) TCR clonotype call output file.


`19130_CART_PMBC_processing.R`: This script reads in the 10X count matricies for all PBMC samples, performs quality control (e.g. filtering and ambient RNA corrections) and other preprocessing steps (e.g. adding cell cycle scores, removing ribosomal and mitochondrial genes), and performs integration. Outputs 3 objects: 1) Seurat integrated object of all PBMC cells, 2) Seurat integrated object of PBMC T cells, 3) TCR clonotype call output file.


`19130_CART_Product_processing.R`: This script reads in the 10X count matricies for all product samples, performs quality control (e.g. filtering and ambient RNA corrections) and other preprocessing steps (e.g. adding cell cycle scores, removing ribosomal and mitochondrial genes), and performs integration. Outputs 2 objects: 1) Seurat integrated object of all product T cells, 2) TCR clonotype call output file.


`19130_CART_Tumor_processing.R`: This script reads in the 10X count matricies for all tumor samples, performs quality control (e.g. filtering and ambient RNA corrections) and other preprocessing steps (e.g. adding cell cycle scores, removing ribosomal and mitochondrial genes), and performs integration. Outputs 4 objects: 1) Seurat integrated object of all tumor cells, 2) Seurat integrated object of tumor immune cells, 3) Seurat integrated object of tumor T cells, 4) TCR clonotype call output file.


`19130_CART_CSF_PBMC_integration.R`: This script integrates all CSF and PBMC samples. Outputs 1 object: 1) Seurat integrated object of all cells.


Many of these scripts use helper functions which can be found here: https://github.com/Banovich-Lab/SingleCellBestPractices/blob/main/scripts/preprocessing_qc_module.R and here: https://github.com/Banovich-Lab/19130_Pediatric_CART/processing_scripts/helper_functions_module.R

### Manuscript figures
`/figure_scripts/` contain code to generate the visualizaitons referenced in the paper. All objects can be found on GEO (ADD INFO HERE).

`Figure_03.R`
- Inputs: CSF (`CSF_all_cells_obj_2023.rds`) and PBMC (`PBMC_all_cells_obj_2023.rds`) integrated objects.
- Output: `Figure_03.pdf`


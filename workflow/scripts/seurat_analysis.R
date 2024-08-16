library(Seurat)

# Load the data
sample <- snakemake@wildcards[['sample']]
input_file <- snakemake@input[['matrix_h5']]
output_file <- snakemake@output[['rds']]

# Load data into Seurat
seurat_data <- Read10X_h5(input_file)
seurat_object <- CreateSeuratObject(counts = seurat_data, project = sample)

# Perform standard pre-processing and clustering
seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object)
seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)
seurat_object <- RunUMAP(seurat_object, dims = 1:10)

# Save the Seurat object
saveRDS(seurat_object, file = output_file)

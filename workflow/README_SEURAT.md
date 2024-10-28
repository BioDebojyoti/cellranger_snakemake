## README for scQCAD

This pipeline takes aggregated cellranger output, and performs:

- Quality control,
- Batch correction (optional),
- Clustering
- Cell-type annotation using [SingleR](https://bioconductor.org/packages/release/bioc/vignettes/SingleR/inst/doc/SingleR.html)
- Differential expression (DE) analysis (optional)  
  using the Seurat package in R. It takes as input a sample sheet specifying metadata and directory paths for each sample and generates a range of outputs for QC, visualization, and downstream analysis.

## Input Requirements

The wrapper will need the [Cell Ranger v8.0.1](https://github.com/10XGenomics/cellranger/releases/tag/cellranger-8.0.1) output folder path or the h5 file following **aggr** step. It will also need a sample sheet (aggregate_csv_file) containing the meta-data for the different samples.
The aggregate_csv_file should be a CSV file with the following columns (it should follow the same order as used in cellranger aggr step). Column names 1,3, and 4 should be same as in the example below:

| sample_id | sample_outs        | donor | origin       | group     |
| :-------- | :----------------- | :---- | :----------- | :-------- |
| "S2"      | "/path/to/outs/S2" | "S2"  | "S2_Disease" | "Disease" |
| "S9"      | "/path/to/outs/S9" | "S9"  | "S9_Normal"  | "Normal"  |
| "S1"      | "/path/to/outs/S1" | "S1"  | "S1_Disease" | "Disease" |
| "S4"      | "/path/to/outs/S4" | "S4"  | "S4_Disease" | "Disease" |

## Output Files

scQCAD will produce figures in pdf format, and tables as CSV files. It will also produce two R objects in output. The R objects can be used later for additional or (re-analysis) analysisis. Here is a brief description of the output files content:

| filename                                                                                                                                                                                                                   | description                                                                                                                                                                        |
| :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Quality control (QC)**                                                                                                                                                                                                   |                                                                                                                                                                                    |
| <span style="color:red">pre_filter_qc_plot</span>                                                                                                                                                                          | pre-filtering QC metrics                                                                                                                                                           |
| <span style="color:green">post_filter_qc_plot</span>                                                                                                                                                                       | post-filtering QC metrics                                                                                                                                                          |
| **Principal Component Analysis (PCA)**                                                                                                                                                                                     |                                                                                                                                                                                    |
| <span style="color:green">**\<project\>**\_pca_plots.pdf</span>                                                                                                                                                            | PCA plots                                                                                                                                                                          |
| **Uniform Manifold Approximation and Projection (UMAP)**                                                                                                                                                                   |                                                                                                                                                                                    |
| <span style="color:orange">**\<project\>**\_no_integration_umap_plots.pdf</span>                                                                                                                                           | UMAP (conditonal on presence of batch variable) before integration                                                                                                                 |
| **Integration (UMAP)**                                                                                                                                                                                                     |                                                                                                                                                                                    |
| <span style="color:green">**\<project\>**\_**\<integration_method\>**\_integrated_umap_plots.pdf</span>                                                                                                                    | UMAP (conditonal on presence of batch variable) after integration                                                                                                                  |
| **Clustering (UMAP)**                                                                                                                                                                                                      |                                                                                                                                                                                    |
| <span style="color:green">**\<project\>**\_seurat_clusters_final_umap_plots.pdf</span> <br>OR</br> <span style="color:green">**\<project\>**\_**\<layer_column\>**\_integrated_seurat_clusters_final_umap_plots.pdf</span> | umap plot(s) of seurat clusters (conditonal on presence of batch variable **\<layer_column\>** after integration)                                                                  |
| **Annotation (UMAP)**                                                                                                                                                                                                      |                                                                                                                                                                                    |
| <span style="color:green">**\<project\>**\_annotated_umap_plots.pdf</span>                                                                                                                                                 | umap plot(s) of **singleR\*** annotated clusters (split by **\<condition_column\>** if present)                                                                                    |
| **Marker Identification**                                                                                                                                                                                                  |                                                                                                                                                                                    |
| <span style="color:green">**\<project\>\_\<cluster_name\>\_all_markers.pdf**</span> <br>OR</br> <span style="color:green">**\<project\>\_\<cluster_name\>\_\<condition_column\>\_all_markers.pdf**</span>                  | FeaturePlot of top1 markers, Heatmap of top10 markers, and DotPlot of top6 markers of **\<cluster_name\>**. FeaturePlot and DotPlots split by **\<condition_column\>** if present. |

**\<project\>** default name is "singleCell"\
**\<integration_method\>** no default integration-method\
**\<layer_column\>** no default batch-variable\
**\<condition_column\>** no default group-variable\
**singleR\*** currently only supports species "human", and "mouse".\
**\<cluster_name\>** "seurat_clusters" and "singleR.labels".

Currently, the following references are used for cell type annotation

```r
# species = "human"
celldex::HumanPrimaryCellAtlasData()
# species = "mouse"
celldex::MouseRNAseqData()
```

The wrapper run saves two R objects

```r
"path/to/<seurat_out_dir>/<project>_raw_seurat.rds"
"path/to/<seurat_out_dir>/<project>_analysed_seurat.rds"
```

The former can be used to re-analyse without having to create the seurat object afresh. The latter object is meant to help perform additional analysis should there be a need for it. Please follow the steps described below to load the object of interest.

```r
library(Seurat)
# for re-analysis
seurat_obj <- readRDS("path/to/<seurat_out_dir>/<project>_raw_seurat.rds")
# for additional analysis
seurat_obj <- readRDS("path/to/<seurat_out_dir>/<project>_analysed_seurat.rds")
```

While the wrapper is expected to be used by command-line invocation it can also be used as a function by first making it available as shown below:

## Make the functions available

```r
source("scQCAD.R")
```

## Run the analysis function

```r
# Run Seurat analysis
seurat_analysis(
data_dir <- "count/filtered_feature_bc_matrix"
data_file <- NULL
project_name <- "test"
seurat_out_dir <- "seurat_out"
min_cells <- 3
min_features <- 100
max_features <- 3000
percent_mt <- NULL
percent_rb <- NULL
aggr_csv_file <- "aggregation.csv"
tcr_file <- "vdj_t/filtered_contig_annotations.csv"
bcr_file <- "vdj_b/filtered_contig_annotations.csv"
layer_column <- "donor"
condition_column <- "health_status"
integration_method <- "RPCAIntegration"
# CCAIntegration, RPCAIntegration, HarmonyIntegration,
# FastMNNIntegration, scVIIntegration
enable_sct <- TRUE
perform_de <- FALSE
species <- "human"
)
```

## CLI options

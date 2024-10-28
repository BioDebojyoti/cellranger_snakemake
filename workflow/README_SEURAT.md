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

In this example:

sample_sheet specifies the path to the sample sheet CSV.
output_directory is where all output files will be saved.
qc_thresholds is a list of QC filtering thresholds.
perform_DE is a logical flag for performing differential expression analysis.
group_by specifies the metadata column used to group cells.
condition_column sets the condition used in DE analysis (e.g., origin). 3. Script Overview
3.1. load_and_aggregate_samples
This function:

Reads each sample listed in the sample sheet.
Loads the sample data into a Seurat object.
Aggregates all samples into a single Seurat object for batch integration and QC.
3.2. seurat_qc_plot
This function:

Creates QC plots (violin and scatter) to inspect metrics like nCount_RNA, nFeature_RNA, percent.mt, and percent.rb.
Saves these plots as qc_plots.pdf in the output directory.
3.3. seurat_normalize_scale_and_pca
This function:

Normalizes and scales the aggregated Seurat object.
Computes principal components and generates PCA plots, saving them as pca_plots.pdf.
3.4. seurat_run_umap
This function:

Runs UMAP on the Seurat object and generates UMAP plots colored by metadata columns (e.g., origin and health_status).
Saves UMAP plots as umap_plots.pdf.
3.5. find_markers_and_DE_analysis
If perform_DE is set to TRUE, this function:

Identifies cluster-specific markers and performs DE analysis between groups specified in group_by and condition_column.
Outputs DE results as CSV files and visualizes them with EnhancedVolcano plots, saved as de_genes.pdf.
3.6. plot_markers
This function:

Generates expression plots for key marker genes from DE analysis.
Saves plots as all_markers.pdf. 4. Output Description
The main output files generated by this script include:

4.1. Quality Control Plots (qc_plots.pdf)
Contains various QC plots, including:

Violin plots for nFeature_RNA, nCount_RNA, percent.mt, and percent.rb.
Density distributions of nCount values.
Scatter plots showing correlations between key QC metrics.
4.2. PCA Plots (pca_plots.pdf)
Displays PCA results, including:

PCA loadings for dimensions 1 and 2.
Scatter plots and elbow plots to visualize the variance captured by each principal component.
4.3. UMAP Plots (umap_plots.pdf)
Presents UMAP visualizations of the data, colored by the specified metadata columns.

4.4. Differential Expression Analysis Results
If perform_DE is enabled, the output will include:

CSV Files (de*results*<group>.csv): Lists DE genes between conditions, including columns for log2 fold change, p-values, and adjusted p-values.
Violin Plots for DE Genes (de_violin_plots.pdf): Violin plots for top DE genes.
Volcano Plots (de_volcano_plots.pdf): Volcano plots showing significant DE genes between specified groups.
4.5. Marker Gene Plots (all_markers.pdf)
Contains:

Feature plots of key marker genes across clusters or conditions.
Heatmaps and dot plots highlighting marker expression patterns. 5. Troubleshooting & Tips
Memory Management: For large datasets, consider increasing R's memory limit.
Customizing DE Analysis: Adjust qc_thresholds or specify alternative condition_column values as needed.
Custom Outputs: Modify plot_markers and seurat_qc_plot to include additional metadata columns if required.
This pipeline allows for the efficient aggregation, QC, normalization, visualization, and DE analysis of single-cell RNA-seq data, yielding high-quality outputs suitable for further exploration or publication.

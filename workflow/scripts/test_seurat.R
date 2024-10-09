#!/opt/sw/bioinfo-tools/sources/anaconda3/envs/d_seurat510/bin/Rscript



# Function to load libraries silently
load_library <- function(x) {
  invisible(
    suppressMessages(
      library(
        x,
        character.only = TRUE,
        quietly = TRUE,
        warn.conflicts = FALSE
      )
    )
  )
}

# List of required libraries
libraries2use <- c(
  "hdf5r", "leidenAlg", "igraph", "patchwork", "ggplot2",
  "Seurat", "optparse", "dplyr", "remotes", "celldex", "SingleR",
  "writexl", "seuratter", "Seurat", "SeuratObject", "harmony",
  "seuratHelper", "this.path"
)

# Load all libraries
invisible(sapply(libraries2use, load_library))

# get directory for the present file
script_path <- this.path::this.path()
print(script_path)
source(paste0(dirname(script_path), "/useful.R", collapse = ""))
# Set options for command-line arguments
option_list <- list(
  make_option(
    c("-f", "--file"),
    type = "character", default = "",
    help = "count matrix data h5 file name",
    dest = "count_matrix_file"
  ),
  make_option(
    c("-a", "--aggregate-csv"),
    type = "character", default = "",
    help = "aggregate csv file", dest = "aggregate_csv_file"
  ),
  make_option(
    c("-d", "--directory"),
    type = "character", default = "",
    help = "count matrix data directory name",
    dest = "count_matrix_dir"
  ),
  make_option(
    c("-o", "--out-directory"),
    type = "character", default = "seurat_out",
    help = "output directory", dest = "outdir"
  ),
  make_option(
    c("--min-cells"),
    type = "integer", default = as.integer(3),
    help = "minimum cells [default= %default]", dest = "min.cells"
  ),
  make_option(
    c("--min-features"),
    type = "integer", default = as.integer(200),
    help = "minimum features [default= %default]", dest = "min.features"
  ),
  make_option(
    c("--max-features"),
    type = "integer", default = as.integer(2500),
    help = "maximum features [default= %default]", dest = "max.features"
  ),
  make_option(
    c("--percent-mt"),
    type = "double", default = as.double(5),
    help = "threshold percent mitochondrial [default= %default]",
    dest = "cutoff.percent.mt"
  ),
  make_option(
    c("--percent-rb"),
    type = "double", default = as.double(5),
    help = "threshold percent ribosomal [default= %default]",
    dest = "cutoff.percent.rb"
  ),
  make_option(
    c("--project"),
    type = "character", default = "singleCell",
    help = "output file name [default= %default]", dest = "project"
  ),
  make_option(
    c("--vdj-t"),
    type = "character", default = "",
    help = "V(D)J-T annotations", dest = "vdj_t"
  ),
  make_option(
    c("--vdj-b"),
    type = "character", default = "",
    help = "V(D)J-B annotations", dest = "vdj_b"
  ),
  make_option(
    c("--layer-column"),
    type = "character",
    default = ""
  ),
  make_option(
    c("--integration-method"),
    type = "character", default = "RPCAIntegration",
    help = "integration method
    (CCAIntegration, RPCAIntegration, HarmonyIntegration,
    FastMNNIntegration, scVIIntegration)",
    dest = "integration_method"
  )
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Define user inputs
project_name <- opt$project
seurat_out_dir <- opt$outdir
min_cells <- opt$min.cells
min_features <- opt$min.features
max_features <- opt$max.features
percent_mt <- opt$cutoff.percent.mt
percent_rb <- opt$cutoff.percent.rb
data_dir <- opt$count_matrix_dir
data_file <- opt$count_matrix_file
aggr_csv_file <- opt$aggregate_csv_file
tcr_file <- opt$vdj_t
bcr_file <- opt$vdj_b
layer_column <- opt$layer_column
integration_method <- opt$integration_method

# Main Seurat analysis function
seurat_analysis <- function(
    data_dir,
    data_file,
    project_name,
    seurat_out_dir,
    min_cells,
    min_features,
    max_features,
    percent_mt,
    percent_rb,
    aggr_csv_file,
    tcr_file,
    bcr_file,
    layer_column,
    integration_method) {
  if (!dir.exists(seurat_out_dir)) {
    dir.create(seurat_out_dir)
  }

  # Load dataset
  data_10_x <- load_data_10X(data_dir, data_file)

  # Load metadata if applicable
  sample_identity <- metadata(data_10_x, aggr_csv_file)

  # Create Seurat object with multiple layers
  seurat_obj <- create_seurat_object(
    data_10_x,
    sample_identity,
    min_cells,
    min_features,
    project_name
  )

  # Add V(D)J-T and B annotations if applicable
  seurat_obj <- add_clonotype(tcr_file, seurat_obj, "t")
  seurat_obj <- add_clonotype(bcr_file, seurat_obj, "b")


  # Quality control
  seurat_obj <- seurat_qc(seurat_obj, ncol, seurat_out_dir, project_name)

  # split if  datasets ran across experimental batches, donors, or condition
  seurat_obj <- split_seurat_by_layer(
    seurat_obj, 
    assay= "RNA", 
    layer_column = layer_column
    )
  seurat_obj <- seurat_normalization(seurat_obj, sct = enable_SCT)
  seurat_obj <- seurat_dimensional_reduction(seurat_obj)

  # Integrative analysis using the specified method
  if (integration_method != "") {
    seurat_obj <- integrate_seurat_layers(
      seurat_obj,
      layer_column,
      integration_method
    )

    seurat_obj[["RNA"]] <- SeuratObject::JoinLayers(
      seurat_obj[["RNA"]]
    )

    seurat_obj <- Seurat::FindNeighbors(
      seurat_obj,
      reduction = method_dict[[integration_method]][2],
      dims = 1:30
    )
    seurat_obj <- Seurat::FindClusters(
      seurat_obj,
      resolution = 0.8,
      random.seed = 234
    )

    seurat_obj <- Seurat::RunUMAP(
      seurat_obj,
      dims = 1:30,
      reduction = method_dict[[integration_method]][2],
      reduction.name = "umap.integrated"
    )
    umap <- Seurat::DimPlot(
      seurat_obj,
      reduction = "umap.integrated",
      group.by = c(layer_column, "seurat_clusters"),
      combine = TRUE
    )

    umap_by_condition <- Seurat::DimPlot(
      seurat_obj,
      reduction = "umap.integrated",
      split.by = layer_column
    )
  } else {
    seurat_obj <- Seurat::FindNeighbors(
      seurat_obj,
      dims = 1:30,
      reduction = "pca"
    )
    seurat_obj <- Seurat::FindClusters(
      seurat_obj,
      resolution = 0.8,
      random.seed = 234,
      cluster.name = "unintegrated_clusters"
    )

    seurat_obj <- Seurat::RunUMAP(
      seurat_obj,
      dims = 1:30,
      reduction = "pca",
      reduction.name = "umap"
    )

    group_by_list <- c("seurat_clusters")

    # if grouping variable exists add it to the list
    if (layer_column != "") {
      group_by_list <- c(layer_column, group_by_list)
    }

    umap <- Seurat::DimPlot(
      seurat_obj,
      reduction = "umap",
      group.by = group_by_list,
      combine = TRUE
    )

    # if grouping variable batch/condition/donor exists
    if (layer_column != "") {
      umap_by_condition <- Seurat::DimPlot(
        seurat_obj,
        reduction = "umap",
        split.by = layer_column,
        ncol = 2,
        combine = TRUE
      )
    }
  }


  # Save Seurat object
  if (integration_method != "") {
    file2save <- paste0(
      seurat_out_dir, "/", project_name, "_integrated_seurat.rds"
    )
  } else {
    file2save <- paste0(seurat_out_dir, "/", project_name, "_seurat.rds")
  }
  saveRDS(
    seurat_obj,
    filename = file2save
  )

  # Visualization and QC plots
  violin_plots <- Seurat::VlnPlot(
    seurat_obj,
    features = c(
      "nFeature_RNA", "nCount_RNA",
      "percent.mt", "percent.rb"
    ),
    ncol = 2
  )
  pdf(paste0(seurat_out_dir, "/", project_name, "_violin_features.pdf"))
  print(violin_plots)
  dev.off()

  umap_plot <- Seurat::DimPlot(
    seurat_obj,
    reduction = "umap",
    combine = TRUE
  )
  pdf(paste0(seurat_out_dir, "/", project_name, "_umap.pdf"))
  print(umap_plot)
  dev.off()

  # Save final integrated Seurat object
  saveRDS(
    seurat_obj,
    file = paste0(
      seurat_out_dir, "/", project_name,
      "_final_integrated_analysis.Rds"
    )
  )
}


data_dir <- "count/sample_filtered_feature_bc_matrix"
data_file <- ""
project_name <- "test"
seurat_out_dir <- "seurat_out"
min_cells <- 3
min_features <- 3
max_features <- 2500
percent_mt <- 5.0
percent_rb <- 5.0
aggr_csv_file <- ""
tcr_file <- "vdj_t/filtered_contig_annotations.csv"
bcr_file <- ""
layer_column <- ""
integration_method <- "" # "RPCAIntegration"
# CCAIntegration, RPCAIntegration, HarmonyIntegration,
# FastMNNIntegration, scVIIntegration

# Run Seurat analysis
seurat_analysis(
  data_dir,
  data_file,
  project_name,
  seurat_out_dir,
  min_cells,
  min_features,
  max_features,
  percent_mt,
  percent_rb,
  aggr_csv_file,
  tcr_file,
  bcr_file,
  layer_column,
  integration_method
)

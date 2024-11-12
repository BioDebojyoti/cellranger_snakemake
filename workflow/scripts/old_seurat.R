#!/opt/sw/bioinfo-tools/sources/anaconda3/envs/d_seurat510/bin/Rscript

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

libraries2use <- c(
  "hdf5r", "leidenAlg", "igraph", "patchwork", "ggplot2",
  "Seurat", "optparse", "dplyr", "remotes", "celldex", "SingleR",
  "writexl", "seuratter"
)

invisible(sapply(libraries2use, load_library))


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
    help = "output directory",
    dest = "outdir"
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
    help = "threshold percent mitrochondrial [default= %default]",
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
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

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




metadata <- function(data_10_x, aggregate_csv_file) {
  # Load metadata
  print(aggregate_csv_file)
  samples <- read.csv(aggregate_csv_file, stringsAsFactors = FALSE)
  cells <- new("seurat", raw.data = data_10_x)
  # Get barcodes and suffix:
  cellcodes <- as.data.frame(cells@raw.data@Dimnames[[2]])
  colnames(cellcodes) <- "barcodes"
  rownames(cellcodes) <- cellcodes$barcodes

  cellcodes$libcodes <- as.factor(
    gsub(pattern = ".+-", replacement = "", cellcodes$barcodes)
  )
  cellcodes$samples <- as.vector(samples$sample_id[cellcodes$libcodes])

  # Create dataframe for meta.data argument and set up object:
  sampleidentity <- cellcodes["samples"]
  return(sampleidentity)
}

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
    bcr_file) {
  if (!dir.exists(seurat_out_dir)) {
    dir.create(seurat_out_dir)
  }

  print(seurat_out_dir)

  # Load the dataset
  if (dir.exists(data_dir)) {
    data_10_x <- Seurat::Read10X(data.dir = data_dir)
  } else {
    data_10_x <- Seurat::Read10X_h5(filename = data_file)
  }

  if (file.exists(aggr_csv_file)) {
    sample_identity <- metadata(data_10_x, aggr_csv_file)
  } else {
    sample_identity <- NULL
  }
  # Initialize the Seurat object with the raw (non-normalized data).
  seurat_obj <- Seurat::CreateSeuratObject(
    data_10_x,
    meta.data = sample_identity,
    min.cells = min_cells,
    min.features = min_features,
    project = project_name
  )

  # Standard pre-processing workflow
  # QC metrics
  # Low-quality/ dying cells often exhibit extensive
  # mitochondrial contamination
  seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(
    seurat_obj,
    pattern = "^MT-"
  )
  seurat_obj[["percent.rb"]] <- Seurat::PercentageFeatureSet(
    seurat_obj,
    pattern = "^RP[SL]"
  )


  violin_plots <- Seurat::VlnPlot(
    seurat_obj,
    features = c(
      "nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"
    ),
    ncol = 2
  )
  pdf(paste0(seurat_out_dir, "/", project_name, "_violin_features.pdf"))
  print(violin_plots)
  dev.off()

  plot1 <- Seurat::FeatureScatter(
    seurat_obj,
    feature1 = "nCount_RNA", feature2 = "nFeature_RNA"
  )
  plot2 <- Seurat::FeatureScatter(
    seurat_obj,
    feature1 = "nCount_RNA", feature2 = "percent.mt"
  )
  plot3 <- Seurat::FeatureScatter(
    seurat_obj,
    feature1 = "nCount_RNA", feature2 = "percent.rb"
  )
  plot4 <- Seurat::FeatureScatter(
    seurat_obj,
    feature1 = "percent.mt", feature2 = "percent.rb"
  )

  plot_features <- cowplot::plot_grid(
    plot1, plot2, plot3, plot4,
    labels = c("A", "B", "C", "D"),
    ncol = 2, align = c("h", "v"),
    label_size = 20
  )

  pdf(paste0(seurat_out_dir, "/", project_name, "_features.pdf"))
  print(plot_features)
  dev.off()

  add_clonotype <- function(filepath, seurat_obj, type = "t") {
    fp <- strsplit(filepath, "/")[[1]]

    vdj_path <- paste0(fp[-length(fp)], collapse = "/")

    tcr <- read.csv(paste0(
      vdj_path,
      "/filtered_contig_annotations.csv",
      collapse = ""
    ))

    # Remove the -1 at the end of each barcode.
    # Subsets so only the first line of each barcode is kept,
    # as each entry for given barcode will have same clonotype.
    tcr <- tcr[!duplicated(tcr$barcode), ]

    # Only keep the barcode and clonotype columns.
    # We'll get additional clonotype info from the clonotype table.
    tcr <- tcr[, c("barcode", "raw_clonotype_id")]
    names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"

    # Clonotype-centric info.
    clono <- read.csv(
      paste0(
        vdj_path,
        "/clonotypes.csv",
        collapse = ""
      )
    )

    # Slap the AA sequences onto our original table by clonotype_id.
    tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")])
    names(tcr)[names(tcr) == "cdr3s_aa"] <- "cdr3s_aa"

    # Reorder so barcodes are first column and set them as rownames.
    tcr <- tcr[, c(2, 1, 3)]
    rownames(tcr) <- tcr[, 1]
    tcr[, 1] <- NULL
    colnames(tcr) <- paste(type, colnames(tcr), sep = "_")
    # Add to the Seurat object's metadata.
    clono_seurat <- Seurat::AddMetaData(object = seurat_obj, metadata = tcr)
    return(clono_seurat)
  }

  # T cell and B cell clonetype information
  if (file.exists(tcr_file)) {
    seurat_obj <- add_clonotype(
      tcr_file,
      seurat_obj,
      "t"
    )
  }

  if (file.exists(bcr_file)) {
    seurat_obj <- add_clonotype(
      bcr_file,
      seurat_obj,
      "b"
    )
  }


  # filter out low quality samples
  seurat_obj <- subset(
    seurat_obj,
    nFeature_RNA > min_features & nFeature_RNA < max_features & percent.mt < percent_mt
  )

  # Normalization
  seurat_obj <- Seurat::NormalizeData(
    seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )

  browser()
  # highly variable feature selection
  seurat_obj <- Seurat::FindVariableFeatures(
    seurat_obj,
    selection.method = "vst",
    nfeatures = 2000
  )

  # Identify the 10 most highly variable genes
  top10 <- head(Seurat::VariableFeatures(seurat_obj), 10)

  # plot variable features with and without labels
  plot1 <- Seurat::VariableFeaturePlot(seurat_obj)
  plot2 <- Seurat::LabelPoints(plot = plot1, points = top10, repel = TRUE)
  combined_plot <- plot1 + plot2
  pdf(paste0(seurat_out_dir, "/", project_name, "_variable_features.pdf"))
  print(combined_plot)
  dev.off()

  # scale data
  all.genes <- rownames(seurat_obj)
  seurat_obj <- Seurat::ScaleData(seurat_obj, features = all.genes)


  # linear dimensional reduction
  seurat_obj <- Seurat::RunPCA(
    seurat_obj,
    features = VariableFeatures(object = seurat_obj)
  )

  # Determine the dimensionality of the dataset
  dimensionality <- determine_dimensionality(
    ElbowPlot(seurat_obj)$data,
    cutoff = 0.95
  )

  print(seurat_obj)
  print(dimensionality)

  pca_plot <- DimPlot(
    seurat_obj,
    reduction = "pca",
    dims = c(1, 2),
    combine = TRUE
  )
  # PCA plot
  pdf(paste0(seurat_out_dir, "/", project_name, "_pca.pdf"))
  print(pca_plot)
  dev.off()

  # Cluster the cells
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:dimensionality)
  seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 0.1)

  # Run non-linear dimensional reduction (UMAP/tSNE)
  seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:dimensionality)

  # UMAP plot
  umap_plot <- Seurat::DimPlot(
    seurat_obj,
    reduction = "umap",
    combine = TRUE
  )
  pdf(paste0(seurat_out_dir, "/", project_name, "_umap.pdf"))
  print(umap_plot)
  dev.off()

  # celldex::
  # celldex::NovershternHematopoieticData
  # celldex::MouseRNAseqData
  # celldex::HumanPrimaryCellAtlasData
  # celldex::ImmGenData
  # celldex::DatabaseImmuneCellExpressionData celldex::MonacoImmuneData
  # celldex::BlueprintEncodeData
  # ref <- celldex::HumanPrimaryCellAtlasData()
  # input4singler <- Seurat::GetAssayData(seurat_obj)

  # pred <- SingleR(test = input4singler, ref = ref, labels = ref$label.main)
  # seurat_obj[["SingleR.labels"]] <- pred$labels

  # writexl::write_xlsx(
  #   x = de_dfs,
  #   path = paste0(
  #     seurat_out_dir, "/", project_name, "_annotated_cluster_DE_results.xlsx"
  #   )
  # )

  # DE between each pair of cluster/annotation groups
  # group_pair <- combn(unique(seurat_obj@meta.data$seurat_clusters), 2, list)
  # de_dfs <- list()
  # for (p in group_pair) {
  #   contrast <- paste0(p[[1]], "_Vs_", p[[2]])
  #   print(contrast)
  #   curr_de <- FindMarkers(
  #     seurat_obj,
  #     ident.1 = p[[1]],
  #     ident.2 = p[[2]],
  #     group.by = "SingleR.labels",
  #     logfc.threshold = log(2)
  #   )
  #   de_dfs[[contrast]] <- de_dfs
  # }

  # writexl::write_xlsx(
  #   x = de_dfs,
  #   path = paste0(
  #     seurat_out_dir, "/", project_name, "_annotated_cluster_DE_results.xlsx"
  #   )
  # )



  saveRDS(
    seurat_obj,
    file = paste0(seurat_out_dir, "/", project_name, "_seurat_analysis.Rds")
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
  bcr_file
)

load_data_10X <- function(data_dir="", data_file=""){
   if (dir.exists(data_dir)) {
    data_10_x <- Seurat::Read10X(data.dir = data_dir)
    return(data_10_x)
  } 
  if (dir.exists(data_file)) {
    data_10_x <- Seurat::Read10X_h5(filename = data_file)
    return(data_10_x)
  }
  print("10X data not found!")
}

# Function to load metadata and integrate cell data
metadata <- function(data_10_x, aggregate_csv_file) {
  if (file.exists(aggregate_csv_file)) {
    samples <- read.csv(aggregate_csv_file, stringsAsFactors = FALSE)
    cells <- new("seurat", raw.data = data_10_x)

    cellcodes <- as.data.frame(cells@raw.data@Dimnames[[2]])
    colnames(cellcodes) <- "barcodes"
    rownames(cellcodes) <- cellcodes$barcodes
    cellcodes$libcodes <- as.factor(
      gsub(pattern = ".+-", replacement = "", cellcodes$barcodes)
    )
    cellcodes$samples <- as.vector(samples$sample_id[cellcodes$libcodes])

    sampleidentity <- cellcodes["samples"]
    return(sampleidentity)
  } else {
    return(NULL)
  }
}

acronym <- function(name_list) {
  acronym_list <- lapply(
    lapply(
      sapply(name_list, function(x) stringr::str_split(x, " ")),
      function(y) stringr::str_sub(y, 1, 1)
    ),
    function(w) paste0(w, collapse = "")
  )
  return(acronym_list)
}
# Generic function to create a Seurat object from Gene Expression and
# add additional assays
create_seurat_object <- function(
    data_10_x, sample_identity, min_cells, min_features, project_name) {
  # Create base Seurat object using Gene Expression
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = data_10_x[["Gene Expression"]],
    meta.data = sample_identity,
    min.cells = min_cells,
    min.features = min_features,
    project = project_name
  )

  # Check if additional assays are present and add them
  additional_assays <- setdiff(names(data_10_x), "Gene Expression")
  acronym_list <- acronym(additional_assays)
  for (assay in additional_assays) {
    current_assay <- SeuratObject::CreateAssay5Object(
      counts = data_10_x[[assay]]
    )
    seurat_obj[[acronym_list[assay][[1]]]] <- current_assay
  }

  return(seurat_obj)
}

seurat_qc <- function(seurat_obj) {
  # Add mitochondrial and ribosomal percentages
  seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(
    seurat_obj,
    pattern = "^MT-"
  )
  seurat_obj[["percent.rb"]] <- Seurat::PercentageFeatureSet(
    seurat_obj,
    pattern = "^RP[SL]"
  )
  return(seurat_obj)
}

split_seurat_by_layer <- function(seurat_obj, assay = "RNA", layer_column = ""){
  if(layer_column == ""){
  seurat_obj[["RNA"]] <- S4Vectors::split(
      seurat_obj[["RNA"]],
      f = seurat_obj@meta.data[[layer_column]]
    )
    }
  return(seurat_obj)   
}

method_dict <- list(
  "CCAIntegration" = c("pca", "integrated.cca"),
  "RPCAIntegration" = c("pca", "integrated.rpca"),
  "HarmonyIntegration" = c("pca", "harmony"),
  "JointPCAIntegration" = c("pca", "integrated.dr")
)

# Function to integrate Seurat layers using specified method
integrate_seurat_layers <- function(seurat_obj, layer_column, method) {
  if(layer_colum ==""){
    return(seurat_obj)
  } else {
  seurat_obj <- Seurat::IntegrateLayers(
    seurat_obj,
    method = method,
    orig.reduction = method_dict[[method]][1],
    new.reduction = method_dict[[method]][2]
  )
  return(seurat_obj)
  }
}

seurat_normalization <- function(seurat_obj, sct = FALSE) {
  if(sct) {
    seurat_obj <- Seurat::SCTransform(seurat_obj, verbose = FALSE)
    return(seurat_obj)
    } else {
    seurat_obj <- Seurat::NormalizeData(seurat_obj)
    seurat_obj <- Seurat::FindVariableFeatures(seurat_obj)
    seurat_obj <- Seurat::ScaleData(seurat_obj)
    return(seurat_obj)
    }
}

seurat_dimensional_reduction <- function(seurat_obj) {
  seurat_obj <- Seurat::RunPCA(seurat_obj, verbose = FALSE)
  seurat_obj <- Seurat::RunUMAP(seurat_obj, verbose = FALSE)
}



add_clonotype <- function(tcr_file, seurat_obj, type = "t") {
  if (tcr_file == "") {
    return(seurat_obj)
  }

  tcr_prefix <- dirname(normalizePath(tcr_file))
  tcr <- read.csv(
    paste(tcr_prefix, "filtered_contig_annotations.csv", sep = "/")
  )

  # Remove the -1 at the end of each barcode.
  # Subsets so only the first line of each barcode is kept,
  # as each entry for given barcode will have same clonotype.
  tcr <- tcr[!duplicated(tcr$barcode), ]

  # Only keep the barcode and clonotype columns.
  # We'll get additional clonotype info from the clonotype table.
  tcr <- tcr[, c("barcode", "raw_clonotype_id")]
  names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"

  # Clonotype-centric info.
  clono <- read.csv(paste(tcr_prefix, "clonotypes.csv", sep = "/"))

  # Slap the AA sequences onto our original table by clonotype_id.
  tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")])
  names(tcr)[names(tcr) == "cdr3s_aa"] <- "cdr3s_aa"

  # Reorder so barcodes are first column and set them as rownames.
  tcr <- tcr[, c(2, 1, 3)]
  rownames(tcr) <- tcr[, 1]
  tcr[, 1] <- NULL
  colnames(tcr) <- paste(type, colnames(tcr), sep = "_")
  # Add to the Seurat object's metadata.
  seurat_obj <- SeuratObject::AddMetaData(object = seurat_obj, metadata = tcr)
  return(seurat_obj)
}


# Plot Functions 
Plot_QC <- function(seurat_obj, ncol, seurat_out_dir, project_name) {
  vp <- Seurat::VlnPlot(
    seurat_obj,
    features = c(
      "nFeature_RNA", "nCount_RNA",
      "percent.mt", "percent.rb"
    ),
    group.by = "orig.ident",
    ncol = ncol,
    pt.size = 0.1
  ) & theme(plot.title = element_text(size = 16))
  fs_1 <- FeatureScatter(
    seurat_obj,
    feature1 = "nCount_RNA",
    feature2 = "percent.mt",
    group.by = "orig.ident"
  )
  fs_2 <- Seurat::FeatureScatter(
    seurat_obj,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA",
    group.by = "orig.ident"
  )
  fs_3 <- Seurat::FeatureScatter(
    seurat_obj,
    feature1 = "nCount_RNA",
    feature2 = "percent.rb",
    group.by = "orig.ident"
  )
  fs_4 <- Seurat::FeatureScatter(
    seurat_obj,
    feature1 = "percent.rb",
    feature2 = "percent.mt",
    group.by = "orig.ident"
  )
  fs_all <- cowplot::plot_grid(
    fs_1, fs_2, fs_3, fs_4,
    labels = c("A", "B", "C", "D"),
    ncol = 2,
    align = c("h", "v"),
    label_size = 20
  )

  pdf(
    file = paste0(seurat_out_dir, "/", project_name, "_QC_plots.pdf"),
    height = 14.4, width = 27.3
  )
  print(vp)
  print(fs_all)
  invisible(dev.off())
}

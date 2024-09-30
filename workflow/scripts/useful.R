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

method_dict <- list(
  "CCAIntegration" = c("pca", "integrated.cca"),
  "RPCAIntegration" = c("pca", "integrated.rpca"),
  "HarmonyIntegration" = c("pca", "harmony"),
  "JointPCAIntegration" = c("pca", "integrated.dr")
)

# Function to integrate Seurat layers using specified method
integrate_seurat_layers <- function(seurat_obj, layer_column, method) {
  seurat_obj <- Seurat::IntegrateLayers(
    seurat_obj,
    method = method,
    orig.reduction = method_dict[[method]][1],
    new.reduction = method_dict[[method]][2]
  )

  return(seurat_obj)
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

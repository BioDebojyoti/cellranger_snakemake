#!/opt/sw/bioinfo-tools/sources/anaconda3/envs/d_seurat510/bin/Rscript

load_library <- function(x) {
    suppressMessages(library(x, character.only = TRUE))
    return(paste0("Loaded ... ", x, collapse = ""))
}

libraries2use <- c("hdf5r", "leidenAlg", "igraph", "patchwork", "ggplot2", "Seurat", "optparse", "dplyr", "remotes", "celldex", "SingleR") # nolint

lapply(libraries2use, load_library)
# remotes::install_github("sotaro0214/seuratter")
library(seuratter)

option_list <- list(
    make_option(c("-f", "--file"), type = "character", default = "", help = "count matrix data h5 file name", dest = "count_matrix_file"), # nolint
    make_option(c("-a", "--aggregate-csv"), type = "character", default = "", help = "aggregate csv file", dest = "aggregate_csv_file"), # nolint
    make_option(c("-d", "--directory"), type = "character", default = "", help = "count matrix data directory name", dest = "count_matrix_dir"), #nolint
    make_option(c("-o", "--out"), type = "character", help = "output directory", dest = "outdir"), # nolint
    make_option(c("--min-cells"), type = "integer", default = as.integer(3), help = "minimum cells [default= %default]", dest = "min.cells"), #nolint
    make_option(c("--min-features"), type = "integer", default = as.integer(200), help = "minimum features [default= %default]", dest = "min.features"), # nolint
    make_option(c("--max-features"), type = "integer", default = as.integer(2500), help = "maximum features [default= %default]", dest = "max.features"), #nolint
    make_option(c("--percent-mt"), type = "double", default = as.double(5), help = "threshold percent mitrochondrial [default= %default]", dest = "cutoff.percent.mt"), #nolint
    make_option(c("--percent-rb"), type = "double", default = as.double(5), help = "threshold percent ribosomal [default= %default]", dest = "cutoff.percent.rb"), #nolint
    make_option(c("--project"), type = "character", default = "singleCell", help = "output file name [default= %default]", dest = "project") #nolint
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
for (p in opt){
    print(str(p)) #nolint
}

project_name = opt$project
out_dir = opt$outdir
min_cells = opt$min.cells 
min_features = opt$min.features 
max_features = opt$max.features
percent_mt = opt$cutoff.percent.mt
percent_rb = opt$cutoff.percent.rb
data_dir =  opt$count_matrix_dir
data_file = opt$count_matrix_file

data_dir = "../gosia_single_cell_runs/output_240227/count_results/aggr_results/outs/count/filtered_feature_bc_matrix"
data_file = ""
project_name = "test_seurat"
out_dir = "test_seurat"
min_cells = 5
min_features = 200
max_features = 2500 
percent_mt = 5.0
percent_rb = 30.0
aggr_csv_file = "../gosia_single_cell_runs/output_240227/count_results/aggregation_count.csv"

seurat_analysis <- function(
    data_dir, 
    data_file,
    project_name,
    out_dir,
    min_cells,
    min_features,
    max_features, 
    percent_mt, 
    percent_rb,
    aggr_csv_file
    ){

    metadata <- function(data10X, aggregate_csv_file){ #nolint
    # Load metadata
    samples <- read.csv(file.path(aggregate_csv_file), stringsAsFactors=F) #nolint

    cells <- new("seurat", raw.data=data10X)
    # Get barcodes and suffix:
    cellcodes <- as.data.frame(cells@raw.data@Dimnames[[2]]) #nolint
    colnames(cellcodes) <- "barcodes"
    rownames(cellcodes) <- cellcodes$barcodes

    cellcodes$libcodes <- as.factor(gsub(pattern=".+-", replacement="", cellcodes$barcodes)) #nolint
    cellcodes$samples <- as.vector(samples$sample_id[cellcodes$libcodes])

    # Create dataframe for meta.data argument and set up object:
    sampleidentity <- cellcodes["samples"] #nolint
    return(sampleidentity)
    } #nolint

    # Load the dataset
    if (data_dir != ""){
        data10X <- Read10X(data.dir = data_dir) # nolint
    } else {
        data10X <- Read10X_h5(filename = data_file)     # nolint
    }

    sample_identity <- metadata(data10X, aggr_csv_file)
    # Initialize the Seurat object with the raw (non-normalized data).
    seurat_obj <- Seurat::CreateSeuratObject(
        data10X,
        meta.data=sample_identity,
        min.cells=min_cells,
        min.features=min_features, 
        project= project_name)    

    # Standard pre-processing workflow
    # QC metrics
    # Low-quality / dying cells often exhibit extensive mitochondrial contamination
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-") #nolint
    seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]") #nolint

    pdf(paste0(out_dir, "/",project_name,"_violin_features.pdf"))
    VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 2)
    invisible(dev.off())

    plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") #nolint
    plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") #nolint
    plot3 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.rb") #nolint
    plot4 <- FeatureScatter(seurat_obj, feature1 = "percent.mt", feature2 = "percent.rb") #nolint

    pdf(paste0(out_dir, "/",project_name,"_features.pdf"))
    cowplot::plot_grid(plot1, plot2, plot3, plot4, labels=c("A", "B", "C", "D"), ncol=2, align=c("h","v"), label_size=20) #nolint
    invisible(dev.off())

    # filter out low quality samples
    seurat_obj <- subset(seurat_obj, nFeature_RNA > min_features & nFeature_RNA < max_features & percent.mt < percent_mt) #nolint

    # Normalization
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000) # nolint

    # highly variable feature selection
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000) # nolint

    # Identify the 10 most highly variable genes
    top10 <- head(VariableFeatures(seurat_obj), 10)

    # plot variable features with and without labels
    plot1 <- VariableFeaturePlot(seurat_obj)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    pdf(paste0(out_dir, "/",project_name,"_variable_features.pdf"))
    plot1 + plot2
    invisible(dev.off())

    # scale data
    all.genes <- rownames(seurat_obj)
    seurat_obj <- ScaleData(seurat_obj, features = all.genes) # nolint


    # linear dimensional reduction
    seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj)) # nolint
    
    # Determine the dimensionality of the dataset
    dimensionality <- determine_dimensionality(ElbowPlot(seurat_obj)$data, cutoff=0.95)

    # PCA plot
    pdf(paste0(out_dir, "/",project_name,"_pca.pdf"))
    DimPlot(seurat_obj, reduction = "pca", dims = c(1, 2))
    invisible(dev.off())

    # Cluster the cells        
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dimensionality)
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.1)    

    # Run non-linear dimensional reduction (UMAP/tSNE)
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:dimensionality)

    # UMAP plot
    pdf(paste0(out_dir, "/",project_name,"_umap.pdf"))  
    DimPlot(seurat_obj, reduction = "umap")
    invisible(dev.off())

    # celldex::
    # celldex::NovershternHematopoieticData      celldex::MouseRNAseqData                   celldex::HumanPrimaryCellAtlasData         celldex::ImmGenData                        
    # celldex::DatabaseImmuneCellExpressionData  celldex::MonacoImmuneData                  celldex::BlueprintEncodeData               
    ref <- celldex::HumanPrimaryCellAtlasData()
    input4singler <- Seurat::GetAssayData(seurat_obj)

    pred <- SingleR(test=input4singler, ref=ref, labels=ref$label.main)
    seurat_obj[["SingleR.labels"]] <- pred$labels

    b_cells <- WhichCells(seurat_obj, expression = SingleR.labels %in% c("B_cell"))
    t_cells <- WhichCells(seurat_obj, expression = SingleR.labels %in% c("T_cells"))

    p1 <- DimPlot(
        seurat_obj, 
        reduction="umap", 
        group.by="SingleR.labels", 
        split.by= "samples",
        cells.highlight=list(b_cells, t_cells), 
        cols.highlight = c("lightblue", "darkred"), cols= "grey") + 
        ggplot2::scale_color_manual(
            labels=c("unselected", "T cells","B cells"), 
            values=c("grey", "lightblue", "darkred")
            )

    # DE between each pair of cluster/annotation groups
    group_pair <- combn(unique(seurat_obj@meta.data$SingleR.labels), 2, list)    #nolint
    de_dfs <- list()
    for (p in group_pair){
        contrast = paste0(p[[1]],"_Vs_",p[[2]])                                                                                 #nolint 
        curr_de = FindMarkers(seurat_obj, ident.1 = p[[1]], ident.2=p[[2]], group.by = "SingleR.labels",logfc.threshold=log(2)) #nolint
        assign(contrast, curr_de)
        de_dfs <- list(de_dfs, contrast)
    }



    saveRDS(seurat_obj, file = paste0(out_dir, "/",project_name,"_seurat_analysis.Rds")) #nolint

    }
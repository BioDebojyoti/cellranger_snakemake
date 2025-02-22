This is a workflow for analyzing single-cell RNA-seq data. 

It starts by performing demultiplexing of bcl files to produce fastq files using
`cellranger <https://www.10xgenomics.com/support/software/cell-ranger/latest>`_ mkfastq. 
Following this fastq files are mapped to the indexed genome using cellranger count. counts from all the samples
are aggregated next using cellranger aggr. Finally, the workflow uses the `seurat <https://satijalab.org/seurat/>`_ (v5) package and does
quality control, filtering, batch correction (optional), clustering, cell-type annotation, and 
differential expression analysis (optional).
The `cellranger count <https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-gex-count>`_ 
takes FASTQ files and performs alignment, filtering, barcode counting, and UMI counting.
It uses the 10x Barcodes to generate feature-barcode matrices, determine clusters, and perform gene expression analysis. 
Web summary of `cellranger_count`_ count (alignment) step for sample: {{ snakemake.wildcards.sample }}
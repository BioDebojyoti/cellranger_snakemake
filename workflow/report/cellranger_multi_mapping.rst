The `cellranger multi <https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-5p-multi>`_ 
takes FASTQ files and performs alignment, filtering, barcode counting, and UMI counting.
It uses the 10x Barcodes to generate feature-barcode matrices, determine clusters, and perform gene expression analysis. 
Web summary of `cellranger_multi`_  multi (alignment) step for sample: {{ snakemake.wildcards.sample }}
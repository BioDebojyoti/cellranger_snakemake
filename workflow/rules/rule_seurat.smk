# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys
import pandas as pd 

# Configuration
configfile: os.path.join("config", "config_seurat.yaml")

count_matrix_file = config['count_matrix_file'] if config['count_matrix_file'] != None else if rules.


# Rule to aggregate libraries (optional)
rule seurat:
    input:
        csv = aggr_input_file if aggr_input_file != None else rules.b4all.output.aggr_input_csv
    output:
        aggr_h5 = "{count_outdir}/aggr_results/outs/count/filtered_feature_bc_matrix.h5"
    resources:
        cores = config["resources"]["localcores"],
        memory = config["resources"]["localmem"]        
    container:
        "docker://biodebojyoti/crossplatformcellranger:8.0.1"
    params:
        aggr_id = config["aggregation_id"],
        aggr_outdir = "{count_outdir}/aggr_results"
    log:
        file = "{count_outdir}/logs/aggr.log"        
    shell:
        """
        cellranger aggr --id={params.aggr_id} --csv={input.csv} --normalize=mapped \
        --localcores={resources.cores} \
        --localmem={resources.memory} \
        2>&1 | tee -a {log.file};
        bash scripts/move_pipestance_count_dir.sh {log.file} {params.aggr_outdir}; 
        """
# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys
import pandas as pd 

configfile: os.path.join("config", "config_mkvdjref.yaml")

def options2pass(genomefasta,annotationgtf, seqfasta):
    if seqfasta != None :
        options2needed = " --seq="+config["vdj_seqs"]
    else:
        options2needed = " --fasta=" + genomefasta + " --genes=" + annotationgtf 

    return options2needed

opts = options2pass(config["genome_fasta"],config["annotation_gtf"],config["vdj_seqs"])
# # Define rules
# rule all:
#     input:
#         os.path.join(config["vdjref_folder"],"fasta","{wildcards}.fa")


# Rule to count features for single library
rule cellranger_mkvdjref:
    output:
        os.path.join(config["vdjref_folder"],"fasta","*.fa")
    params:
        options_vdjref = opts,
        folder = config["vdjref_folder"]
    resources:
        cores = config["resources"]["localcores"],
        memory = config["resources"]["localmem"]        
    container:
        "docker://biodebojyoti/crossplatformcellranger:8.0.1"
    log:
        os.path.join("logs","mkvdjref.log") 
    benchmark:
        os.path.join("benchmarks", "benchmarks_mkvdjref.txt")        
    shell:
        """
        cellranger mkvdjref --genome={params.folder} {params.options_vdjref} \
        --localcores={resources.cores} --localmem={resources.memory}     
        """



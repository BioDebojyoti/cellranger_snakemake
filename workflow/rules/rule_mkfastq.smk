# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys, re
import pandas as pd 

# Configuration
configfile: os.path.join("config", "config_mkfastq.yaml")

# Define paths based on configuration
bcl_folder = config["bcl_folder"] 
fastq_folder_root = config['flowcell_id'] 
outdir = config['outputdir']

def find_keyword_line(filepath, keyword):
    with open(filepath, 'r') as file:
        for i, line in enumerate(file):
            if re.search(keyword, line):
                return i + 1  # Return the line number after the keyword
    raise ValueError(f"Keyword '{keyword}' not found in the file.")

# Load samples from CSV
if (config["IEM_samplesheet"] == True):
    keyword_line = find_keyword_line(config['samples_csv'], '(\[Cloud_Data\]|\[Data\])')
    samples_df = pd.read_csv(config['samples_csv'], skiprows=keyword_line)
    samples = samples_df['Sample_ID'].tolist()
    project = samples_df['Sample_Project'].tolist()[0]
else:
    samples_df = pd.read_csv(config['samples_csv'])
    samples = samples_df['Sample'].tolist()
    project = fastq_folder_root



opt_sample_sheet = "--samplesheet="+config["samples_csv"] if (config["IEM_samplesheet"] == True) else "--csv="+config["samples_csv"]

# Define rules
# rule all:
#     input:
#         expand(os.path.join(outdir, project,"{sample}"), sample=samples)   

# Rule to generate FASTQ files if starting from BCL files
rule mkfastq:
    input:
        os.path.abspath(bcl_folder)
    output:
        directory(os.path.join(outdir, project,"{sample}"))
    resources:
        cores = config["resources"]["localcores"],
        memory = config["resources"]["localmem"]    
    params:
        sampleinfo = opt_sample_sheet,
        outdir2use = outdir
    container:
        "docker://litd/docker-cellranger:v8.0.1" 
    log:
        os.path.join("logs","mkfastq_{sample}.log")    
    benchmark:
        os.path.join("benchmarks", "benchmarks_{sample}_mkfastq.txt")
    shell:
        """
        cellranger mkfastq --run={input} \
        --output-dir={params.outdir2use} \
        {params.sampleinfo} \
        --localcores={resources.cores} \
        --localmem={resources.memory} 
        """


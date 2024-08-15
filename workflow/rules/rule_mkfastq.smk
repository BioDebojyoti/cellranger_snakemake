# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys, re
import pandas as pd 

# Configuration
configfile: os.path.join("config", "config_mkfastq.yaml")

# Define paths based on configuration
bcl_folder = config["bcl_folder"] 

outdir = config['outputdir']
add_args = [config['additional_arguments'] if config['additional_arguments'] != None else ""]

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
    samples = samples_df['Sample_Name'].tolist()
else:
    samples_df = pd.read_csv(config['samples_csv'])
    samples = samples_df['Sample'].tolist()



opt_sample_sheet = " --samplesheet="+config["samples_csv"] if (config["IEM_samplesheet"] == True) else " --csv="+config["samples_csv"]

# Define rules
# rule all:
#     input:
#         expand(os.path.join(outdir, project,"{sample}"), sample=samples)   

# Rule to generate FASTQ files if starting from BCL files
rule mkfastq:
    input:
        os.path.abspath(bcl_folder)
    output:
        flag = os.path.join(outdir, "mkfastq.sucess"),
        all = expand(os.path.join(outdir, "{sample}_S*_R*.gz"), sample=samples)
    resources:
        cores = config["resources"]["localcores"],
        memory = config["resources"]["localmem"]    
    params:
        sampleinfo = opt_sample_sheet,
        args2add = add_args[0],
        outdir2use = lambda wc, output: os.path.splitext(os.path.dirname(output.all[0]))[0],
        file2create = lambda wc, output: os.path.abspath(output.flag)
    container:
        "docker://litd/docker-cellranger:v8.0.1" 
    log:
        os.path.join("logs","mkfastq.log")    
    benchmark:
        os.path.join("benchmarks", "benchmarks_mkfastq.txt")
    shell:
        """
        cellranger mkfastq --run={input} \
        --output-dir={params.outdir2use} \
        {params.sampleinfo} \
        {params.args2add} \
        --localcores={resources.cores} \
        --localmem={resources.memory} \
        2>&1 | tee {log} && touch {params.file2create} && \
        find {params.outdir2use} -iname *gz | grep -v "Undetermined" | xargs -I {{}} mv {{}} {params.outdir2use};
        """


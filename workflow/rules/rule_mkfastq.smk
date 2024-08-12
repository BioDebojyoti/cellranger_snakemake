# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys, re
import pandas as pd 

# def get_snakefile_dir(arguments):
#     intersection = [ arguments[ia + 1] for ia, a in enumerate(arguments) if a in ["--snakefile", "-s"]]
#     if len(intersection) == 0:
#         return os.getcwd()
#     else:
#         return os.path.dirname(os.path.abspath(intersection[0]))


# SNAKEFILE_DIR = get_snakefile_dir(sys.argv)

# Configuration
# configfile: os.path.join(SNAKEFILE_DIR, "config.yaml")
configfile: os.path.join("config", "config_mkfastq.yaml")


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
else:
    samples_df = pd.read_csv(config['samples_csv'])
    samples = samples_df['Sample'].tolist()

# Define paths based on configuration
bcl_folder = config["bcl_folder"] 
fastq_folder = os.path.join(config["flowcell_id"],"outs","fastq_path") if config["flowcell_id"] else os.path.join("*","outs","fastq_path")  

opt_sample_sheet = "--samplesheet="+config["samples_csv"] if (config["samples_csv"] != None) else ""
opt_id_args = "--id="+config["flowcell_id"] if (config["flowcell_id"] != None) else ""
# # Define rules
# rule all:
#     input:
#         expand(os.path.join(fastq_folder, "{sample}","*_R{read}_*.fastq.gz"), sample=samples, read=[1, 2])

# Rule to generate FASTQ files if starting from BCL files
rule mkfastq:
    input:
        bcl_folder
    output:
        expand(os.path.join(fastq_folder, "{sample}","*_R{read}_*.fastq.gz"), sample=samples, read=[1, 2])
    resources:
        cores = config["resources"]["localcores"],
        memory = config["resources"]["localmem"]    
    params:
        samplesheet = opt_sample_sheet,
        flowcellid = opt_id_args
    container:
        "docker://biodebojyoti/crossplatformcellranger:8.0.1"
    log:
        os.path.join("logs","mkfastq.log")    
    benchmark:
        os.path.join("benchmarks", "benchmarks_mkfastq.txt")
    shell:
        """
        cellranger mkfastq --run={input} {params.samplesheet} {params.flowcellid} \
        --localcores={resources.cores} \
        --localmem={resources.memory} 
        """


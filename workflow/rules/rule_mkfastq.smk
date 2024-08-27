# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys, re
import pandas as pd 

# Configuration
configfile: os.path.join("config", "config_mkfastq.yaml")

# Define paths based on configuration
bcl_folder = config["bcl_folder"] 

outdir = config['outputdir']
library_type = config["library_type"]

add_args = [config['additional_arguments'] if config['additional_arguments'] != None else ""]

mkfastq_cores = config["resources"]["localcores"]
mkfastq_memory = config["resources"]["localmem"]  

# Define the upper limits
max_cores = config["resources"]["max_cores"]
max_memory = config["resources"]["max_memory"]

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
    
    sample_names = samples_df['Sample_ID'].tolist()
    sample_names = list(map(str.upper, sample_names))
    
    samples_df['Sample_ID'] = sample_names
    
    lanes = samples_df['Lane'].tolist()
    lanes = list(map(str, lanes))

    samples_df["lanes"] = lanes
    samples_df["read1"] = samples_df["Sample_Name"] + "_" + samples_df["Sample_ID"] + "_L00" + samples_df["lanes"] + "_R1_001.fastq.gz"
    samples_df["read2"] = samples_df["Sample_Name"] + "_" + samples_df["Sample_ID"] + "_L00" + samples_df["lanes"] + "_R2_001.fastq.gz"
else:
    samples_df = pd.read_csv(config['samples_csv'])
    samples = samples_df['Sample'].tolist()
    sample_names = list(map(str.upper, samples))
    lanes = samples_df['Lane'].tolist()
    lanes = list(map(str, lanes))




opt_sample_sheet = " --samplesheet="+config["samples_csv"] if (config["IEM_samplesheet"] == True) else " --csv="+config["samples_csv"]

# Define rules
# rule all:
#     input:
#         expand(os.path.join(outdir, project,"{sample}"), sample=samples)   

# Rule to generate FASTQ files if starting from BCL files
rule cellranger_mkfastq:
    input:
        os.path.abspath(bcl_folder)
    output:
        flag = os.path.join(outdir, "mkfastq.sucess.csv")
    resources:
        cores = lambda wc, attempt: min(mkfastq_cores * attempt, max_cores),
        memory = lambda wc, attempt: min(mkfastq_memory * attempt, max_memory)   
    params:
        sampleinfo = opt_sample_sheet,
        args2add = add_args[0],
        outdir2use = lambda wc, output: os.path.splitext(os.path.dirname(output.flag))[0],
        file2create = lambda wc, output: os.path.abspath(output.flag),
        lib_type = library_type
    container:
        "docker://litd/docker-cellranger:v8.0.1" 
    log:
        os.path.join(outdir, "logs","mkfastq.log")    
    benchmark:
        os.path.join(outdir, "benchmarks", "benchmarks_mkfastq.csv")
    shell:
        """
        cellranger mkfastq --run={input} \
        --output-dir={params.outdir2use} \
        {params.sampleinfo} \
        {params.args2add} \
        --localcores={resources.cores} \
        --localmem={resources.memory} \
        2>&1 | tee -a {log}; 
        #find {params.outdir2use} -iname "*gz" | grep -Ev "Undetermined|\_I1_|\_I2_" | xargs -I {{}} mv {{}} {params.outdir2use};
        bash scripts/get_fastq_csv.sh {params.outdir2use} "\""{params.lib_type}"\"" > {params.file2create}; \
        bash scripts/move_pipestance_mkfastq_dir.sh {log} {params.outdir2use};
        """

# get quality control of fastq files 
# rule fastqc:
#     input:
#         fasq_path = 
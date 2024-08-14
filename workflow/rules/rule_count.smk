# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys
import pandas as pd 

# Configuration
configfile: os.path.join("config", "config_count.yaml")


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
fastq_folder = os.path.join(config["flowcell_id"],"outs","fastq_path") if config["flowcell_id"] else os.path.join("*","outs","fastq_path")  


def options4count(feature_ref, libraries_csv):
options2use = ""
    if feature_ref != None :
        options2use = options2use +  " --feature-ref " + feature_ref
    if libraries_csv != None
        options2use = options2use + " --libraries " + libraries_csv 

    return options2use

opts4count = options4count(config["feature_ref"],config["libraries_csv"])
add_arguments = config["additional_arguments"]

count_outdir = config["output_count"]
# # Define rules
rule all:
    input:
        expand(os.path.join(,"outs","filtered_feature_bc_matrix.h5"), sample=samples)


# Rule to count features for single library
rule cellranger_count:
    input:
        fastq_folder = fastq_folder,
        transcriptome = config["transcriptome"]
    output:
        matrix_h5 = os.path.join("{sample}","outs","filtered_feature_bc_matrix.h5")
    params:
        sample = lambda sample: "{sample}",
        feature_ref = config.get("feature_ref", None),
        libraries_csv = config.get("libraries_csv", None),

    resources:
        cores = config["resources"]["localcores"],
        memory = config["resources"]["localmem"]      
    
    container:
        "docker://biodebojyoti/crossplatformcellranger:8.0.1"
    log:
        os.path.join("logs","count_{sample}.log")  
    shell:
        """
        if [ "{params.use_feature_ref}" == "True" ] && [ -n "{params.feature_ref}" ]; then
        cellranger count --id={params.sample} \
            --transcriptome={input.transcriptome} \
            $libraries_param \
            $feature_ref_param \
            --create-bam={params.create_bam} \
            --localcores={resources.cores} \
            --localmem={resources.memory}; 
        else
        cellranger count --id={params.sample} \
            --transcriptome={input.transcriptome} \
            --fastqs={input.fastq_folder} \
            --sample={params.sample} \
            --create-bam={params.create_bam} \
            --localcores={resources.cores} \
            --localmem={resources.memory} \
            $expect_cells_param \
            $force_cells_param;
        fi;
        """


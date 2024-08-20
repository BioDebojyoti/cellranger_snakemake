# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys
import pandas as pd 

# Configuration
configfile: os.path.join("config", "config_count.yaml")

fastq_df = pd.read_csv(config["paths2fastq_file"])
samples_seqs = fastq_df["sample"].tolist()
samples_paths = fastq_df["fastq"].tolist()
fastq_dict = { s: samples_paths[i] for i, s in enumerate(samples_seqs)}

# def options4count(feature_ref, libraries_csv):
# options2use = ""
#     if feature_ref != None :
#         options2use = options2use +  " --feature-ref " + feature_ref
#     if libraries_csv != None
#         options2use = options2use + " --libraries " + libraries_csv 

#     return options2use

# opts4count = options4count(config["feature_ref"],config["libraries_csv"])
# add_arguments = config["additional_arguments"]

count_outdir = config["output_count"]
# # Define rules
# rule all:
#     input:
#         expand(os.path.join(count_outdir,"{sample}","outs","filtered_feature_bc_matrix.h5"), 
#         sample=list(fastq_dict.keys()))

# Rule to count features for single library
rule cellranger_count:
    input:
        fastq_folder = lambda wc: fastq_dict[wc.sample],
        transcriptome = config["transcriptome"]
    output:
        matrix_h5 = os.path.join(count_outdir, "{sample}","outs","filtered_feature_bc_matrix.h5")
    params:
        sample = lambda wc: wc.sample,
        count_outdir2use = lambda wc: os.path.join(count_outdir, wc.sample),
        add_arguments = config["additional_arguments"]
        # feature_ref = config.get("feature_ref", None),
        # libraries_csv = config.get("libraries_csv", None)
    resources:
        cores = config["resources"]["localcores"],
        memory = config["resources"]["localmem"]      
    container:
        "docker://litd/docker-cellranger:v8.0.1"      
    log:
        os.path.join("logs","count_{sample}.log")
    benchmark:
        os.path.join(outdir, "benchmarks", "benchmarks_{sample}_count.csv")             
        # if [ "{params.use_feature_ref}" == "True" ] && [ -n "{params.feature_ref}" ]; then
        # cellranger count --id={params.sample} \
        #     --transcriptome={input.transcriptome} \
        #     $libraries_param \
        #     $feature_ref_param \
        #     --create-bam={params.create_bam} \
        #     --localcores={resources.cores} \
        #     --localmem={resources.memory} \
        #     add_arguments; 
        # else
    shell:
        """        
        cellranger count --id={params.sample} \
            --transcriptome={input.transcriptome} \
            --fastqs={input.fastq_folder} \
            --sample={params.sample} \
            --output-dir={params.count_outdir2use} \
            --localcores={resources.cores} \
            --localmem={resources.memory} \
            {params.add_arguments};
        """
        # fi;


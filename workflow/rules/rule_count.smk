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
        Permoleculereadinformation = "{count_outdir}/{sample}_count/outs/molecule_info.h5"
        # summarycsv = "{count_outdir}/{sample}_count/outs/metrics_summary.csv",
        # BAM = "{count_outdir}/{sample}_count/outs/possorted_genome_bam.bam",
        # BAMindex = "{count_outdir}/{sample}_count/outs/possorted_genome_bam.bam.bai",
        # FilteredfeaturebarcodematricesHDF5 = "{count_outdir}/{sample}_count/outs/filtered_feature_bc_matrix.h5",
        # UnfilteredfeaturebarcodematricesHDF5 = "{count_outdir}/{sample}_count/outs/raw_feature_bc_matrix_h5.h5",
        # SecondaryanalysisoutputCSV = "{count_outdir}/{sample}_count/outs/analysis",
        # Permoleculereadinformation = "{count_outdir}/{sample}_count/outs/molecule_info.h5",
        # LoupeBrowserfile = "{count_outdir}/{sample}_count/outs/cloupe.cloupe" 
# - Run summary HTML:                         /home/debojyoti/Projects/core_facility/cellranger_snakemake/workflow/pipestance_14G/outs/web_summary.html
# - Run summary CSV:                          /home/debojyoti/Projects/core_facility/cellranger_snakemake/workflow/pipestance_14G/outs/metrics_summary.csv
# - BAM:                                      /home/debojyoti/Projects/core_facility/cellranger_snakemake/workflow/pipestance_14G/outs/possorted_genome_bam.bam
# - BAM BAI index:                            /home/debojyoti/Projects/core_facility/cellranger_snakemake/workflow/pipestance_14G/outs/possorted_genome_bam.bam.bai
# - BAM CSI index:                            null
# - Filtered feature-barcode matrices MEX:    /home/debojyoti/Projects/core_facility/cellranger_snakemake/workflow/pipestance_14G/outs/filtered_feature_bc_matrix
# - Filtered feature-barcode matrices HDF5:   /home/debojyoti/Projects/core_facility/cellranger_snakemake/workflow/pipestance_14G/outs/filtered_feature_bc_matrix.h5
# - Unfiltered feature-barcode matrices MEX:  /home/debojyoti/Projects/core_facility/cellranger_snakemake/workflow/pipestance_14G/outs/raw_feature_bc_matrix
# - Unfiltered feature-barcode matrices HDF5: /home/debojyoti/Projects/core_facility/cellranger_snakemake/workflow/pipestance_14G/outs/raw_feature_bc_matrix.h5
# - Secondary analysis output CSV:            /home/debojyoti/Projects/core_facility/cellranger_snakemake/workflow/pipestance_14G/outs/analysis
# - Per-molecule read information:            /home/debojyoti/Projects/core_facility/cellranger_snakemake/workflow/pipestance_14G/outs/molecule_info.h5
# - Loupe Browser file:                       /home/debojyoti/Projects/core_facility/cellranger_snakemake/workflow/pipestance_14G/outs/cloupe.cloupe        
    params:
        id2use = "pipestance_{sample}",
        sample = "{sample}",
        count_outdir2use = "{count_outdir}/{sample}_count",
        add_arguments = config["additional_arguments"],
        # feature_ref = config.get("feature_ref", None),
        # libraries_csv = config.get("libraries_csv", None)
    resources:
        cores = config["resources"]["localcores"],
        memory = config["resources"]["localmem"]      
    container:
        "docker://litd/docker-cellranger:v8.0.1"      
    log:
        file = "{count_outdir}/logs/count_{sample}.log"
    benchmark:
        "{count_outdir}/benchmarks/benchmarks_{sample}_count.csv"           
    shell:
        """   
        cellranger count \
        --id={params.id2use} \
        --transcriptome={input.transcriptome} \
        --fastqs={input.fastq_folder} \
        --sample={params.sample} \
        --localcores={resources.cores} \
        --localmem={resources.memory} \
        {params.add_arguments} \
        2>&1 | tee -a {log.file}; \
        bash scripts/move_pipestance_count_dir.sh {log.file} {params.count_outdir2use}; 
        """
        # --output-dir={params.count_outdir2use} \

        # Remove existing directory if it exists
        # if [ -d {params.count_outdir2use} ]; then
        #     rm -rf {params.count_outdir2use}
        # fi;


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
        # fi;


# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys
import pandas as pd 

def get_snakefile_dir(arguments):
    intersection = [ arguments[ia + 1] for ia, a in enumerate(arguments) if a in ["--snakefile", "-s"]]
    if len(intersection) == 0:
        return os.getcwd()
    else:
        return os.path.dirname(os.path.abspath(intersection[0]))


SNAKEFILE_DIR = get_snakefile_dir(sys.argv)

# Configuration
configfile: os.path.join(SNAKEFILE_DIR, "config.yaml")


def find_keyword_line(filepath, keyword):
    with open(filepath, 'r') as file:
        for i, line in enumerate(file):
            if keyword in line:
                return i + 1  # Return the line number after the keyword
    raise ValueError(f"Keyword '{keyword}' not found in the file.")

# Load samples from CSV
keyword_line = find_keyword_line(config['samples_csv'], 'Cloud_Data')
samples_df = pd.read_csv(config['samples_csv'], skiprows=keyword_line)
samples = samples_df['Sample_ID'].tolist()

# Define paths based on configuration
bcl_folder = config["bcl_folder"] if config["start_with_bcl"] else None
fastq_folder = os.path.join(config["flowcell_id"],"outs","fastq_path") if config["start_with_bcl"] else config["fastq_path"]

# # Define rules
rule all:
    input:
        expand(os.path.join(fastq_folder, "{sample}","*_R{read}_*.fastq.gz"), sample=samples, read=[1, 2]),
        expand(os.path.join("{sample}","outs","filtered_feature_bc_matrix.h5"), sample=samples),
        os.path.join(config["aggregation_id"],"outs","filtered_feature_bc_matrix.h5") if config.get("aggregate", False) else [],
        os.path.join("reanalyzed","outs","filtered_feature_bc_matrix.h5") if config.get("reanalyze", False) else []

# Rule to generate FASTQ files if starting from BCL files
rule mkfastq:
    input:
        bcl_folder
    output:
        fastq_folder
    resources:
        cores = config["resources"]["localcores"],
        memory = config["resources"]["localmem"]    
    params:
        samplesheet = config.get("bcl_samplesheet",None),
        flowcell_id = config.get("flowcell_id", None)
        # flowcell_id = lambda wildcards: config.get("flowcell_id", None)
    container:
        "docker://biodebojyoti/crossplatformcellranger:8.0.1"
    log:
        os.path.join("logs","mkfastq.log")    
    shell:
        """
        if [ -n "{params.samplesheet}" ]; then
        sample_sheet="--samplesheet={params.samplesheet}";
        else
        sample_sheet="";
        fi

        if [ -n "{params.flowcell_id}" ]; then
        id_arg="--id={params.flowcell_id}";
        else
        id_arg="";
        fi

        cellranger mkfastq --run {input} $id_arg --output-dir={output} \
        $sample_sheet \
        --localcores={resources.cores} \
        --localmem={resources.memory} 
        """

# Rule to count features for single library
rule cellranger_count:
    input:
        fastq_folder = fastq_folder,
        transcriptome = config["transcriptome"]
    output:
        matrix_h5 = os.path.join("{sample}","outs","filtered_feature_bc_matrix.h5")
    params:
        sample = lambda sample: "{sample}",
        expect_cells = config.get("expect_cells", None),
        pass_expect_cells = config.get("pass_expect_cells", False),
        force_cells = config.get("force_cells", None),
        pass_force_cells = config.get("pass_force_cells", False),
        feature_ref = config.get("feature_ref", None),
        use_feature_ref = config.get("use_feature_ref", False),
        libraries_csv = config.get("libraries_csv", None),
        use_libraries = config.get("use_libraries", False),
        create_bam = config.get("create_bam", False)        

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
        feature_ref_param="--feature-ref={params.feature_ref}";
        fi

        if [ "{params.use_libraries}" == "True" ] && [ -n "{params.libraries_csv}" ]; then
        libraries_param="--libraries={params.libraries_csv}";
        fi;

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

# Rule to count features for single library
rule cellranger_mkvdjref:
    input:
        genome_fasta = config.get("genome_fasta", None),
        annotation_gtf = config.get("genome_annotation_gtf", None)
    output:
        os.path.join(config["vdjref_folder"],"fasta","{wildcards}.fa")
    params:
        make_vdjref = config.get("make_vdjref", False),
        vdjref_folder = config["vdjref_folder"],
        vdj_seqs = config.get("vdj_seqs", None) 

    resources:
        cores = config["resources"]["localcores"],
        memory = config["resources"]["localmem"]        
    
    container:
        "docker://biodebojyoti/crossplatformcellranger:8.0.1"
    log:
        lambda wildcards: os.path.join("logs","mkvdjref_{wildcards}.log") 
    shell:
        """
        if [ "{params.make_vdjref}" == "True" ]; then
            if [ -n "{input.genome_fasta}" ] && [ -n "{input.annotation_gtf}" ]; then
                cellranger mkvdjref --genome={params.vdjref_folder} \
                    --fasta={input.genome_fasta} \
                    --genes={input.annotation_gtf} \
                    --localcores={resources.cores} \
                    --localmem={resources.memory}; 
            else
                cellranger mkvdjref --genome={params.vdjref_folder} \
                    --seqs={params.vdj_seqs} \
                    --localcores={resources.cores} \
                    --localmem={resources.memory}  ;       
            fi;
        fi
        """


# Rule to aggregate libraries (optional)
rule cellranger_aggr:
    input:
        csv = "aggregation.csv",
    output:
        matrix_h5 = os.path.join(config["aggregation_id"],"outs","filtered_feature_bc_matrix.h5")
    resources:
        cores = config["resources"]["localcores"],
        memory = config["resources"]["localmem"]        
    container:
        "docker://biodebojyoti/crossplatformcellranger:8.0.1"
    params:
        aggr_id = config["aggregation_id"] 
    log:
        os.path.join("logs","aggr.log")         
    shell:
        """
        cellranger aggr --id={params.aggr_id} --csv={input.csv} --normalize=mapped \
        --localcores={resources.cores} \
        --localmem={resources.memory}
        """
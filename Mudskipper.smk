from importlib.resources import path
import pandas as pd
import os

container: "docker://continuumio/miniconda3:23.10.0-1"
workdir: config['WORKDIR']
locals().update(config)
config['UNINFORMATIVE_READ'] = 3 - int(INFORMATIVE_READ)

manifest = pd.read_table(MANIFEST, index_col = False, sep = ',')
print(manifest)
barcode_df = pd.read_csv(config['barcode_csv'], header = None, sep = ':', names = ['barcode', 'RBP'])
# basic checking
assert not barcode_df['barcode'].duplicated().any()
assert not barcode_df['RBP'].duplicated().any() # cannot have any duplicated RBP names
assert not barcode_df['RBP'].str.contains(' ').any() # DO NOT CONTAIN white space lah
if READ_TYPE == "paired":
    assert not manifest['fastq1'].duplicated().any()
    assert not manifest['fastq2'].duplicated().any()
elif READ_TYPE == "single":
    assert not manifest['fastq'].duplicated().any()
assert not manifest['libname'].str.contains(' ').any()
libnames = manifest['libname'].tolist() 


config['libnames'] = libnames
experiments = manifest['experiment'].tolist()
config['experiments'] = experiments
rbps = barcode_df['RBP'].tolist()
config['rbps'] = rbps

print(f'RBPs: {rbps}',
    f'experiments:{experiments}',
    f'libnames:{libnames}')

try:
    external_normalization = config['external_bam']
    print(external_normalization)
    print('External normalization libs:',list(external_normalization.keys()))
except:
    external_normalization = None

if config['RBP_TO_RUN_MOTIF'] is None:
    config['RBP_TO_RUN_MOTIF'] = []

if config['AS_INPUT'] is None:
    config['AS_INPUT'] = []

if len(rbps)==1:
    singleplex = True
else:
    singleplex = False

# making the error files directory
try:
    os.mkdir('error_files')
except:
    pass

# making the stdout directory
try:
    os.mkdir('stdout')
except:
    pass

if READ_TYPE == "paired":
    module preprocess:
        snakefile:
            "rules/pe_preprocess.smk"
        config: config
elif READ_TYPE == "single":
    module preprocess:
        snakefile:
            "rules/se_preprocess.smk"
        config: config

module QC:
    snakefile:
        "rules/QC.smk"
    config:
        config

module DMM_BBM:
    snakefile:
        "rules/DMM_BBM.smk"
    config:
        config

module finemap:
    snakefile:
        "rules/finemap.smk"
    config:
        config

module repeat:
    snakefile:
        "rules/repeat.smk"
    config:
        config

module repeat_DMM:
    snakefile:
        "rules/repeat_DMM.smk"
    config:
        config

module bedgraphs_n_bws:
    snakefile:
        "rules/bedgraphs_n_bws.smk"
    config:
        config

module merge_bw:
    snakefile:
        "rules/merge_bw.smk"
    config:
        config

module run_homer:
    snakefile:
        "rules/run_homer.smk"
    config:
        config

module get_counts:
    snakefile:
        "rules/get_counts.smk"
    config:
        config

include: "generate_output.py"
rule all:
    input:
        get_output(config['DMM'], config['BBM'])
    
############## PREPROCESS #################
use rule * from preprocess as pre_*

############## QUALITY CONTROL #################

use rule * from QC as qc_*

if READ_TYPE == "paired":
    use rule gather_fastqc_report from QC as fastqc_gather with:
        input:
            expand("{libname}/fastqc/ultraplex_demux_{sample_label}_Rev.Tr_fastqc/fastqc_data.txt", libname = libnames, sample_label = rbps)+
            expand("{libname}/fastqc/ultraplex_demux_{sample_label}_Fwd.Tr_fastqc/fastqc_data.txt", libname = libnames, sample_label = rbps)

    use rule gather_fastqc_report from QC as fastqc_gather_initial with:
        input:
            expand("{libname}/fastqc/initial_{read}_fastqc/fastqc_data.txt",
            libname = libnames,
            read = ['r1', 'r2'])
        output:
            basic='QC/fastQC_initial_basic_summary.csv',
            passfail='QC/fastQC_initial_passfail.csv'

elif READ_TYPE == "single":
    use rule gather_fastqc_report from QC as fastqc_gather with:
        input:
            expand("{libname}/fastqc/ultraplex_demux_{sample_label}.rev_fastqc/fastqc_data.txt", libname = libnames, sample_label = rbps)
    
    use rule gather_fastqc_report from QC as fastqc_gather_initial with:
        input:
            expand("{libname}/fastqc/initial_{read}_fastqc/fastqc_data.txt", libname = libnames, read = ['r1'])
        output:
            basic='QC/fastQC_initial_basic_summary.csv',
            passfail='QC/fastQC_initial_passfail.csv' 

    use rule count_demultiplex_ultraplex from QC with:
        input:
            fq1=expand("{libname}/fastqs/ultraplex_demux_{sample_label}.fastq.gz", 
            libname = libnames, sample_label = rbps)

############## counting and finemapping #################
use rule * from get_counts as get_counts_*
use rule * from finemap as fine_*
use rule * from repeat as re_*

############## DMM and the beta binomial mixture model (BBM) #################
use rule * from DMM_BBM as dmm_bbm_*
use rule * from repeat_DMM as redmm_*

############## BEDGRAPHS & BIGWIGS #################
if READ_TYPE == "paired":
    use rule extract_read_two from bedgraphs_n_bws as extract_r1 with: # oligoPE truncation is in read1
        input:
            bam="{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam"
        output:
            read2="{libname}/bw/{sample_label}.r2.bam",
            read1="{libname}/bw/{sample_label}.r1.bam"
    
    use rule CITS_bam_to_bedgraph from bedgraphs_n_bws as CITS_bedgraph_r1 with:
        input:
            bam="{libname}/bw/{sample_label}.r1.bam"
        output:
            pos="{libname}/bw/CITS/{sample_label}.pos.bedgraph",
            neg="{libname}/bw/CITS/{sample_label}.neg.bedgraph"
    
    use rule COV_bam_to_bedgraph from bedgraphs_n_bws as COV_bedgraph_r1 with:
        input:
            bam="{libname}/bw/{sample_label}.r1.bam"
        output:
            pos="{libname}/bw/COV/{sample_label}.pos.bedgraph",
            neg="{libname}/bw/COV/{sample_label}.neg.bedgraph"
    
    use rule CITS_bam_to_bedgraph from bedgraphs_n_bws as CITS_bedgraph_external with:
        input:
            bam=lambda wildcards: ancient(config['external_bam'][wildcards.external_label]['file'])
        output:
            pos="external_bw/CITS/{external_label}.pos.bedgraph",
            neg="external_bw/CITS/{external_label}.neg.bedgraph"
        params:
            run_time="1:00:00",
            error_out_file = "error_files/coverage_bedgraph",
            out_file = "stdout/CITS_bedgraph.{external_label}",
            cores = 1,
            memory = 40000,
    use rule COV_bam_to_bedgraph from bedgraphs_n_bws as COV_bedgraph_external with:
        input:
            bam=lambda wildcards: ancient(config['external_bam'][wildcards.external_label]['file'])
        output:
            pos="external_bw/COVuse rule * from bedgraphs_n_bws/{external_label}.pos.bedgraph",
            neg="external_bw/COV/{external_label}.neg.bedgraph"
        params:
            run_time="1:00:00",
            error_out_file = "error_files/coverage_bedgraph",
            out_file = "stdout/CITS_bedgraph.{external_label}",
            cores = 1,
            memory = 40000,

    use rule * from bedgraphs_n_bws

elif READ_TYPE == "single":

    use rule CITS_bam_to_bedgraph from bedgraphs_n_bws as CITS_bedgraph with:
        input:
            bam="{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam"
        output:
            pos=temp("{libname}/bw/CITS/{sample_label}.pos.bedgraph"),
            neg=temp("{libname}/bw/CITS/{sample_label}.neg.bedgraph")
    
    use rule COV_bam_to_bedgraph from bedgraphs_n_bws as COV_bedgraph with:
        input:
            bam="{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam"
        output:
            pos=temp("{libname}/bw/COV/{sample_label}.pos.bedgraph"),
            neg=temp("{libname}/bw/COV/{sample_label}.neg.bedgraph")
    
    use rule CITS_bam_to_bedgraph from bedgraphs_n_bws as CITS_bedgraph_external with:
        input:
            bam=lambda wildcards: ancient(config['external_bam'][wildcards.external_label]['file'])
        output:
            pos=temp("external_bw/CITS/{external_label}.pos.bedgraph"),
            neg=temp("external_bw/CITS/{external_label}.neg.bedgraph")
        params:
            run_time="1:00:00",
            error_out_file = "error_files/CIT_bedgraph.{external_label}",
            out_file = "stdout/CITS_bedgraph.{external_label}",
            cores = 1,
            memory = 32000,
    
    use rule COV_bam_to_bedgraph from bedgraphs_n_bws as COV_bedgraph_external with:
        input:
            bam=lambda wildcards: ancient(config['external_bam'][wildcards.external_label]['file'])
        output:
            pos=temp("external_bw/COV/{external_label}.pos.bedgraph"),
            neg=temp("external_bw/COV/{external_label}.neg.bedgraph")
        params:
            run_time="1:00:00",
            error_out_file = "error_files/COV_bedgraph.{external_label}",
            out_file = "stdout/COV_bedgraph.{external_label}",
            cores = 1,
            memory = 32000,

    use rule bedgraph_to_bw from bedgraphs_n_bws

########## MERGE BW ############
use rule * from merge_bw

########## HOMER ############
use rule * from run_homer
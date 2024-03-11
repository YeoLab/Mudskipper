import pandas as pd

container: "docker://continuumio/miniconda3:23.10.0-1"
workdir: config['WORKDIR']

locals().update(config)

config['UNINFORMATIVE_READ'] = 3 - int(INFORMATIVE_READ) # whether read 1 or read 2 is informative

manifest = pd.read_table(MANIFEST, index_col = False, sep = ',')
print(manifest)
barcode_df = pd.read_csv(config['barcode_csv'], header = None, sep = ':', names = ['barcode', 'RBP'])
# basic checking
assert not barcode_df['barcode'].duplicated().any()
assert not barcode_df['RBP'].duplicated().any() # cannot have any duplicated RBP names
assert not barcode_df['RBP'].str.contains(' ').any() # DO NOT CONTAIN white space lah
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

module preprocess:
    snakefile:
        "rules/se_preprocess.smk"
    config: config

module QC:
    snakefile:
        "rules/QC.smk"
    config:
        config

module normalization:
    snakefile:
        "rules/skipper.smk"
    config:
        config

module DMN:
    snakefile:
        "rules/normalization_DMN.smk"
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

module repeat_dmn:
    snakefile:
        "rules/repeat_DMN.smk"
    config:
        config

module make_track:
    snakefile:
        "rules/make_track.smk"
    config:
        config

module merge_bw:
    snakefile:
        "rules/merge_bw.smk"
    config:
        config

module analysis:
    snakefile:
        "rules/analysis.smk"
    config:
        config

module clipper:
    snakefile:
        "rules/clipper.smk"
    config:
        config

module clipper_analysis:
    snakefile:
        "rules/clipper_analysis.smk"
    config:
        config


module compare:
    snakefile:
        "rules/compare.smk"
    config:
        config


include: "generate_output.py"
rule all:
    input:
        get_output(config['run_clipper'], config['run_skipper'], config['run_comparison'])
    
############## PREPROCESS #################
use rule * from preprocess as pre_*

############## QUALITY CONTROL #################
use rule * from QC as qc_*

use rule gather_fastqc_report from QC as fastqc_gather with:
    input:
        expand("{libname}/fastqc/ultraplex_demux_{sample_label}.rev_fastqc/fastqc_data.txt", 
        libname = libnames, 
        sample_label = rbps)

use rule gather_fastqc_report from QC as fastqc_gather_initial with:
    input:
        expand("{libname}/fastqc/initial_{read}_fastqc/fastqc_data.txt",
        libname = libnames,
        read = ['r1'])
    output:
        basic='QC/fastQC_initial_basic_summary.csv',
        passfail='QC/fastQC_initial_passfail.csv'

use rule count_demultiplex_ultraplex from QC with:
    input:
        fq1=expand("{libname}/fastqs/ultraplex_demux_{sample_label}.fastq.gz", 
        libname = libnames, sample_label = rbps)
    

############## SKIPPER: GENOME #################
use rule * from normalization as skipper_*
use rule * from finemap as fine_*
############## SKIPPER: REPEAT #################
use rule * from repeat as re_*

############## DMN #################
use rule * from DMN as dmn_*
use rule * from repeat_dmn as redmn_*

############## DMN #################
use rule * from clipper as clipper_*
use rule * from clipper_analysis as clipper_analysis_*

############## BIGWIGS #################
use rule CITS_bam_to_bedgraph from make_track as CITS_bedgraph with:
    input:
        bam="{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam"
    output:
        pos=temp("{libname}/bw/CITS/{sample_label}.pos.bedgraph"),
        neg=temp("{libname}/bw/CITS/{sample_label}.neg.bedgraph")

use rule COV_bam_to_bedgraph from make_track as COV_bedgraph with:
    input:
        bam="{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam"
    output:
        pos=temp("{libname}/bw/COV/{sample_label}.pos.bedgraph"),
        neg=temp("{libname}/bw/COV/{sample_label}.neg.bedgraph")

use rule CITS_bam_to_bedgraph from make_track as CITS_bedgraph_external with:
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
        memory = 40000,

use rule COV_bam_to_bedgraph from make_track as COV_bedgraph_external with:
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
        memory = 40000,

use rule bedgraph_to_bw from make_track

########## MERGE BW ############
use rule * from merge_bw

########## HOMER ############

use rule * from analysis


######## COMPARE ###########
use rule * from compare
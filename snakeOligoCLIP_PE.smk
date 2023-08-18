
from importlib.resources import path
import pandas as pd
import os
#snakemake -s snakeOligoCLIP_PE.smk -j 12 --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo -e {params.error_out_file} -o {params.out_file}" --configfile config/preprocess_config/oligope_iter7_splint.yaml --use-conda --conda-prefix /home/hsher/snakeconda -n
#snakemake -s snakeOligoCLIP_PE.smk -j 12 --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo -e {params.error_out_file} -o {params.out_file}" --configfile config/preprocess_config/oligope_iter4.yaml --use-conda --conda-prefix /home/hsher/snakeconda -n
#snakemake -s snakeOligoCLIP_PE.smk -j 30 --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo -e {params.error_out_file} -o {params.out_file}" --configfile config/preprocess_config/oligope_iter6_reseq.yaml --use-conda --conda-prefix /home/hsher/snakeconda  QC/summary.csv -n

workdir: config['WORKDIR']
MANIFEST=config['MANIFEST']
SCRIPT_PATH=config['SCRIPT_PATH']
UNINFORMATIVE_READ = 3 - int(config['INFORMATIVE_READ']) # whether read 1 or read 2 is informative
CHROM_SIZES = config['CHROM_SIZES']
R_EXE = config['R_EXE']
DB_FILE=config['DB_FILE']
GENOME_dir=config['GENOME_dir']
GENOMEFA=config['GENOMEFA']

manifest = pd.read_table(MANIFEST, index_col = False, sep = ',')
print(manifest)
barcode_df = pd.read_csv(config['barcode_csv'], header = None, sep = ':', names = ['barcode', 'RBP'])
# basic checking
assert not barcode_df['barcode'].duplicated().any()
assert not barcode_df['RBP'].duplicated().any() # cannot have any duplicated RBP names
assert not barcode_df['RBP'].str.contains(' ').any() # DO NOT CONTAIN white space lah
assert not manifest['fastq1'].duplicated().any()
assert not manifest['fastq2'].duplicated().any()
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
        "rules/pe_preprocess.smk"
    config: config

module mapr1:
    snakefile:
        "rules/map_r1.smk"
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

module repeat:
    snakefile:
        "rules/repeat.smk"
    config:
        config

module make_track:
    snakefile:
        config['MAKE_TRACK']
    config:
        config

module analysis:
    snakefile:
        "rules/analysis.smk"
    config:
        config

module repeat_dmn:
    snakefile:
        "rules/repeat_DMN.smk"
    config:
        config

# module multimap:
#     snakefile:
#         "rules/multimap.smk"
#     config:
#         config

module merge_bw:
    snakefile:
        "rules/merge_bw.smk"
    config:
        config

module finemap:
    snakefile:
        "rules/finemap.smk"
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


include: "generate_output.py"
rule all:
    input:
        get_output(config['run_clipper'], config['run_skipper'], config['run_comparison'])
    
use rule * from preprocess as pre_*
use rule * from analysis as analysis_*

############# Quality control #################
use rule * from QC as qc_*
use rule gather_fastqc_report from QC as fastqc_gather with:
    input:
        expand("{libname}/fastqc/ultraplex_demux_{sample_label}_Rev.Tr_fastqc/fastqc_data.txt", libname = libnames, sample_label = rbps)+
        expand("{libname}/fastqc/ultraplex_demux_{sample_label}_Fwd.Tr_fastqc/fastqc_data.txt", libname = libnames, sample_label = rbps)

### region caller ###
use rule * from normalization as skipper_*

### repeat caller ###
use rule * from repeat as re_*

############## DMN #################
use rule * from DMN as dmn_*
use rule * from repeat_dmn as redmn_*

############## DMN #################
use rule * from clipper as clipper_*
use rule * from clipper_analysis

########## BIGWIGS ############
use rule extract_read_two from make_track as extract_r1 with: # oligoPE truncation is in read1
    input:
        bam="{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam"
    output:
        read2="{libname}/bw/{sample_label}.r2.bam",
        read1="{libname}/bw/{sample_label}.r1.bam"

use rule strand_specific_bam from make_track as strand_specific_bam with:
    input:
        bam = "{libname}/bw/{sample_label}.r1.bam"
    output:
        pos_bam = "{libname}/bw/{sample_label}.pos.bam",
        neg_bam = "{libname}/bw/{sample_label}.neg.bam"

use rule CITS_bam_to_bedgrah from make_track as CITS_bedgraph with:
    input:
        bam="{libname}/bw/{sample_label}.r1.bam"
    output:
        pos="{libname}/bw/CITS/{sample_label}.pos.bedgraph",
        neg="{libname}/bw/CITS/{sample_label}.neg.bedgraph"
use rule COV_bam_to_bedgraph from make_track as COV_bedgraph with:
    input:
        bam="{libname}/bw/{sample_label}.r1.bam"
    output:
        pos="{libname}/bw/COV/{sample_label}.pos.bedgraph",
        neg="{libname}/bw/COV/{sample_label}.neg.bedgraph"

use rule bedgraph_to_bw from make_track
########## MERGE BW ############
use rule * from merge_bw

########## DEBUG: FILTER MULTIMAPPED ############

use rule * from finemap
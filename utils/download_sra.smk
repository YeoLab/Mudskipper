import pandas as pd
# This file downloads ABC K562_rep4/rep6 data from SRA.
# To download other replicates, see https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205536
workdir: '/tscc/nfs/home/hsher/scratch/ABC_SRA_data'

"""
snakemake -s download_sra.smk \
    -j 2 \
    --cluster "sbatch -t {params.run_time} -e {params.error_file} -o {params.out_file} -p gold -q hcg-csd792 -A csd792 --mem {params.memory} --tasks-per-node {params.cores} -J {rule}" \
    --conda-prefix /tscc/nfs/home/hsher/snakeconda \
    --use-conda \
    -n
"""
try:
    os.mkdir('error_files')
except:
    pass

try:
    os.mkdir('stdout')
except:
    pass

rule all:
    input:
        expand("{accession}.fastq.gz", accession = ['SRR19547041','SRR19547040']) # K562 files

rule get_fastq_se_gz:
    output:
        # the wildcard name must be accession, pointing to an SRA number
        "{accession}.fastq.gz",
    log:
        "error_files/{accession}.gz.log"
    params:
        extra="--skip-technical --split-3 --format=fastq",
        error_file = "error_files/smkwrap.{accession}",
        out_file = "stdout/smkwrap.{accession}",
        run_time = "3:40:00",
        cores = "1",
        memory = 160000,
    threads: 6  # defaults to 6
    resources:
        tmpdir = "~/scratch",
    wrapper:
        "v3.3.3/bio/sra-tools/fasterq-dump"


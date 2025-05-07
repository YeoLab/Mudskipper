import pandas as pd
locals().update(config)
rbps = config['rbps']
experiments = config['experiments']
libnames = config['libnames']
manifest = pd.read_table(config['MANIFEST'], index_col = False, sep = ',')


def libname_to_experiment(libname):
    return manifest.loc[manifest['libname']==libname, 'experiment'].iloc[0]
def experiment_to_libname(experiment):
    libnames = manifest.loc[manifest['experiment']==experiment, 'libname'].tolist()
    assert len(libnames)>0
    return libnames

############ counting ################
# count reads in each region for each library
# line 0 is the {libname}.{sample_label}
rule partition_bam_reads:
    input:
        chrom_size = config['CHROM_SIZES'],
        bam = "{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam",        
        region_partition = config['PARTITION'],
    output:
        counts= "counts/genome/vectors/{libname}.{sample_label}.counts",
    params:
        error_out_file = "error_files/partition_bam_reads.{libname}.{sample_label}.err",
        out_file = "stdout/partition_bam_reads.{libname}.{sample_label}.out",
        run_time = "20:00",
        cores = "1",
        memory = 10000,
        replicate_label = "{libname}.{sample_label}",
        uninformative_read = config['UNINFORMATIVE_READ']
    benchmark: "benchmarks/counts/unassigned_experiment.{libname}.{sample_label}.partition_bam_reads.txt"
    conda:
        "envs/bedtools.yaml"
    shell:
        "bedtools bamtobed -i {input.bam} | awk '($1 != \"chrEBV\") && ($4 !~ \"/{params.uninformative_read}$\")' | bedtools flank -s -l 1 -r 0 -g {input.chrom_size} -i - | bedtools shift -p -1 -m 1 -g {input.chrom_size} -i - | bedtools coverage -counts -s -a {input.region_partition} -b - | cut -f 7 | awk 'BEGIN {{print \"{params.replicate_label}\"}} {{print}}' > {output.counts};"

# concat all reps of the same experiment into 1 table with annotation
# outputs columns: [annotation] [repcounts]
rule make_genome_count_table:
    input:
        partition=config['PARTITION'],
        replicate_counts = lambda wildcards: expand("counts/genome/vectors/{libname}.{sample_label}.counts", 
            libname = experiment_to_libname(wildcards.experiment), 
            sample_label = [wildcards.sample_label]),
    output:
        count_table = "counts/genome/tables/{experiment}.{sample_label}.tsv.gz",
    params:
        error_out_file = "error_files/{experiment}.{sample_label}.make_count_table.err",
        out_file = "stdout/{experiment}.{sample_label}.make_count_table.out",
        run_time = "00:05:00",
        cores = "1",
        memory = 200,
    benchmark: "benchmarks/counts/{experiment}.{sample_label}.all_replicates.make_genome_count_table.txt"
    shell:
        "paste <(zcat {input.partition} | awk -v OFS=\"\\t\" 'BEGIN {{print \"chr\\tstart\\tend\\tname\\tscore\\tstrand\"}} {{print $1,$2,$3,$4,$5,$6}}' ) {input.replicate_counts} | gzip -c > {output.count_table}"

rule calc_partition_nuc:
    input:
        partition = config['PARTITION'],
        genome = config['GENOMEFA']
    output:
        nuc = config['PARTITION'].replace(".bed", ".nuc")
    params:
        error_out_file = "stderr/calc_partition_nuc.err",
        out_file = "stdout/calc_partition_nuc.out",
        run_time = "00:10:00",
        memory = 1000,
    benchmark: "benchmarks/partition_nuc.txt"
    container:
        "docker://howardxu520/skipper:bigwig_1.0"
    shell:
        "bedtools nuc -s -fi {input.genome} -bed {input.partition} > {output.nuc}"

####################### CC: complementary control ######################
rule sum_all_other_background_as_CC:
    input:
        lambda wildcards: expand("counts/genome/vectors/{libname}.{sample_label}.counts",
        libname = ["{libname}"],
        sample_label = list(set(rbps)-set([wildcards.clip_sample_label])-set(config['AS_INPUT']))
        )
    output:
        counts= "counts_CC/genome/vectors/{libname}.{clip_sample_label}.counts",
    params:
        error_out_file = "error_files/sum_reads.{libname}.{clip_sample_label}",
        out_file = "stdout/sum_reads.{libname}.{clip_sample_label}",
        run_time = "20:00",
        cores = "1",
        memory = 10000,
        replicate_label = "{libname}.internal"
    benchmark: "benchmarks/counts/unassigned_experiment.{libname}.{clip_sample_label}.sum_read.txt"
    shell:
        """
        awk '{{arr[FNR]+=$1}}END{{for(i=2;i<=FNR;i+=1){{print arr[i]}} }}' {input} | awk 'BEGIN {{print \"{params.replicate_label}\"}} {{print}}' > {output}
        """

rule combine_ip_to_CC:
    input:
        count_table = "counts/genome/tables/{experiment}.{clip_sample_label}.tsv.gz",
        bg_counts = lambda wildcards: expand("counts_CC/genome/vectors/{libname}.{clip_sample_label}.counts", 
            libname = experiment_to_libname(wildcards.experiment), 
            clip_sample_label = [wildcards.clip_sample_label])
    output:
        combined_count_table = "counts_CC/genome/bgtables/internal/{experiment}.{clip_sample_label}.tsv.gz"
    params:
        error_out_file = "error_files/combine.CC.{experiment}.{clip_sample_label}.err",
        out_file = "stdout/combine.CC{experiment}.{clip_sample_label}.out",
        run_time = "1:00:00",
        cores = "1",
        memory = 10000,
    benchmark: "benchmarks/combine_table/{experiment}.internal.{clip_sample_label}.combine.txt"
    shell:
        """
        paste <(zcat {input.count_table}) {input.bg_counts}| gzip -c > {output.combined_count_table}
        """


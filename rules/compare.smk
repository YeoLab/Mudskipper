import pandas as pd
import os
locals().update(config)
manifest = pd.read_table(config['MANIFEST'], index_col = False, sep = ',')
rbps = config['rbps']
experiments = config['experiments']
libnames = config['libnames']




def libname_to_experiment(libname):
    return manifest.loc[manifest['libname']==libname, 'experiment'].iloc[0]
def experiment_to_libname(experiment):
    libnames = manifest.loc[manifest['experiment']==experiment, 'libname'].tolist()
    assert len(libnames)>0
    return libnames

rule piranha_internal:
    input:
        counts = lambda wildcards: f"counts_CC/genome/bgtables/internal/"+libname_to_experiment(wildcards.libname)+f".{wildcards.clip_sample_label}.tsv.gz",
    output:
        ip = temp("comparison/piranha/CC/{libname}.{clip_sample_label}.IP.bed"),
        filtered_ip = temp("comparison/piranha/CC/{libname}.{clip_sample_label}.IP.filtered.bed"),
        cc = temp("comparison/piranha/CC/{libname}.{clip_sample_label}.CC.bed"),
        out="comparison/piranha/CC/{libname}.{clip_sample_label}.bed",
    params:
        error_out_file = "error_files/prianha.{libname}.{clip_sample_label}.err",
        out_file = "stdout/piranha.{libname}.{clip_sample_label}.out",
        run_time = "06:10:00",
        memory = 10000,
        cores = 1,
    conda:
        "envs/piranha.yaml"
    benchmark: "benchmarks/run_piranha.{libname}.{clip_sample_label}.txt"
    shell:
        """  
        zcat {input.counts} | awk -f {SCRIPT_PATH}/awk_by_colname.txt -v cols=chr,start,end,name,{wildcards.libname}.{wildcards.clip_sample_label},strand > {output.ip}
        zcat {input.counts} | awk -f {SCRIPT_PATH}/awk_by_colname.txt -v cols=chr,start,end,name,{wildcards.libname}.internal,strand > {output.cc}
        awk '{{if ($5>0) {{print}}}}' {output.ip} > {output.filtered_ip}
        Piranha -s -u 0 {output.filtered_ip} {output.cc} > {output.out};
        """

rule piranha_nobg:
    input:
        counts = lambda wildcards: f"counts_CC/genome/bgtables/internal/"+libname_to_experiment(wildcards.libname)+f".{wildcards.clip_sample_label}.tsv.gz",
    output:
        ip = temp("comparison/piranha/nobg/{libname}.{clip_sample_label}.IP.bed"),
        filtered_ip = temp("comparison/piranha/nobg/{libname}.{clip_sample_label}.IP.filtered.bed"),
        out="comparison/piranha/nobg/{libname}.{clip_sample_label}.bed",
    params:
        error_out_file = "error_files/prianha.{libname}.{clip_sample_label}.err",
        out_file = "stdout/piranha.{libname}.{clip_sample_label}.out",
        run_time = "06:10:00",
        memory = 10000,
        cores = 1,
    conda:
        "envs/piranha.yaml"
    benchmark: "benchmarks/run_piranha.{libname}.{clip_sample_label}.txt"
    shell:
        """  
        zcat {input.counts} | awk -f {SCRIPT_PATH}/awk_by_colname.txt -v cols=chr,start,end,name,{wildcards.libname}.{wildcards.clip_sample_label},strand > {output.ip}
        awk '{{if ($5>0) {{print}}}}' {output.ip} > {output.filtered_ip}
        Piranha -s -u 0 {output.filtered_ip} > {output.out};
        """

rule uniquely_mapped_reads_for_omni_and_pure:
    input:
        bam_ip="{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam",
        bam_in="CLIPper/CC_bams/{libname}.{sample_label}.sorted.bam",
    output:
        bam_ip_umap = temp("{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.umap.bam"),
        bam_in_umap = temp("CLIPper/CC_bams_umap/{libname}.{sample_label}.sorted.bam"),
    params:
        error_out_file = "error_files/uniquemap.{libname}.{sample_label}.err",
        out_file = "stdout/uniquemap.{libname}.{sample_label}.out",
        run_time = "0:30:00",
        memory = 10000,
        cores = 1,
    conda:
        "envs/bamtools.yaml"
    shell:
        """
        bamtools filter -in {input.bam_ip} -out {output.bam_ip_umap} -mapQuality ">3"
        bamtools filter -in {input.bam_in} -out {output.bam_in_umap} -mapQuality ">3"
        samtools index {output.bam_ip_umap}
        samtools index {output.bam_in_umap}
        """

# The below rules don't work
rule pureclip_internal:
    input:
        bam_ip_umap = "{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.umap.bam",
        bam_in_umap = "CLIPper/CC_bams_umap/{libname}.{sample_label}.sorted.bam",
    output:
        csln = "comparison/pureclip/{libname}.{sample_label}.csln.bed",
        bind = "comparison/pureclip/{libname}.{sample_label}.bind.bed",
        param = "comparison/pureclip/{libname}.{sample_label}.param"
    params:
        error_out_file = "error_files/pureclip.{libname}.{sample_label}.err",
        out_file = "stdout/pureclip.{libname}.{sample_label}.out",
        run_time = "60:10:00",
        memory = 10000,
        cores = 8,
    conda:
        "envs/pureclip.yaml"
    benchmark: "benchmarks/pureclip.{libname}.{sample_label}.txt"
    shell:
        """
        pureclip -i {input.bam_ip_umap} -bai {input.bam_ip_umap}.bai -g {GENOMEFA} \
            -ibam {input.bam_in_umap} -ibai {input.bam_in_umap}.bai \
            -o {output.csln} \
            -or {output.bind} \
            -p {output.param} \
            --nta 8 \
            --nt 8
        """


rule omniclip_parse:
    input:
        bam_ip_umap = "{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.umap.bam",
        bam_in_umap = "CLIPper/CC_bams_umap/{libname}.{sample_label}.sorted.bam",
    output:
        cc_data = "comparison/omniCLIP/{libname}.{sample_label}.cc.h5py",
        ip_data = "comparison/omniCLIP/{libname}.{sample_label}.ip.h5py"
    params:
        error_out_file = "error_files/omniclip_parse.{libname}.{sample_label}.err",
        out_file = "stdout/omniclip_parse.{libname}.{sample_label}.out",
        run_time = "16:10:00",
        memory = 10000,
        cores = 1,
    conda:
        "envs/omniclip.yaml"
    benchmark: "benchmarks/omniclip/parse.{libname}.{sample_label}.txt"
    shell:
        """
        omniCLIP parsingBG --db-file {DB_FILE} --bg-files {input.bam_in_umap} --genome-dir {GENOME_dir} --out-file {output.cc_data}
        omniCLIP parsingCLIP --db-file {DB_FILE} --clip-files {input.bam_ip_umap} --genome-dir {GENOME_dir} --out-file {output.ip_data}
        """

rule omniclip_call:
    input:
        cc_data = "comparison/omniCLIP/{libname}.{sample_label}.cc.h5py",
        ip_data = "comparison/omniCLIP/{libname}.{sample_label}.ip.h5py",
    output:
        outdir = directory("comparison/omniCLIP/output/{libname}.{sample_label}"),
        done = "comparison/omniCLIP/output/{libname}.{sample_label}.omniclip_done.txt"
    params:
        error_out_file = "error_files/omniclip_call.{libname}.{sample_label}.err",
        out_file = "stdout/omniclip_call.{libname}.{sample_label}.out",
        run_time = "16:10:00",
        memory = 10000,
        cores = 12,
    conda:
        "envs/omniclip.yaml"
    benchmark: "benchmarks/omniclip/call.{libname}.{sample_label}.txt"
    shell:
        """
        omniCLIP run_omniCLIP \
            --db-file {DB_FILE} --bg-dat {input.cc_data} --clip-dat {input.ip_data} \
            --out-dir {output.outdir} > {output.done}
        """

rule ctk:
    input:
        bam_ip_umap = "{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.umap.bam"
    output:
        mdtag = temp("{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.umap.md.sam.gz"),
        mutation_file = "{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.mutation.txt",
        tagbed = "{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.tag.bed.gz",
        
    params:
        error_out_file = "error_files/ctk.{libname}.{sample_label}.err",
        out_file = "stdout/ctk.{libname}.{sample_label}.out",
        run_time = "2:10:00",
        memory = 10000,
        cores = 1,
        peak = "comparison/CTK/{libname}.{sample_label}.uniq.peak.sig.bed",
        peak_bd = "comparison/CTK/{libname}.{sample_label}.uniq.peak.sig.boundary.bed",
        peak_PH = "comparison/CTK/{libname}.{sample_label}.uniq.peak.sig.halfPH.bed",
    conda:
        "envs/ctk.yaml"
    shell:
        """
        samtools fillmd {input.bam_ip_umap} {GENOMEFA} | gzip -c > {output.mdtag}
        parseAlignment.pl \
            -v --map-qual 1 \
            --min-len 18 \
            --mutation-file {output.mutation_file} \
            {output.mdtag} - | gzip -c > {output.tagbed}

        tag2peak.pl -big -ss \
            -v --valley-seeking -p 0.05 --valley-depth 0.9 \
            --multi-test --dbkey hg38 \
            {output.tagbed} \
            {params.peak} \
            --out-boundary {params.peak_bd} \
            --out-half-PH {params.peak_PH}
        """

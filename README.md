# Mudskipper: Multiplex CLIP processing pipeline from fastq.gz to binding sites and motifs
- [Link to original ABC paper](https://www.nature.com/articles/s41592-022-01708-8): use `snakeABC_SE.smk`
- Yeolab paired-end protocol: use `snakeOligoCLIP_PE.smk`

# Installation
- Main environment: Snakemake 7.3.8 and scipy:  
    - [Snakemake Installation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
    - install snakemake 7.3.8 using `rules/envs/snakemake.yaml`.
    - Snakemake 8 has different command line options that will need modification in `--profile`
- Singularity 3.11: [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/build_a_container.html).
    - If setting up Mudskipper on a server/computing cluster, it is reccomended that you ask a system admin to install singularity. This is to prevent potential permission issues. 
    - It is also possible to install singularity via conda, though this is not reccomended: [install via conda](https://anaconda.org/conda-forge/singularity)
- Download this repository by `git clone https://github.com/YeoLab/Mudskipper.git`.


# How to run. (Using ABC as an example)
1. Download data from [SRA](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205536). 
    - `utils/download_sra.smk` this script can be used to complete the download.
2. Prepare config and manifest `PATH_TO_YOUR_CONFIG`. Example inputs:
    - config file: `config/preprocess_config/oligose_k562.yaml`
    - manifest: `config/fastq_csv/ABC_2rep.csv`
    - barcode csv: `config/barcode_csv/ABC_barcode.csv`
3. Adjust profile for your cluster and computing resource:
    - see profiles/tscc2 as an example
    - for each option, [see documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html)
4. Run snakemake
    ```
    snakemake -s snakeABC_SE.smk \
        --configfile config/preprocess_config/oligose_k562_noalt.yaml \
        --profile profiles/tscc2
    ```
    - `-s`: if using YeoLab internal pair-end protocol, use `snakeOligoCLIP_PE.smk`. if you did ABC-CLIP data, use `snakeABC_SE.smk`.
    - `--configfile`: yaml file to specify your inputs, including where are the fastqs, what are the barcode, what reference genome...etc.
    - Add `-n` to the command in order to complete a dry run (sets up snakemake architecture without actually running). 
    - the rest of the options are in `--profile`. Adjust as needed. [see documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html)

Refer to the following sections to learn what to include in your configuration.

## Example config files for different use cases.
- Multiplex Example: 
    - Yeo lab internal pair-end protocol: `config/preprocess_config/oligope_iter5.yaml` 
    - ABC single-end protocol: `config/preprocess_config/ABC_2rep.yaml` 
- Singleplex Example: 
    - ABC single-end protocol: `config/preprocess_config/oligose_single_rbfox2_hek.yaml`
    - Yeo lab internal paired-end protocol: `config/preprocess_config/oligope_v5_nanos2.yaml`
    - Process 1 type of singleplex per 1 manifest.

# Options for Input files

## `MANIFEST`: a CSV file specifying FASTQ file locations and replicate information
All FASTQ files in the same manifest should have the same combinations of barcodes. For example, if you performed two singleplex experiments—one with barcode 1 and the other with barcode 2—they should be specified in two separate config files.
- Example manifest files: 
    - Multiplex Example:
        - Yeo lab paired-end: `config/fastq_csv/katie_pe_iteration5.csv`
        - ABC: `config/fastq_csv/ABC_2rep.csv`
    - Singleplex Example:
        - Yeo labe paired-end: `config/fastq_csv/V5_NANOS2.csv`
        - ABC: `/config/fastq_csv/ABC_SLBP_singleplex.csv`
- columns:
    - `fastq1`&`fastq2`: *.fastq.gz file for read1 and read 2
    - `libname`: unique names for each library. Should not contain space, special characters such as #,%,*
    - `experiment`: unique names for experiment. **Rows with the same `experiment` will be treated as replicates.** Should not contain space, special characters such as #,%,*

## `barcode_csv`: specifies the barcode sequencing information for each Antibody/RBP.
- Example: `config/barcode_csv/iter5.csv`
- Notebook to generate this file (Yeolab internal user): `utils/generate barcode-iter5.ipynb`
- delimiter: `:`
- columns:
    - 1st column: barcode sequence
        - YeoLab internal protocol: read 2 starts with this sequence. Double check with `zcat READ2.fastq.gz | grep BARCODE_SEQ`. This sequence is reverse complement to the adapter sequence (see notebook for detail)
        - ABC: read starts with this sequence.
    - 2nd column: Antibody/RBP name, Should not contain space, special characters such as #,%,*.

# Options to Control Output
- `WORKDIR`: output directory
- `RBP_TO_RUN_MOTIF`: list of RBP names to run motif analysis. Must be one of the rows in `barcode_csv`.

Although the Mudskipper paper describes the Dirichlet Multinomial Mixture (DMM) algorithm, this pipeline supports all other peak callers and background models benchmarked in the paper. Below are options for running them. For everyday use, it is recommended to set all of these options to False unless you are specifically trying to regenerate all results.
- `run_clipper`: True if you want CLIPper outputs (works, but slow)
- `run_skipper`: True if you want to run Skipper. (usually doesn't work in ABC)
- `run_comparison`: True if you want to run other peak caller such as Piranha and OmniCLIP
    - `DB_FILE` and `GENOME_dir`: are the omniclip files. check out their [github](https://github.com/philippdre/omniCLIP) on how to run them
- `debug`: True if you want to debug. This tries to blast the unmapped reads.

# Options to Choose Backgrounds
By default, if the options below are left blank, the pipeline runs the Dirichlet Multinomial Mixture (DMM) algorithm for multiplex datasets, where RBPs are explicitly compared against each other. DMM is the most effective model for multiplex datasets. There is an option to compare to 'internal control' such as a spike-in or IgG etc with a barcode as background.

Unfortunately, DMM is not compatible with singleplex datasets. Calling binding sites in singleplex datasets requires an 'external control' (see details below). Without an external control, the process will stop at the read counting stage. The best external control according to Mudskipper benchmark is an eCLIP SMInput.

The backgrounds can be applied to other peak callers such as Skipper, CLIPper, and the Beta-Binomial Mixture model.

If you want to include a background library, here's how to do it:

## "Internal control: A barcode that measures background signal, present in the same FASTQ file.
- `AS_INPUT`: If you have an IgG antibody/spike-in/bead-only control in the multiplex experiment that will serve as the normalization reference, enter its name here. It must be one of the rows in the barcode_csv file. 

## "External control": a library that is NOT in the same fastq as your oligoCLIP/ABC to serve as the background.
- specify them in `external_bam` with name of the library (first line, ex `eCLIP_SLBP_SMInput`), followed by  `file:` and `INFORMATIVE_READ`
    ```
    # For example:
    eCLIP_SLBP_SMInput: 
        file: /tscc/nfs/home/hsher/ps-yeolab5/ENCODE_k562_noalt/output/bams/dedup/genome_R2/SLBP_IN_1.genome.Aligned.sort.dedup.R2.bam
        INFORMATIVE_READ: 2
    ```
- This can be an eCLIP SMInput, total RNA-seq, IgG pull down from another experiment, bead control, spike-ins
- How to generate them? the bams must be processed with the exact same STAR index as `STAR_DIR`, and is recommended to be processed with the same/similar mapping parameters as this repo or skipper.


# Preprocessing Options:
- `adaptor_fwd`,`adaptor_rev`: adapter sequence to trim. Do not include barcode
- `tile_length`: we tile adapter sequences of this length so that indels do not interfere with trimming
- `QUALITY_CUTOFF`: default 15. cutadapt params.
- `umi_length`: Length of unique molecular identifier (UMI).
- `STAR_DIR`: directory of the STAR index.

# Annotation Options:
- skipper annotations: [follow skipper instructions](https://github.com/YeoLab/skipper#prerequisites) or generate with [skipper_utils](https://github.com/algaebrown/skipper_utils)
    - Yeolab internal users: Annotations can be found in `/tscc/projects/ps-yeolab4/software/skipper/1.0.0/bin/skipper/annotations/`.
- `CHROM_SIZES`
- `GENOMEFA`

# Output files
## Trimmed fastqs, bams, bigwigs:
These are in the `EXPERIMENT_NAME` folders. For example, in your manifest.csv, there are two experiments, "GN_1019" and "GN_1020", then, under the `GN_1019/` folder you would see the following:
1. `fastqs`: The trimmed and the demultiplexed fastqs.
    - `all.Tr.fq*` is the adapter and UMI removed ones
    - `ultraplex*RBP.fastq.gz` are the demultiplexed.
2. `bams`: bam files.
    - `*rmDup.Aligned.sortedByCoord.out.bam` is the PCR deduplicated bams. This is usually the bam you would want to use for other analysis!
    - `*Aligned.sortedByCoord.out.bam` is before deduplication.
3. `bw` contains bigwig files!
    - `CITS` calculates the number of 5' read stops per nucleotide position. 5' read stop is presumed to be the crosslinking site.
    - `COV` calculates read coverage.
4. `bw_bg` contains bigwigs of "complementary control(CC)", namely, summing all the other RBPs together.

## QC: Quality control files. 
1. A summary of all QC statistics can be found in `QC/summary.csv`.
2. `cutadapt_stat.csv` contains how many times reads contain adapter, and how many reads are too short after adapter trimming.
    - One of my favorite column is "Reads/Pairs that were too short". When reads are too short after trimming, it means they are probably adapter dimer, or you fragment it too much.
3. `fastQC*` are the summary from fastqs that has been trimmed. 
    - CLIP reads are usually bad at GC content. So it is normal to see failed GC.
    - I will always look at "adapter content" in `fastQC_passfail.csv` to make sure adapters are all gone. If this column fails, maybe you input the wrong adapter sequence, or gets contaminated by some other stuffs.
4. Read count after demultiplexing. `demux_read_count.txt`.
5. Duplication: `dup_level.csv`. This file contains the number of reads after and before PCR deduplication. If the column "unique reads" is very low, you might have amplified too much!
6. Mapping statistics `mapping_stat.csv`.
    - "Unique mapped reads%" are the percentage of reads that map to only 1 position in the genome. These were called "Usable reads" in the old eCLIP terminology.
    - "% of reads mapped to multiple loci" are the reads that map to 2~100 positions. These are mostly ribosomal RNAs. It is normal to see quite some (30-50%) multimapping reads in eCLIP
    - Then it is the reason why the rest are unmapped:
        - "% of reads mapped to multiple loci", "% of reads unmapped: too short", and "other reasons". These will help your bioinformatician debug what is going wrong. Common reasons to have lots of unmapped reads:
            - Adapter trimming is bad. The read still contains adapter sequences when they enter mapping algorithms.
            - Cells is contaminated. The baterial/fungal genomes gets sequences and does not map to human genome.
7. Read level metrices in `read_count/` folder:
    - `*cosine_similarity.csv`: Here we construct a vector, contains read counts per genomic bin for each RBP. Cosine similarity measure globally how two RBP are similar with each other. If you see big numbers, it means the two RBPs are very similar. In the Antibody barcoded CLIP paper we published, the cosine similarity is typically between 0.4~0.7. If all of the RBPs are similar, it can indicate loss of specificity.
    - `*gene_name.csv`: This contains how many read per gene for each RBP.
    - `*gene_type.csv`: contains how many read per gene/transcript type for each RBP.
    - `*region.csv`: This contains how many read per region for each RBP.
    - All of the above can help you see if there is the right enrichment and specificity for your multiplex CLIP.

## Peaks, Binding Sites
This pipeline tries to integrate multiple peak callers/binding site finders and orchestrate secondary analysis (peaks per region, motifs etc).

### Mudskipper: Mixture Modelling in `beta-mixture_*` and `DMM`.
- `DMM`(Dirichlet Multinomial Mixture) considers the distribution of reads among all RBPs without summing the rest into CC. This model detects shared binding site better than beta-mixture model, but is slower. 
- `beta-mixture*` Considers enrichment of RBP reads against complementary control (CC) or an internal library(IgG)! Compared to the DMM model, it struggles to identify shared binding sites, but it runs faster.
- The folder output/folder structure is the same for both methods: 
    - For most researchers, `*enriched_windows.tsv` will be the core output of interest. This table contains the presumed binding sites.
        - column `logLR` measures confidence. This number represents the log likelihood ratio(LR) of a window being a binding site versus not, aka, how likely is it to observe the data if it is bound, versus it being not bound. A logLR of 2 means you are about 10x more likely to observe a result like this  from binding rather than from random chance. In short, a Higher logLR means higher confidence in the binding site
        - `p_bar`,`p_raw`, `fc_raw`, `fc_bar`: Effect sizes
            - `*_bar` is the estimate from the model.
            - `*_raw` is directly calculated from counts:
            - p= (read in RBP)/(read in window)
            - fc = ((read in RBP)/(read in window))/((all reads in RBP)/total reads)
    - Secondary analysis: 
        - `feature_logistic.pdf` and `feature_ridge.pdf` contains how likely each types of regions is bound.
        - `*summary.csv`: What gene types/transcript types/region types are bound.
        - `homer/` folder contains motifs generated via [HOMER](http://homer.ucsd.edu/homer/motif/)
        - `finemapping` contains finemapped binding sites. It isn't great for splice site binding proteins.
    - Modelling outputs:
        - `fit.rda` contains everything including models of various numbers of components(K).
        - `*alpha.tsv` contains the parameter alpha/beta for beta-binomial distribution for each component.
        - `*null_alpha.tsv` contains the parameter alpha/beta for a single component beta-binomial distribution.
        - `*goodness_of_fit.pdf`: AIC, BIC, log likelihood for model selection.
        - `*label_component.csv`: The labels (bound vs not bound) for each component.
        - `*mixture_weight.tsv`: This is namely badly. What is contains is $E[z_ik]$ which is "how likely each window belong to a cluster".
        - `*weights`: The "mixture weight".

### Skipper: in `skipper*`
- `skipper_CC` models RBP versus complementary control.
- `skipper_{INTERNAL_CONTROL_NAME}` models RBP versus an internal control library, e.g. IgG.
- `skipper_external/{EXTERNAL_CONTROL_NAME}` models RBP versus an external control library, e.g. RNA-seq or SMInput.
- For skipper outputs, see [skipper's documentation](https://github.com/YeoLab/skipper)!

### CLIPper: in `CLIPper*`
- `CLIPper` only uses the IP to find local read enrichment.
- `CLIPper_CC` contains local read enrichment, "normalized to" complementary control using chi-square or fisher exact test. This is what we publised in the original paper.
- `CLIPper-{EXTERNAL_CONTROL_NAME}`: contains peaks "normalized to" an external control library. (SMInput or total RNA-seq)
- `CLIPper.{INTERNAL_CONTROL_NAME}`: contains peaks "normalized to" an internal control library. (IgG)
- `*normed.compressed.annotate.bed` is the final output. See the ENCODE pipeline for columns specification

### Comparison: `comparison/` Only if you want to run Piranha, OmniCLIP and PureCLIP.
- See their respective documentation for detail.


# For developers:
The majority of the code for the pipeline is in `rules/`:
- `se_preprocess` and `pe_preprocess` takes fastq --> trim -> demultiplex -> deduplicate -> bams
- `QC` contains rules to assemble quality control statistics, and some additional debugging rules such as investigating unmapped reads and those without barcode.
- rules for bigwig generation is in `make_track.smk`
- `merge_bw.smk` sums up bigwigs to make complementary control.
- `normalization_DMN`, `repeat_DMN` contains Mudskipper code, which does mixture model/generative clustering in genomic windows and repeat windows.
- `skipper.smk`, `repeat.smk` is entirely stolen from skipper
- `clipper.smk` runs CLIPper and the chi-square things. Stolen from ENCODE pipeline.
- `analysis.smk` and `finemap.smk`: runs finemapping, motif detection from Skipper and MudSkipper

Some rules to help you debug
- `map_r1.smk` and `multimap.smk`






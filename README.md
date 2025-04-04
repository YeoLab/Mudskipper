# Mudskipper: Multiplex CLIP processing pipeline from fastq.gz to binding sites and motifs
- [Link to original ABC paper](https://www.nature.com/articles/s41592-022-01708-8)

# Installation
There are 2 methods for installing/running mudskipper. A containerized version using singularity, and a more basic version using conda. We encourage users to adopt the Singularity version of Mudskipper because it streamlines setup by requiring only Singularity, STAR, Scipy, and Snakemake. For users without Singularity or superuser privileges, we provide an option for a purely conda based installation. 

## Using Singularity
- Download this repository using `git clone https://github.com/YeoLab/Mudskipper.git`.
- Create a conda environment to run Mudskipper in.
    - `conda create -n mudskipper_s snakemake==7.32.4 scipy==1.15.2 sra-tools==2.11.0 star==2.7.10b`
    - Note: Mudskipper depends on Snakemake version 7.3.8 to run, and is incompatible with Snakemake version > 8.
- Singularity 3.11: [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/build_a_container.html).
    - If setting up Mudskipper on a server/computing cluster, it is reccomended that you ask a system admin to install singularity. This is to prevent potential permission issues. 
    - It is also possible to install singularity via conda, though this is not reccomended: [install via conda](https://anaconda.org/conda-forge/singularity)

## Using Conda
Work in progress. 

# How to run. (Using ABC-CLIP data from original paper as an example)

## Yeo lab internal example (only for Yeo lab members)
Please see yeo_lab_internal.md (work in progress). 

## Example with singularity:
1. Activate the conda environment.
   - `conda activate mudskipper_s`
2. Download data from [SRA](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205536) using the fasterq-dump command from SRA toolkit:
   - `fasterq-dump SRR19547040 -O path/to/save/data`
   - `fasterq-dump SRR19547041 -O path/to/save/data`
   - `cd path/to/save/data`
   - `gzip SRR19547040.fastq SRR19547041.fastq` 
3. Prepare config file `example/example.yaml`. For this example, you should only have to edit two lines of the config:
    - `REPO_PATH`: This must be changed to the absolute path to the Mudskipper repository that you cloned during the installation step. 
    - `WORKDIR`: This must be changed to the absolute path to the directory where you want Mudskipper to save the results/intermediate files. 
   When running Mudskipper with your own data, you will be required to edit more of the configuration.
4. Prepare manifest file `example/example_manifest.csv`:
    - `fastq`: replace `path/to/` with the absolute path to the directory you saved the fastq files to in step 3.
5. Download/create additional annotation files (must be done inside of Mudskipper's annotations folder).
    - Download the repeat table:
        - `rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz ./repeatmasker.grch38.tsv.gz`
    - Download the genome reference:
        - `curl -L -O "https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz"`
        - `gunzip GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz`
    - Generate the STAR reference:
        - `STAR --runThreadN 8 --runMode genomeGenerate --genomeDir genome_ref --genomeFastaFiles GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta --sjdbGTFtagExonParentTranscript gencode.v38.annotation.gff3.gz`
6. Adjust profile for your cluster and computing resource:
    - An example profile is provided in `profiles/tscc2` based on the San Diego Supercomputer Center’s Triton Shared Computing Cluster [(TSCC)](https://www.sdsc.edu/systems/tscc/user_guide.html). How much this example needs to be edited will depend on the similarity between your system and TSCC. 
    - for more profile options, [see documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html)
7. Run snakemake:
    ```
    snakemake -s Mudskipper.smk \
        --configfile example/example.yaml \
        --profile profiles/tscc2
    ```
    - Add `-n` to the command in order to complete a dry run (sets up snakemake architecture without actually running). 
    - the rest of the options are in `--profile`.
      
When running Mudskipper with your own data, you will be required to edit more of the configuration. Refer to the Configuration file section to learn what to include in your configuration.

## Example with Conda:
work in progress

# Configuration file:

## Options for Input files

- `MANIFEST`: a CSV file specifying FASTQ file locations and replicate information
    - All FASTQ files in the same manifest should have the same combinations of barcodes.
    - Example: `example/example_manifest.csv`
    - columns:
        - `fastq1`&`fastq2`: *.fastq.gz file for read1 and read 2
        - `libname`: unique names for each library. Should not contain space, special characters such as #,%,*
        - `experiment`: unique names for experiment. **Rows with the same `experiment` will be treated as replicates.** Should not contain space, special characters such as #,%,*
- `barcode_csv`: specifies the barcode sequencing information for each Antibody/RBP.
    - Example: `example/example_barcode.csv`
    - delimiter: `:`
    - columns:
        - 1st column: barcode sequence
        - 2nd column: Antibody/RBP name, Should not contain space, special characters such as #,%,*.
- `WORKDIR`: output directory.
- `SCRIPT_PATH`: Path to the location of Mudskipper scripts. 

## Background Options. 
By default, the options below are left blank andthe pipeline runs the Dirichlet Multinomial Mixture (DMM) algorithm for multiplex datasets, where RBPs are explicitly compared against each other. DMM is the most effective model for multiplex datasets. There is an option to compare to 'internal control' such as a spike-in or IgG etc with a barcode as background.

Unfortunately, DMM is not compatible with singleplex datasets. Calling binding sites in singleplex datasets requires an 'external control' (see details below). Without an external control, the process will stop at the read counting stage. The best external control according to Mudskipper benchmark is an eCLIP SMInput.

If you want to include a background library, utilize one of the options below:

### "Internal control: A barcode that measures background signal, present in the same FASTQ file.
- `AS_INPUT`: If you have an IgG antibody/spike-in/bead-only control in the multiplex experiment that will serve as the normalization reference, enter its name here. It must be identical to one of the entires in the first column of the barcode_csv file. 

### "External control": a library that is NOT in the same fastq as your oligoCLIP/ABC to serve as the background.
- specify them in `external_bam` with name of the library (first line, ex `eCLIP_SLBP_SMInput`), followed by  `file:` and `INFORMATIVE_READ`
    ```
    # For example:
    eCLIP_SLBP_SMInput: 
        file: path/to/SMInput/SLBP_IN_1.genome.Aligned.sort.dedup.R2.bam
        INFORMATIVE_READ: 2
    ```
- This can be an eCLIP SMInput, total RNA-seq, IgG pull down from another experiment, bead control, or spike-ins.
- The bams used must be processed with the exact same STAR index as used in `STAR_DIR`, and it is recommended to be processed with the same/similar mapping parameters as this repo or skipper.

## Analysis options:
- `DMM`,`BBM`: specify whether to run the Dirichlet multinomial mixture (DMM) model (reccomended) or the beta-binomial mixture (BBM) model. If both options are left blank, then mudskipper just completes pre-processing without further analysis.
- `FINEMAPPING`: specify whether or not to complete fineapping analysis. Note that this is required to perform motif analysis.
- `RBP_TO_RUN_MOTIF`: list of RBP names to run motif analysis. Must be one of the RBPs listed in `barcode_csv`. Note that this requires finemapping.
- `READ_TYPE`: Either "paired" for paired end or "single" for single end. 
- `INFORMATIVE_READ`: if using single end data, enter 1. if using paired-end data, enter read (1 or 2) corresponding to crosslink site
- `SEED`: specify the seed to be used for model fitting. 

## Preprocessing Options:
- `adaptor_fwd`,`adaptor_rev`: adapter sequence to trim. Do not include barcode. Note, if using single end then only one adaptor sequence is necessary. 
- `tile_length`: we tile adapter sequences of this length so that indels do not interfere with trimming.
- `QUALITY_CUTOFF`: default 15. cutadapt params.
- `umi_pattern`: pattern of unique molecular identifier (UMI).

## Annotation Options:
Annotation data for the human reference genome [GRCh38](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/) can be found in `/annotations`. It may be necessary to download additional annotation data depending on cell type or species used (or based on user preference). 
- `REPEAT_TABLE`: Coordinates of repetitive elements, available from UCSC Genome Browser
- `GENOMEFA`: A reference genome in fasta format. 
- `STAR_DIR`: Directory of the STAR index.
- `CHROM_SIZES`: A table of chromosome sizes generated from the STAR index. 
- `PARTITION`: Gzipped BED file of windows to test.
- `FEATURE_ANNOTATIONS`: Gzipped TSV file with the following columns: chrom,start,end,name,score,strand,feature_id,feature_bin,feature_type_top,feature_types,gene_name,gene_id, transcript_ids,gene_type_top,transcript_type_top,gene_types,transcript_types.

# Output files

## Mudskipper: Mixture Modelling in `beta-binomial-mixture_* (BBM)` and `DMM`.
- `DMM/`(Dirichlet Multinomial Mixture) considers the distribution of reads among all RBPs without summing the rest into CC. This model detects shared binding site better than beta-mixture model, but is slower. 
- `beta-binomial-mixture/*` Considers enrichment of RBP reads against complementary control (CC) or an internal library(IgG). Compared to the DMM model, it struggles to identify shared binding sites, but it is much faster.
- The folder output/folder structure is the same for both methods: 
    - For most researchers, the `*enriched_windows.tsv` files will be the core output of interest, as these tables contains the presumed binding sites. As such, they are saved directly into the output folders. 
        - column `logLR` measures confidence. This number represents the log likelihood ratio(LR) of a window being a binding site versus not, aka, how likely is it to observe the data if it is bound, versus it being not bound. A logLR of 2 means you are about 10x more likely to observe a result like this  from binding rather than from random chance. In short, a Higher logLR means higher confidence in the binding site
        - `p_bar`,`p_raw`, `fc_raw`, `fc_bar`: Effect sizes
            - `*_bar` is the estimate from the model.
            - `*_raw` is directly calculated from counts:
            - p= (read in RBP)/(read in window)
            - fc = ((read in RBP)/(read in window))/((all reads in RBP)/total reads)
    - Secondary analysis: 
        - `finemapping/` contains finemapped binding sites. 
        - `homer/` folder contains motifs generated via [HOMER](http://homer.ucsd.edu/homer/motif/)
        - `*cluster_summary.csv`: A summary of each of the clusters (binding site patterns) found 
    - `intermediates/`
        - `*alpha.tsv` contains the parameter alpha/beta for beta-binomial distribution for each component.
        - `*null_alpha.tsv` contains the parameter alpha/beta for a single component beta-binomial distribution.
        - `*label_component.csv`: The labels (bound vs not bound) for each component.
        - `*mixture_weight.tsv`: Contains $E[z_ik]$ values which measure how likely each window belongs to a specific cluster.
        - `*weights`: The "mixture weight".
    - `plots/`
        - `feature_logistic.pdf` and `feature_ridge.pdf` contains how likely each types of regions is bound.
        - `*goodness_of_fit.pdf`: AIC, BIC, log likelihood for model selection.

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

# For developers:
The majority of the code for the pipeline is in `rules/`:
- `se_preprocess` and `pe_preprocess` takes fastq --> trim -> demultiplex -> deduplicate -> bams
- `QC` contains rules to assemble quality control statistics, and some additional debugging rules such as investigating unmapped reads and those without barcode.
- rules for bigwig generation is in `bedgraphs_n_bws.smk`
- `merge_bw.smk` sums up bigwigs to make complementary control.
- `DMM_BBM`, `get_counts`, `repeat_DMN` contains Mudskipper code, which does mixture model/generative clustering in genomic windows and repeat windows.
- `repeat.smk` Ggther's information on repetive elements. 
- `finemap.smk` and `run_homer.smk`: runs finemapping and motif detection.






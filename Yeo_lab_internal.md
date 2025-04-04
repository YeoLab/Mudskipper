# Yeo-lab external example (currently a work in progress). 
1. Download data from [SRA](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205536). 
    - `utils/download_sra.smk` this script can be used to complete the download.
2. Prepare config file `example/example.yaml`. For this example, you should only have to edit two lines of the config.
    - `REPO_PATH`: this must be changed to the absolute path to the Mudskipper repository that you cloned during the installation step. 
    - `WORKDIR`
   When running Mudskipper with your own data, you will be required to edit more of the configuration.
3. Adjust 
4. Adjust profile for your cluster and computing resource:
    - see profiles/tscc2 as an example.
        - Note that tscc 
    - for each option, [see documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html)
5. Run snakemake
    ```
    snakemake -s Mudskipper.smk \
        --configfile example.yaml \
        --profile profiles/tscc2
    ```
    - `--configfile`: yaml file to specify your inputs, including where are the fastqs, what are the barcode, what reference genome, is the run paired or single end,...etc.
    - Add `-n` to the command in order to complete a dry run (sets up snakemake architecture without actually running). 
    - the rest of the options are in `--profile`. Adjust as needed. [see documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html)







# These configs contain the most successful ABC samples that can be served as references
- notice STAR parameters can change the QC metrics


# Command to run
```
conda activate snakemake738
module load singularitypro
cd /tscc/nfs/home/hsher/Mudskipper/
snakemake -s snakeABC_SE.smk \
    --configfile /tscc/nfs/home/hsher/Mudskipper/config/standard_preprocess_config/oligose_k562_noalt.yaml \
    --profile profiles/tscc2 \
    -n
```
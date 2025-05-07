# Yeo-lab internal example. 
1. Create a long but lightweight interactive job:
    - e.g. `srun -N 1 -c 1 -t 48:00:00 -p gold -q hcg-csd792 -A csd792 --mem 4G --pty /bin/bash`
2. Downloading additional data is not necessary. 
3. Load the mudskipper module (replaces activating the conda environment):
    - `module load mudskipper`
4. Prepare config file `yeo_lab_internal_example/oligose_k562.yaml`. For this example, you should only have to edit two lines of the config.
    - `WORKDIR`: Please change the value of the WORKDIR to a path somewhere in your own scratch folder. 
    - Please change all instances of `/tscc/nfs/home/kflanagan/projects/Mudskipper` to the absolute path to the Mudskipper repository that you cloned during the installation step. (Note, later versions of this example will just direct you to the module installation of mudskipper)
5. Downloading additional annotation files is not necessary. 
6. Profile adjustments are not necessary. 
7. Run snakemake
    ```
    snakemake -s Mudskipper.smk \
        --configfile yeo_lab_internal_example/oligose_k562.yaml \
        --profile profiles/tscc2
    ```
    - Add `-n` to the command in order to complete a dry run (sets up snakemake architecture without actually running). 







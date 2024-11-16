conda activate clipper_helper
#install this https://github.com/kopardev/clipperhelper/tree/main and its dependency
cd /tscc/nfs/home/hsher/scratch/clipper_data
GTF=/tscc/nfs/home/hsher/gencode_coords/gencode.v38.primary_assembly.annotation.gtf
SPECIES=gencode_v38
python /tscc/nfs/home/hsher/bin/clipperhelper/util/make_custom_species_files.py --gtf $GTF --species $SPECIES
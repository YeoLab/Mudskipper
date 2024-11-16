srun --qos=hcp-csd792 --partition platinum --nodes=1 --tasks-per-node=1 --cpus-per-task=1 --mem=64G --time=5:00:00 -A csd792 --pty /bin/bash -i
cd ~/scratch/olig_PE_iter12
module load star
STAR --alignEndsType Local --genomeDir /tscc/projects/ps-yeolab3/bay001/annotations/GRCh38/star_2_7_gencode40_sjdb \
    --genomeLoad NoSharedMemory --outBAMcompression 10 --outFileNamePrefix oligoCLIP_Multiplex_1028_Rep1/bams/CSTF2.r1. \
    --winAnchorMultimapNmax 100 --outFilterMultimapNmax 100 --outFilterMultimapScoreRange 1 --outSAMmultNmax 1 --outMultimapperOrder Random \
    --outFilterScoreMin 10 --outFilterType BySJout --limitOutSJcollapsed 5000000 --outReadsUnmapped Fastx \
    --outSAMattrRGline ID:CSTF2 --outSAMattributes All --outSAMmode Full --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within --outStd Log --readFilesIn oligoCLIP_Multiplex_1028_Rep1/fastqs/ultraplex_demux_CSTF2_Rev.Tr.fastq.gz \
    --readFilesCommand zcat --runMode alignReads --runThreadN 8
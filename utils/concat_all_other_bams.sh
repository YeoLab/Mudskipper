module load samtools
OUTPUT=/tscc/nfs/home/hsher/scratch/ABC_2rep_K562_rep6.bam
SORTED_OUTPUT=/tscc/nfs/home/hsher/scratch/ABC_2rep_K562_rep6.sorted.bam
samtools merge $OUTPUT /tscc/nfs/home/hsher/ps-yeolab5/ABC_2rep_rerun/K562_rep6/bams/*rmDup*.out.bam
samtools sort $OUTPUT | samtools view -Sb - > $SORTED_OUTPUT
samtools index $SORTED_OUTPUT


OUTPUT=/tscc/nfs/home/hsher/scratch/ABC_2rep_K562_rep4.bam
SORTED_OUTPUT=/tscc/nfs/home/hsher/scratch/ABC_2rep_K562_rep4.sorted.bam
samtools merge $OUTPUT /tscc/nfs/home/hsher/ps-yeolab5/ABC_2rep_rerun/K562_rep4/bams/*rmDup*.out.bam
samtools sort $OUTPUT | samtools view -Sb - > $SORTED_OUTPUT
samtools index $SORTED_OUTPUT

# exclude RBFOX2
FILES=$(find /tscc/nfs/home/hsher/ps-yeolab5/ABC_2rep_rerun/K562_rep4/bams/ -name "*rmDup*.out.bam"  ! -regex "RBFOX2*.bam")
OUTPUT=/tscc/nfs/home/hsher/scratch/ABC_2rep_K562_rep4_no_rbfox.bam
SORTED_OUTPUT=/tscc/nfs/home/hsher/scratch/ABC_2rep_K562_rep4_no_rbfox.sorted.bam
samtools merge $OUTPUT $FILES
samtools sort $OUTPUT | samtools view -Sb - > $SORTED_OUTPUT
samtools index $SORTED_OUTPUT

FILES=$(find /tscc/nfs/home/hsher/ps-yeolab5/ABC_2rep_rerun/K562_rep4/bams/ -name "*rmDup*.out.bam"  ! -regex "RBFOX2*.bam")
OUTPUT=/tscc/nfs/home/hsher/scratch/ABC_2rep_K562_rep4_no_rbfox.bam
SORTED_OUTPUT=/tscc/nfs/home/hsher/scratch/ABC_2rep_K562_rep4_no_rbfox.sorted.bam
samtools merge $OUTPUT $FILES
samtools sort $OUTPUT | samtools view -Sb - > $SORTED_OUTPUT
samtools index $SORTED_OUTPUT
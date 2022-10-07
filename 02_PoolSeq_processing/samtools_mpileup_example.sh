PATH=$PATH:~/software/samtools-0.1.18

samtools index SRR869609_aligned_MAPQ20_rdup.bam
samtools mpileup SRR869609_aligned_MAPQ20_rdup.bam -f ~/datasets/dataset_collection/stickleback/reference_genome/stickleback_v5_assembly.fa > SRR869609.pileup

mv *.pileup pileup

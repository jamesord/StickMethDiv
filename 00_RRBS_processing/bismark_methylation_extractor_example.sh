module load UHTS/Analysis/samtools/1.10
PATH=$PATH:~/software/Bismark-0.22.3

# run Bismark methylation extractor
for i in `cat ../SRR_Acc_List_AT_filt.txt | tail -n 2`; do
bismark_methylation_extractor --comprehensive --parallel 10 --bedGraph ../mapped/${i}.sra_trimmed_bismark_bt2.bam
rm CHH_context_${i}.sra_trimmed_bismark_bt2.txt
rm CHG_context_${i}.sra_trimmed_bismark_bt2.txt

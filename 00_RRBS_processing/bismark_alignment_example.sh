module load vital-it/7
module load UHTS/Aligner/bowtie2/2.3.4.1
module load UHTS/Analysis/samtools/1.10
PATH=$PATH:~/software/Bismark-0.22.3

for i in `cat ../SRR_Acc_List_AT.txt`; do

bismark --parallel 10 --genome /storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/stickleback/reference_genome/ \
../trimmed/${i}.sra_trimmed.fq.gz

done

mv *report.txt reports

module load UHTS/Quality_control/fastqc/0.11.9
module load UHTS/Quality_control/cutadapt/2.5
module load Python/3.8.6-GCCcore-10.2.0
PATH=$PATH:/storage/homefs/jo20n766/software/TrimGalore-0.6.6

for i in `cat ../SRR_Acc_List_AT.txt | tail -n 3`; do

trim_galore -j 6 --fastqc ../raw/${i}.sra.fastq.gz

done

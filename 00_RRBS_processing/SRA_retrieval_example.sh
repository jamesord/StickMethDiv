module load vital-it/7
module load UHTS/Analysis/sratoolkit/2.9.6.1
PATH=$PATH:~/software/pigz-2.6

for i in `cat SRR_Acc_List_AT.txt`; do

prefetch ${i}
fasterq-dump ~/ncbi/public/sra/${i}.sra -e 10
rm ~/ncbi/public/sra/${i}.sra

pigz ${i}.sra.fastq -p 10

done

mkdir raw
mv *.fastq.gz raw

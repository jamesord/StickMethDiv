module load vital-it/7
module load UHTS/Analysis/sratoolkit/2.9.6.1
module load UHTS/Analysis/trimmomatic/0.36
module load UHTS/Aligner/bowtie2/2.3.4.1
module load UHTS/Analysis/samtools/1.10
module load UHTS/Analysis/sambamba/0.7.1
PATH=$PATH:~/software/pigz-2.6

for i in SRR869609; do

# fetch data from SRA
prefetch ${i}
fasterq-dump ~/ncbi/public/sra/${i}.sra -e 15
rm ~/ncbi/public/sra/${i}.sra

pigz ${i}.sra_1.fastq -p 15
pigz ${i}.sra_2.fastq -p 15

# run trimmomatic
trimmomatic PE -threads 15 -phred33  raw/${i}_1.fastq.gz raw/${i}_2.fastq.gz ${i}.sra_1_trimmed_paired.fastq.gz ${i}.sra_1_unpaired.fastq.gz ${i}.sra_2_trimmed_paired.fastq.gz ${i}.sra_2_unpaired.fastq.gz SLIDINGWINDOW:4:20
rm ${i}.sra_1_unpaired.fastq.gz ${i}.sra_2_unpaired.fastq.gz

mv ${i}.sra_1_trimmed_paired.fastq.gz ~/datasets/dataset_collection/stickleback/PRJNA204958_479509/trimmed
mv ${i}.sra_2_trimmed_paired.fastq.gz ~/datasets/dataset_collection/stickleback/PRJNA204958_479509/trimmed

# MAPPING
bowtie2 -x /storage/homefs/jo20n766/datasets/dataset_collection/stickleback/reference_genome/stickleback_v5_assembly \
-1 /storage/homefs/jo20n766/datasets/dataset_collection/stickleback/PRJNA204958_479509/trimmed/${i}.sra_1_trimmed_paired.fastq.gz \
-2 /storage/homefs/jo20n766/datasets/dataset_collection/stickleback/PRJNA204958_479509/trimmed/${i}.sra_2_trimmed_paired.fastq.gz \
-p 10 | samtools sort -@ 4 -o ${i}_aligned_sorted.bam

# quality filtering and duplicate removal
sambamba view -F "mapping_quality >= 20" -t 15 -f bam -o ${i}_aligned_MAPQ20.bam ${i}_aligned_sorted.bam
sambamba markdup -r -t 15 ${i}_aligned_MAPQ20.bam ${i}_aligned_MAPQ20_rdup.bam

# cleanup
rm ${i}.sra_1.fastq.gz
rm ${i}.sra_2.fastq.gz
rm ${i}_aligned_sorted.bam
rm ${i}_aligned_MAPQ20.bam
mv ${i}_aligned_MAPQ20_rdup.bam mapped/MAPQ20_rdup

done

#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=james.ord@vetsuisse.unibe.ch
#SBATCH --job-name="TK_popoolation_mean_SDS_methdiff_bins_190422.sh"
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G

for i in SRR869609:10 SRR7470095:12; do
i1=$(echo $i | cut -d ":" -f 1)
i2=$(echo $i | cut -d ":" -f 2)

for k in pi theta D; do

perl /storage/homefs/jo20n766/software/popoolation_1.2.2/Variance-at-position.pl --pool-size ${i2} --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
--pileup /storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/stickleback/PRJNA204958_479509/mapped/MAPQ20_rdup/pileup/${i1}.pileup \
--gtf pop_SDS_bins_190422.gtf --output ${i1}_temp.${k} --measure ${k}

sed -i "s/$/\t${i1}/" ${i1}_temp.${k}

done

paste ${i1}_temp.pi <(cut -f 4 ${i1}_temp.theta) <(cut -f 4 ${i1}_temp.D) > ${i1}_SDS_temp.div

done

cat *_SDS_temp.div > pop_SDS_bins_190422.div
rm *_temp.*

for i in SRR869609:10 SRR7470095:12; do
i1=$(echo $i | cut -d ":" -f 1)
i2=$(echo $i | cut -d ":" -f 2)

for k in pi theta D; do

perl /storage/homefs/jo20n766/software/popoolation_1.2.2/Variance-at-position.pl --pool-size ${i2} --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
--pileup /storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/stickleback/PRJNA204958_479509/mapped/MAPQ20_rdup/pileup/${i1}.pileup \
--gtf pop_methdiff_bins_190422.gtf --output ${i1}_temp.${k} --measure ${k}

sed -i "s/$/\t${i1}/" ${i1}_temp.${k}

done

paste ${i1}_temp.pi <(cut -f 4 ${i1}_temp.theta) <(cut -f 4 ${i1}_temp.D) > ${i1}_methdiff_temp.div

done

cat *_methdiff_temp.div > pop_methdiff_bins_190422.div
rm *_temp.*
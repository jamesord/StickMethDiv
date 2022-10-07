#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=james.ord@vetsuisse.unibe.ch
#SBATCH --job-name="TK_popoolation_induc_excl_200422.sh"
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G

for i in SRR869609:10 SRR7470095:12; do
i1=$(echo $i | cut -d ":" -f 1)
i2=$(echo $i | cut -d ":" -f 2)

for k in pi theta D; do

perl /storage/homefs/jo20n766/software/popoolation_1.2.2/Variance-at-position.pl --pool-size ${i2} --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
--pileup /storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/stickleback/PRJNA204958_479509/mapped/MAPQ20_rdup/pileup/${i1}.pileup \
--gtf induc_excl_200422.gtf --output ${i1}_temp.${k} --measure ${k}

sed -i "s/$/\t${i1}/" ${i1}_temp.${k}

done
done

cat *_temp.theta > temp.theta
cat *_temp.pi > temp.pi
cat *_temp.D > temp.D

paste temp.pi <(cut -f 4 temp.theta) <(cut -f 4 temp.D) > pop_induc_excl_200422.div

rm *.theta *.pi *.D *.params
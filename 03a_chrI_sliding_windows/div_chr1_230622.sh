#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=james.ord@vetsuisse.unibe.ch
#SBATCH --job-name="div_chr1_230622.sh"
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G

#grep 'chrI' /storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/stickleback/PRJNA204958_479509/mapped/MAPQ20_rdup/pileup/SRR869609.pileup | grep -v 'chrII' | grep -v 'chrIII' | grep -v 'chrIV' | grep -v 'chrIX' > SRR869609_chrI.pileup

#grep 'chrI' /storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/stickleback/PRJNA204958_479509/mapped/MAPQ20_rdup/pileup/SRR7470095.pileup | grep -v 'chrII' | grep -v 'chrIII' | grep -v 'chrIV' | grep -v 'chrIX' > SRR7470095_chrI.pileup

for i in SRR869609:10 SRR7470095:12; do
i1=$(echo $i | cut -d ":" -f 1)
i2=$(echo $i | cut -d ":" -f 2)

perl /storage/homefs/jo20n766/software/popoolation_1.2.2/Variance-sliding.pl --input ${i1}_chrI.pileup --output ${i1}_chrI.pi \
--measure pi --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --min-qual 20 --pool-size ${i2} --fastq-type sanger

perl /storage/homefs/jo20n766/software/popoolation_1.2.2/Variance-sliding.pl --input ${i1}_chrI.pileup --output ${i1}_chrI.D \
--measure D --window-size 10000 --step-size 10000 --min-count 1 --min-coverage 3 --min-qual 20 --pool-size ${i2} --fastq-type sanger \
--dissable-corrections on

sed -i "s/$/\t${i1}/" ${i1}_chrI.pi
sed -i "s/$/\t${i1}/" ${i1}_chrI.D

done

cat *.pi > chrI.pi
cat *.D > chrI.D

paste chrI.pi <(cut -f 5 chrI.D) > chrI.div

rm *.params
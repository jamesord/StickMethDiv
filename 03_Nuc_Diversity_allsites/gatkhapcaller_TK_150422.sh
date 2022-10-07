#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=james.ord@vetsuisse.unibe.ch
#SBATCH --job-name="gatkhapcaller_TK_150422.sh"
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=30G

today="$(date +'%d%m%y')"

module load vital-it/7
module load UHTS/Analysis/samtools/1.10
module load UHTS/Analysis/GenomeAnalysisTK/4.2.0.0
module load UHTS/Analysis/picard-tools/2.21.8 
module load UHTS/Analysis/HTSlib/1.10.1

# where i1 is the population and i2 is the ploidy...
for i in SRR869609:20 SRR7470095:24; do
i1=$(echo $i | cut -d ":" -f 1)
i2=$(echo $i | cut -d ":" -f 2)

# filter the bam file to get reads containing the positions of interest
samtools view -b -h -L pop_DMCs_150422.bed \
/storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/stickleback/PRJNA204958_479509/mapped/MAPQ20_rdup/${i1}_aligned_MAPQ20_rdup.bam \
> ${i1}_aligned_MAPQ20_rdup_sitesNsamples.bam

# add read groups
picard-tools AddOrReplaceReadGroups I=${i1}_aligned_MAPQ20_rdup_sitesNsamples.bam O=${i1}_aligned_MAPQ20_rdup_sitesNsamples_RG.bam RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=woof

# then of course they need to be indexed again with samtools index...
samtools index ${i1}_aligned_MAPQ20_rdup_sitesNsamples_RG.bam

# GATK HaplotypeCaller
# make sure genome is INDEXED with samtools faidx and a dictionary file has been created with CreateSequenceDictionary
# in case of problems, try re-generating the index and dictionary files

# HaplotypeCaller seems to run on 4 threads by default, so we can request 4 threads in the SLURM parameters without specifying any threads in the call to GATK

GenomeAnalysisTK HaplotypeCaller -I ${i1}_aligned_MAPQ20_rdup_sitesNsamples_RG.bam \
--R /storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/stickleback/reference_genome/stickleback_v5_assembly.fa \
--sample-ploidy ${i2} -O ${i1}_sitesNsamples_read_variants.g.vcf

# Now, the vcf files generated above contain all the SNPs called in the READS that were extracted in the first step. Therefore, I want to extract only the sites of interest:
# first bgzip and tabix the vcfs
bgzip ${i1}_sitesNsamples_read_variants.g.vcf
tabix ${i1}_sitesNsamples_read_variants.g.vcf.gz

# then filter the VCFs using a 1-based regions file of the sitesNsamples, also filtering to retain only biallelic SNPs
cut -f1,1,2 pop_DMCs_for_bcftools_150422.txt > temp.txt
bcftools view ${i1}_sitesNsamples_read_variants.g.vcf.gz -m2 -M2 -v snps --regions-file temp.txt > ${i1}_sitesNsamples_variants.g.vcf.gz

# now print out the table of biallelic SNPs
zgrep -v "^#" ${i1}_sitesNsamples_variants.g.vcf.gz | awk -v OFS='\t' '{ print $1,$2,$4,$5 } ' > ${i1}_sitesNsamples_biallelic.txt
# add the population
sed -i "s/$/\t${i1}/" ${i1}_sitesNsamples_biallelic.txt

# end loop
done

# concatentate SNP lists from the two populations
cat *biallelic.txt > sitesNsamples_biallelic_${today}.txt
# remove temporary files
rm *biallelic.txt *_variants.g.vcf.gz* *.bam* *temp.txt *.idx
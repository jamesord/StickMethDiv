module load UHTS/Analysis/BEDTools/2.29.2
PATH=$PATH:~/software/samtools-0.1.18
module load R/latest

# create joint mpileup file (assuming samtools index has already been run for each BAM file)
samtools mpileup -B \
/storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/stickleback/PRJNA204958_479509/mapped/MAPQ20_rdup/SRR869609_aligned_MAPQ20_rdup.bam \
/storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/stickleback/PRJNA204958_479509/mapped/MAPQ20_rdup/SRR7470095_aligned_MAPQ20_rdup.bam \
> WS_MAS.mpileup

# create sync file
perl /software/EcologyEvolution/popoolation2/1.201/bin/mpileup2sync.pl --fastq-type sanger --min-qual 20 --input WS_MAS.mpileup --output WS_MAS.sync

# convert the sync file into bed by duplicating the 2nd column
cat WS_MAS.sync | awk '{print $1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' > WS_MAS.sync.bed

#make a bed of our sites too. we keep everything 1-based.
cat pop_DMCs_150422.gtf | awk '{print $1"\t"$4"\t"$5}' > pop_DMCs_150422.1.bed

# run bedtools intersect
bedtools intersect -f 1 -a WS_MAS.sync.bed -b pop_DMCs_150422.1.bed > WS_MAS.sync.popDMC.bed

# remove the third column of the resulting file to convert it back into a sync file
cat WS_MAS.sync.popDMC.bed | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' > WS_MAS.popDMC.sync
rm *.bed

# calculate 'gene-wise' FST using each GTF file in the directory

for i in `ls *.gtf | rev | cut -c 12- | rev`;do

# forget the get-genewise-sync.pl script...this is way faster!
cat ${i}*.gtf | awk '{print $1"\t"$4"\t"}' > sitelocs
cat ${i}*.gtf | awk '{print $10}' | tr -d '"'| tr -d ';' > siteIDs
paste sitelocs siteIDs > sitelocs_IDs.txt
Rscript get_genewise_sync_FAST.R WS_MAS.popDMC.sync sitelocs_IDs.txt WS_MAS_pop_DMCwise.sync
rm sitelocs siteIDs sitelocs_IDs.txt

# calculate fst - this part is very quick
perl /software/EcologyEvolution/popoolation2/1.201/bin/fst-sliding.pl --min-count 2 --min-coverage 3 --pool-size 22 --min-covered-fraction 0 \
--max-coverage 1000 --window-size 1000000 --step-size 1000000 --input WS_MAS_pop_DMCwise.sync --output WS_MAS_pop_DMCs.fst
rm WS_MAS_pop_DMCwise.sync

# make a neater fst file
awk '{print $6}' WS_MAS_pop_DMCs.fst | cut -d "=" -f2 > fstvals
awk '{print $1}' WS_MAS_pop_DMCs.fst > fstcats
paste fstcats fstvals > WS_MAS_${i}.fstlite
rm fstcats fstvals

# final cleanup
rm WS_MAS_pop_DMCs.fst *.params

done

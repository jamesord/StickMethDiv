# only if BAMs have not yet been sorted (include more time if needed):
#module load UHTS/Analysis/samtools/1.10
#for i in `cat ../SRR_Acc_List_AT_filt.txt`; do
#samtools sort -O bam -o ${i}.sra_trimmed_bismark_bt2_sorted.bam ../mapped/${i}.sra_trimmed_bismark_bt2.bam
#done
#module unload UHTS/Analysis/samtools/1.10

# assumes sorted BAM files in same folder

PATH=$PATH:~/software/BS-Snper
PATH=$PATH:~/software/BS-Snper/samtools-0.1.19

for i in `cat ../SRR_Acc_List_AT_filt.txt`; do
BS-Snper.pl /storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/stickleback/PRJNA324599/BS-SNPer/${i}.sra_trimmed_bismark_bt2_sorted.bam \
--fa /storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/stickleback/reference_genome/stickleback_v5_assembly.fa \
--output /storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/stickleback/PRJNA324599/BS-SNPer/${i}.snp_result_file \
--methcg /storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/stickleback/PRJNA324599/BS-SNPer/${i}.meth_cg_result_file \
--methchg /storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/stickleback/PRJNA324599/BS-SNPer/${i}.meth_chg_result_file \
--methchh /storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/stickleback/PRJNA324599/BS-SNPer/${i}.meth_chh_result_file \
--mincover 5 >${i}.SNP.vcf 2>${i}.ERR.log
done

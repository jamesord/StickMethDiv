# REMOVE C to T and G to A SNPs from .cov.gz files
for i in `cat ../SRR_Acc_List_AT_filt.txt`; do
# get SNP IDs (in format 'chr_loc')
awk '{print $1"_"$2}' ../BS-SNPer/${i}_CtoT_GtoA_SNPs.txt > ${i}_CtoT_GtoA_SNPIDs.txt
# make new .cov file with added site IDs (in format 'chr_loc')
zcat ${i}.sra_trimmed_bismark_bt2.bismark.cov.gz | awk '{print $1"_"$2, $0}' > ${i}.sra_trimmed_bismark_bt2.bismark.SNPfilt.cov
# remove lines from the new file with site ID matching SNP ID
awk 'NR == FNR {a[$1]; next} !($1 in a)' ${i}_CtoT_GtoA_SNPIDs.txt ${i}.sra_trimmed_bismark_bt2.bismark.SNPfilt.cov > tmp && mv tmp ${i}.sra_trimmed_bismark_bt2.bismark.SNPfilt.cov
# remove site ID column to convert back to original file format; command to the right of the pipe converts the file back to tab-separated format
awk '{print $2,$3,$4,$5,$6,$7}' ${i}.sra_trimmed_bismark_bt2.bismark.SNPfilt.cov | awk -v OFS='\t' '{ $1=$1; print }' > tmp && mv tmp ${i}.sra_trimmed_bismark_bt2.bismark.SNPfilt.cov
gzip ${i}.sra_trimmed_bismark_bt2.bismark.SNPfilt.cov
# remove site ID file
rm ${i}_CtoT_GtoA_SNPIDs.txt
done

# REMOVE C to T and G to A SNPs from .bedGraph.gz files
for i in `cat ../SRR_Acc_List_AT_filt.txt`; do
# get SNP IDs (in format 'chr_loc')
awk '{print $1"_"$2}' ../BS-SNPer/${i}_CtoT_GtoA_SNPs.txt > ${i}_CtoT_GtoA_SNPIDs.txt
# make new .cov file with added site IDs (in format 'chr_loc')
zcat ${i}.sra_trimmed_bismark_bt2.bedGraph.gz | awk '{print $1"_"$3, $0}' > ${i}.sra_trimmed_bismark_bt2.SNPfilt.bedGraph
# remove lines from the new file with site ID matching SNP ID
awk 'NR == FNR {a[$1]; next} !($1 in a)' ${i}_CtoT_GtoA_SNPIDs.txt ${i}.sra_trimmed_bismark_bt2.SNPfilt.bedGraph > tmp && mv tmp ${i}.sra_trimmed_bismark_bt2.SNPfilt.bedGraph
# remove site ID column to convert back to original file format; command to the right of the pipe converts the file back to tab-separated format
awk '{print $2,$3,$4,$5}' ${i}.sra_trimmed_bismark_bt2.SNPfilt.bedGraph | awk -v OFS='\t' '{ $1=$1; print }' > tmp && mv tmp ${i}.sra_trimmed_bismark_bt2.SNPfilt.bedGraph
gzip ${i}.sra_trimmed_bismark_bt2.SNPfilt.bedGraph
# remove site ID file
rm ${i}_CtoT_GtoA_SNPIDs.txt
done

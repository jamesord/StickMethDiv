for i in `cat ../SRR_Acc_List_AT_filt.txt`; do
awk '$7=="PASS" && $4=="C" && $5=="T" {print $1, $2}' ${i}.SNP.vcf > ${i}.CtoT.txt
awk '$7=="PASS" && $4=="G" && $5=="A" {print $1, $2}' ${i}.SNP.vcf > ${i}.GtoA.txt
cat ${i}.CtoT.txt ${i}.GtoA.txt > ${i}_CtoT_GtoA_SNPs.txt
done

rm *CtoT.txt *GtoA.txt
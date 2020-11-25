#GSE136103_cirrhosis 17
#GSE112271_multiple 12
#GSE149614_HCC 12
#GSE157783_IPDCO 10
#integrated_HCC 13

pheno=HCC
sc_data=GSE149614_HCC
n_cluster=12
pop=EAS


#data path
data_path=/net/mulan/disk2/yasheng/test/LDSC_test/

gs_path=${data_path}${sc_data}/gene_list/
gc_path=${data_path}${sc_data}/gene_coord/
pf_path=${data_path}1000G_Phase3_${pop}_plinkfiles/
hmp_path=${data_path}hapmap3_snps/
ldsc_outpath=${data_path}${sc_data}/ldsc_out/${pop}/

#software
ldsc=/net/mulan/home/yasheng/comparisonProject/program/ldsc/ldsc.py
mk_annot=/net/mulan/home/yasheng/comparisonProject/program/ldsc/make_annot.py

source activate ldsc


###ldsc
for num in `seq 0 12`
do
for chr in `seq 1 22`
do

gs_file=${gs_path}gene_list_${num}.txt
gc_file=${gc_path}gene_coord_all.txt
b_ref=${pf_path}1000G.${pop}.QC.${chr}
hmp_file=${hmp_path}hm.${chr}.snp

#make annot
python2 ${mk_annot} \
		--gene-set-file ${gs_file} \
		--gene-coord-file ${gc_file} \
		--windowsize 100000 \
		--bimfile ${b_ref}.bim \
		--annot-file ${ldsc_outpath}${num}_${chr}.gz

#ldsc
python2 ${ldsc} \
		--l2 --bfile ${b_ref} \
		--ld-wind-cm 1 \
		--annot ${ldsc_outpath}${num}_${chr}.gz \
		--thin-annot \
		--out ${ldsc_outpath}${num}_${chr} \
		--print-snps ${hmp_file}
done
done

#GWAS 
GWAS_input=/net/mulan/disk2/yasheng/test/rolypoly/GWAS_data/${pheno}_${pop}/
GWAS_output=${data_path}GWAS_path/
python2 ${munge_sumstats} \
		--sumstats ${GWAS_input}rs_GWAS.txt \
		--out ${GWAS_output}GWAS_${pheno}_${pop} \
		--chunksize 500000 \
		--N 452264 \
		--merge-alleles ${data_path}w_hm3.snplist

###LDSC-cts
baseline_path=${data_path}1000G_${pop}_Phase3_baselineLD/
cts_path=${data_path}${sc_data}/
outpath=${cts_path}outcome/${pop}/
weight_path=${data_path}weights_hm3_no_hla/

python2 ${ldsc} \
--h2-cts ${GWAS_output}GWAS_${pheno}_${pop}.sumstats.gz \
--ref-ld-chr ${baseline_path}baselineLD. \
--out ${outpath}${pheno} \
--ref-ld-chr-cts ${cts_path}${sc_data}_${pop}_ldcts \
--w-ld-chr ${weight_path}weights.









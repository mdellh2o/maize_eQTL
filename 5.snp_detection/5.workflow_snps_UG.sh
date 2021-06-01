#!/bin/bash

###########################################################################
# WORKFLOW SNP CALL RIL MAIS
# temporary UnifiedGenotyper
############################################################################
REF=/home/m.miculan/reference/zea_mays/Zea_mays.AGPv4.dna.toplevel.fa
ANN=/home/m.miculan/reference/zea_mays/annotation/Zea_mays.AGPv4.37.chr.gtf

GATK=/home/m.miculan/sw/gatk/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
java=/home/m.miculan/sw/jre1.8.0_151/bin/java



# launch GATK --> UnifiedGenotyper
###################################################
# call SNPs with UnifiedGenotyper only in chr, only in parental population, divided in two subsets

cd /home/m.miculan/zea_mays_projects/eQTL/scripts/snps/unifiedgenotyper
./unifiedgenotyper_tmp_parents1.sh


# COORDINATES
###################################################
# extract coordinates from the two subsets and extend the interval (+10-10 from the SNP coordinate)
cd /home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/intervals/

egrep -v '#' ../snps.UG.parents2.vcf | cut -f 1,2 | sed 's/\t/:/' > coordinates.snp.parents.2.intervals
egrep -v '#' ../snps.UG.parents3.vcf | cut -f 1,2 | sed 's/\t/:/' > coordinates.snp.parents.3.intervals
#cat coordinates.snp.parents.2.intervals coordinates.snp.parents.3.intervals | sed 's/:/\t/' | awk 'BEGIN{OFS="\t"} {print $1,$2-10, $2+10; }' | sort -k1,1 -k2,2n | bedtools merge | sed 's/\t/:/' | sed 's/\t/-/' > coordinates.snp.parents.extended.intervals
#cat coordinates.snp.parents.2.intervals coordinates.snp.parents.3.intervals | sed 's/:/\t/' | awk 'BEGIN{OFS="\t"} {print $1,$2,$2}' | sort -k1,1 -k2,2n | bedtools merge | cut -f 1,2 | sed 's/\t/:/' | sed 's/\t/-/' > coordinates.snp.parents.intervals


# launch GATK --> UnifiedGenotyper all together with  COORDINATES - FINAL RAW SET
######################################################################################################
tmp=parents
tmp_dir=/home/m.miculan/zea_mays_projects/eQTL/tmp_dir/${tmp}
VCF=/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/raw.snps.UG.${tmp}.vcf

cd /home/m.miculan/zea_mays_projects/eQTL
mkdir -p ${tmp_dir}
cd ${WORKDIR}

${java} -Xmx5g -Djava.io.tmpdir=${tmp_dir} -jar ${GATK} -T UnifiedGenotyper -R ${REF} -o ${VCF} \
-I mappings/A632_repl1/A632_repl1.realigned_L.bam \
-I mappings/A632_repl2/A632_repl2.realigned_L.bam \
-I mappings/B73_repl1/B73_repl1.realigned_L.bam \
-I mappings/B73_repl2/B73_repl2.realigned_L.bam \
-I mappings/B73_repl3/B73_repl3.realigned_L.bam \
-I mappings/CML91_repl1/CML91_repl1.realigned_L.bam \
-I mappings/CML91_repl2/CML91_repl2.realigned_L.bam \
-I mappings/F7_repl1/F7_repl1.realigned_L.bam \
-I mappings/F7_repl2/F7_repl2.realigned_L.bam \
-I mappings/H99_repl1/H99_repl1.realigned_L.bam \
-I mappings/H99_repl2/H99_repl2.realigned_L.bam \
-I mappings/HP301_repl1/HP301_repl1.realigned_L.bam \
-I mappings/HP301_repl2/HP301_repl2.realigned_L.bam \
-I mappings/Mo17_repl1/Mo17_repl1.realigned_L.bam \
-I mappings/Mo17_repl2/Mo17_repl2.realigned_L.bam \
-I mappings/W153R_repl1/W153R_repl1.realigned_L.bam \
-I mappings/W153R_repl2/W153R_repl2.realigned_L.bam \
--genotype_likelihoods_model SNP \
--heterozygosity 0.01 \
--reference_sample_name Zea_mays.AGPv4.dna.toplevel \
-L /home/m.miculan/zea_mays_projects/eQTL/SNPs/intervals/unifiedgenotyper/coordinates.snp.parents.intervals \
-dcov 2000 -nt ${NPROC} >${WORKDIR}/scripts/logs/unifiedgenotyper/${tmp}.out 2>${WORKDIR}/scripts/logs/unifiedgenotyper/${tmp}.err

rm -rf ${tmp_dir}


# launch GATK --> UnifiedGenotyper IN DIFFERENT POPULATIONS
######################################################################################################
cd /home/m.miculan/zea_mays_projects/eQTL/scripts/snps

./select_and_call_multiparental.sh 
./select_and_call_biparental.sh

# in these scripts are selected the coordinates for biparental and multiparental populations. SNPs called in subsets (for cluster limits) and then put together
# output: two raw.snp, one for each population
# /home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/raw.snps.multipop.UG.vcf 
# TOT= 4,651,361
# /home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/raw.snps.bipop.UG.vcf 
# TOT= 1,181,182

cd /home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper
mkdir vcf_tmp
mv raw.snps.multi*tmp* vcf_tmp/



######################################################################################################
# FILTER THE FINAL DATASET MULTIPARENTAL
######################################################################################################
# biallelic snps, good quality. at least 80% genotyped
cd /home/m.miculan/zea_mays_projects/eQTL/SNPs
VCF_input=/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/raw.snps.multipop.UG.vcf
VCF_output=/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/snps.multipop.UG.goodQUAL.08.vcf
threads=16
tmp_dir=/home/m.miculan/zea_mays_projects/eQTL/tmp_dir/select; mkdir ${tmp_dir}

${java} -Xmx5g -Djava.io.tmpdir=${tmp_dir} -jar ${GATK} -T SelectVariants -R ${REF} -V ${VCF_input} -o ${VCF_output} -mvq 100 -restrictAllelesTo BIALLELIC --maxNOCALLfraction 0.2 -nt ${threads}; rm -rf ${tmp_dir}
#TOT= 525,513

# filter, output=genotypes.chr.multi.filtered.txt
cd /home/m.miculan/zea_mays_projects/eQTL/scripts/snps
Rscript filter_snps_multi.r
# TOT= 377,783


# create vcf only selected (FOR MANEL)
cd /home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper
VCF_input=/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/snps.multipop.UG.goodQUAL.08.vcf
VCF_output=/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/snps.multipop.UG.goodQUAL.08.filtered.vcf
threads=16
tmp_dir=/home/m.miculan/zea_mays_projects/eQTL/tmp_dir/select; mkdir ${tmp_dir}

awk 'BEGIN{OFS="\t"} { split($1,p, "_"); print p[1], p[2]}' genotypes.chr.multi.filtered.txt | sed 's/\t/:/' | grep -v 'SNP' > genotypes.chr.multi.filtered.intervals
${java} -Xmx5g -Djava.io.tmpdir=${tmp_dir} -jar ${GATK} -T SelectVariants -R ${REF} -V ${VCF_input} -o ${VCF_output} -L genotypes.chr.multi.filtered.intervals -nt ${threads}; rm -rf ${tmp_dir}


# subsample MULTIPARENTAL
cd /home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper
VCF_input=snps.multipop.UG.goodQUAL.08.filtered.vcf
VCF_output=snps.multipop.UG.subset.vcf
tmp_dir=/home/m.miculan/zea_mays_projects/eQTL/tmp_dir/select; mkdir -p ${tmp_dir}
${java} -Xmx5g -Djava.io.tmpdir=${tmp_dir} -jar ${GATK} -T SelectVariants -R ${REF} -V ${VCF_input} -o ${VCF_output} --maxNOCALLfraction 0.05 -fraction 0.2; rm -rf ${tmp_dir}



# VALIDATE MULTIPARENTAL (to have the genotype) remove FP
# at the same time EXTRACT only BIPARENTAL POPULATION and filte 
######################################################################################################
cd /home/m.miculan/zea_mays_projects/eQTL/scripts/snps

./select_and_call_multi_in_biparental.sh

Rscript ./filter_snps_all.r



# EXTRACT to have VCF file of ALL POPULATIONS TOGETHER 
######################################################################################################
cd /home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper
VCF_input=/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/vcf_tmp/snps.all.UG.goodQUAL.08.filtered.tmp2.vcf
VCF_output=/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/snps.all.UG.goodQUAL.08.filtered.biallelic.vcf
geno_file=/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/genotypes.chr.all.filtered.txt
tmp_dir=/home/m.miculan/zea_mays_projects/eQTL/tmp_dir/select

awk 'BEGIN{OFS="\t"} { split($1,p, "_"); print p[1], p[2]}' ${geno_file} | sed 's/\t/:/' | grep -v 'SNP' > intervals/genotypes.chr.all.filtered.intervals
mkdir -p ${tmp_dir}
${java} -Xmx5g -Djava.io.tmpdir=${tmp_dir} -jar ${GATK} -T SelectVariants -R ${REF} -V ${VCF_input} -o ${VCF_output} -L intervals/genotypes.chr.all.filtered.intervals; rm -rf ${tmp_dir}


# subsample ALLPOP dataset
cd /home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper
VCF_input=snps.all.UG.goodQUAL.08.filtered.biallelic.vcf
VCF_output=snps.all.UG.subset.vcf
tmp_dir=/home/m.miculan/zea_mays_projects/eQTL/tmp_dir/select; mkdir -p ${tmp_dir}
${java} -Xmx5g -Djava.io.tmpdir=${tmp_dir} -jar ${GATK} -T SelectVariants -R ${REF} -V ${VCF_input} -o ${VCF_output} --maxNOCALLfraction 0.05 -fraction 0.19; rm -rf ${tmp_dir}



# NEW!!!
# EXTRACT  BIPARENTAL POPULATION (to have the genotype)
######################################################################################################
cd /home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper
VCF_input=/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/vcf_tmp/snps.all.UG.goodQUAL.08.filtered.tmp2.vcf
VCF_output=/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/snps.bipop.UG.goodQUAL.08.filtered2.vcf
tmp_dir=/home/m.miculan/zea_mays_projects/eQTL/tmp_dir/select; mkdir -p ${tmp_dir}

awk 'BEGIN{OFS="\t"} { split($1,p, "_"); print p[1], p[2]}' genotypes.chr.bi.filtered2.txt | sed 's/\t/:/' | grep -v 'SNP' > intervals/genotypes.chr.bi.filtered2.intervals
${java} -Xmx5g -Djava.io.tmpdir=${tmp_dir} -jar ${GATK} -T SelectVariants -R ${REF} -V ${VCF_input} -o ${VCF_output} -L intervals/genotypes.chr.bi.filtered2.intervals -sn "B73" -sn "H99" -se "_RIL" ; rm -rf ${tmp_dir}

#TOT=173,541

# subsample
VCF_input=snps.bipop.UG.goodQUAL.08.filtered2.vcf
VCF_output=snps.bipop.UG.subset.vcf
tmp_dir=/home/m.miculan/zea_mays_projects/eQTL/tmp_dir/select; mkdir -p ${tmp_dir}
${java} -Xmx5g -Djava.io.tmpdir=${tmp_dir} -jar ${GATK} -T SelectVariants -R ${REF} -V ${VCF_input} -o ${VCF_output} --maxNOCALLfraction 0.05 -fraction 0.36; rm -rf ${tmp_dir}



# FILTER THE FINAL DATASET BIPARENTAL
######################################################################################################
## biallelic snps, good quality. at least 80% genotyped
#cd /home/m.miculan/zea_mays_projects/eQTL/SNPs
#VCF_input=/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/raw.snps.bipop.UG.vcf
#VCF_output=/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/snps.bipop.UG.goodQUAL.08.vcf
#threads=16
#tmp_dir=/home/m.miculan/zea_mays_projects/eQTL/tmp_dir/select; mkdir ${tmp_dir}
#
#${java} -Xmx5g -Djava.io.tmpdir=${tmp_dir} -jar ${GATK} -T SelectVariants -R ${REF} -V ${VCF_input} -o ${VCF_output} -mvq 100 -restrictAllelesTo BIALLELIC --maxNOCALLfraction 0.2 -nt ${threads}; rm -rf ${tmp_dir}
## TOT = 247,703
#
#
## filter, output= genotypes.chr.bi.filtered.txt
#Rscript filter_snps_bi.r
## TOT= 219,965
#
#
## create vcf only selected (FOR MANEL)
#cd /home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper
#VCF_input=/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/snps.bipop.UG.goodQUAL.08.vcf
#VCF_output=/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/snps.bipop.UG.goodQUAL.08.filtered.vcf
#threads=16
#tmp_dir=/home/m.miculan/zea_mays_projects/eQTL/tmp_dir/select; mkdir ${tmp_dir}
#
#awk 'BEGIN{OFS="\t"} { split($1,p, "_"); print p[1], p[2]}' genotypes.chr.bi.filtered.txt | sed 's/\t/:/' | grep -v 'SNP' > genotypes.chr.bi.filtered.intervals
#${java} -Xmx5g -Djava.io.tmpdir=${tmp_dir} -jar ${GATK} -T SelectVariants -R ${REF} -V ${VCF_input} -o ${VCF_output} -L genotypes.chr.bi.filtered.intervals -nt ${threads}; rm -rf ${tmp_dir}





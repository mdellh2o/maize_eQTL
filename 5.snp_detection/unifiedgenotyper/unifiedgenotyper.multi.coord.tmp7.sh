#!/bin/bash

#######################################################################################################################################################################################
# SNP calling ZEA MAYS TRANSCRIPTOME
# may 2018
################################################################################################################################################################################

GATK=/home/m.miculan/sw/gatk/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
#GATK-3.3=/home/m.miculan/sw/gatk/GenomeAnalysisTK-3.3.0/GenomeAnalysisTK.jar
java=/home/m.miculan/sw/jre1.8.0_151/bin/java

tmp=multi.tmp7.UG
NPROC=1

REF=/home/m.miculan/reference/zea_mays/Zea_mays.AGPv4.dna.toplevel.fa
ANN=/home/m.miculan/reference/zea_mays/annotation/Zea_mays.AGPv4.37.chr.gtf
tmp_dir=/home/m.miculan/zea_mays_projects/eQTL/tmp_dir/${tmp}
WORKDIR=/home/m.miculan/zea_mays_projects/eQTL

VCF=/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/raw.snps.${tmp}.vcf
intervals=/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/coordinates.snp.parents.intervals

cd /home/m.miculan/zea_mays_projects/eQTL

mkdir -p ${tmp_dir}
cd ${WORKDIR}

${java} -Xmx10g -Djava.io.tmpdir=${tmp_dir} -jar ${GATK} -T UnifiedGenotyper -R ${REF} -o ${VCF} \
-I mappings/NG-6860_R_51/NG-6860_R_51.realigned_L.bam \
-I mappings/NG-6860_R_52/NG-6860_R_52.realigned_L.bam \
-I mappings/NG-7084_R_63/NG-7084_R_63.realigned_L.bam \
-I mappings/NG-7084_R_64/NG-7084_R_64.realigned_L.bam \
-I mappings/NG-7084_R_65/NG-7084_R_65.realigned_L.bam \
-I mappings/NG-7084_R_66/NG-7084_R_66.realigned_L.bam \
-I mappings/NG-7084_R_67/NG-7084_R_67.realigned_L.bam \
-I mappings/NG-7084_R_68/NG-7084_R_68.realigned_L.bam \
-I mappings/NG-7084_R_69/NG-7084_R_69.realigned_L.bam \
-I mappings/NG-7084_R_70/NG-7084_R_70.realigned_L.bam \
--genotype_likelihoods_model SNP \
--heterozygosity 0.01 \
--reference_sample_name Zea_mays.AGPv4.dna.toplevel \
-L ${intervals} \
-dcov 2000 -nt ${NPROC} >${WORKDIR}/scripts/logs/unifiedgenotyper/${tmp}.out 2>${WORKDIR}/scripts/logs/unifiedgenotyper/${tmp}.err

rm -rf ${tmp_dir}
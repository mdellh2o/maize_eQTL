#!/bin/bash

#######################################################################################################################################################################################
# SNP calling ZEA MAYS TRANSCRIPTOME
# may 2018
################################################################################################################################################################################

GATK=/home/m.miculan/sw/gatk/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
#GATK-3.3=/home/m.miculan/sw/gatk/GenomeAnalysisTK-3.3.0/GenomeAnalysisTK.jar
java=/home/m.miculan/sw/jre1.8.0_151/bin/java

tmp=multi.tmp1.UG
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
-I mappings/NG-6860_R_1/NG-6860_R_1.realigned_L.bam \
-I mappings/NG-6860_R_2/NG-6860_R_2.realigned_L.bam \
-I mappings/NG-6860_R_3/NG-6860_R_3.realigned_L.bam \
-I mappings/NG-6860_R_4/NG-6860_R_4.realigned_L.bam \
-I mappings/NG-6860_R_5/NG-6860_R_5.realigned_L.bam \
-I mappings/NG-6860_R_6/NG-6860_R_6.realigned_L.bam \
-I mappings/NG-6860_R_7/NG-6860_R_7.realigned_L.bam \
-I mappings/NG-6860_R_8/NG-6860_R_8.realigned_L.bam \
-I mappings/NG-6860_R_9/NG-6860_R_9.realigned_L.bam \
--genotype_likelihoods_model SNP \
--heterozygosity 0.01 \
--reference_sample_name Zea_mays.AGPv4.dna.toplevel \
-L ${intervals} \
-dcov 2000 -nt ${NPROC} >${WORKDIR}/scripts/logs/unifiedgenotyper/${tmp}.out 2>${WORKDIR}/scripts/logs/unifiedgenotyper/${tmp}.err

rm -rf ${tmp_dir}

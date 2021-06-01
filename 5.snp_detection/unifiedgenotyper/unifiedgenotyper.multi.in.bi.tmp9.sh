#!/bin/bash

#######################################################################################################################################################################################
# SNP calling ZEA MAYS TRANSCRIPTOME
# may 2018
################################################################################################################################################################################

GATK=/home/m.miculan/sw/gatk/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
#GATK-3.3=/home/m.miculan/sw/gatk/GenomeAnalysisTK-3.3.0/GenomeAnalysisTK.jar
java=/home/m.miculan/sw/jre1.8.0_151/bin/java

tmp=multi.in.bi.tmp9.UG
NPROC=1

REF=/home/m.miculan/reference/zea_mays/Zea_mays.AGPv4.dna.toplevel.fa
ANN=/home/m.miculan/reference/zea_mays/annotation/Zea_mays.AGPv4.37.chr.gtf
tmp_dir=/home/m.miculan/zea_mays_projects/eQTL/tmp_dir/${tmp}
WORKDIR=/home/m.miculan/zea_mays_projects/eQTL

VCF=/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/vcf_tmp/snps.${tmp}.goodQUAL.08.filtered.vcf
intervals=/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/intervals/genotypes.chr.multi.filtered.intervals

cd /home/m.miculan/zea_mays_projects/eQTL

mkdir -p ${tmp_dir}
cd ${WORKDIR}

${java} -Xmx10g -Djava.io.tmpdir=${tmp_dir} -jar ${GATK} -T UnifiedGenotyper -R ${REF} -o ${VCF} \
-I mappings/NG-6564_RIL113_lib19912_1418_5/NG-6564_RIL113_lib19912_1418_5.realigned_L.bam \
-I mappings/NG-6564_RIL114_lib19913_1418_5/NG-6564_RIL114_lib19913_1418_5.realigned_L.bam \
-I mappings/NG-6564_RIL115_lib19914_1418_5/NG-6564_RIL115_lib19914_1418_5.realigned_L.bam \
-I mappings/NG-6564_RIL116_lib19915_1418_5/NG-6564_RIL116_lib19915_1418_5.realigned_L.bam \
-I mappings/NG-6564_RIL117_lib19916_1418_5/NG-6564_RIL117_lib19916_1418_5.realigned_L.bam \
-I mappings/NG-6564_RIL118_lib19917_1418_6/NG-6564_RIL118_lib19917_1418_6.realigned_L.bam \
-I mappings/NG-6564_RIL119_lib19918_1418_6/NG-6564_RIL119_lib19918_1418_6.realigned_L.bam \
-I mappings/NG-6564_RIL120_lib19919_1418_6/NG-6564_RIL120_lib19919_1418_6.realigned_L.bam \
-I mappings/NG-6564_RIL121_lib19920_1418_6/NG-6564_RIL121_lib19920_1418_6.realigned_L.bam \
-I mappings/NG-6564_RIL122/NG-6564_RIL122.realigned_L.bam \
--genotype_likelihoods_model SNP \
--heterozygosity 0.01 \
--reference_sample_name Zea_mays.AGPv4.dna.toplevel \
-L ${intervals} \
-dcov 2000 --output_mode EMIT_ALL_SITES -nt ${NPROC} >${WORKDIR}/scripts/logs/unifiedgenotyper/${tmp}.out 2>${WORKDIR}/scripts/logs/unifiedgenotyper/${tmp}.err

rm -rf ${tmp_dir}

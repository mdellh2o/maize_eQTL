#!/bin/bash

#######################################################################################################################################################################################
# SNP calling ZEA MAYS TRANSCRIPTOME
# PARENTS
# may 2018
################################################################################################################################################################################

GATK=/home/m.miculan/sw/gatk/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
#GATK-3.3=/home/m.miculan/sw/gatk/GenomeAnalysisTK-3.3.0/GenomeAnalysisTK.jar
java=/home/m.miculan/sw/jre1.8.0_151/bin/java
REF=/home/m.miculan/reference/zea_mays/Zea_mays.AGPv4.dna.toplevel.fa
ANN=/home/m.miculan/reference/zea_mays/annotation/Zea_mays.AGPv4.37.chr.gtf
WORKDIR=/home/m.miculan/zea_mays_projects/eQTL

NPROC=1


# FIRST 4 SAMPLES - PARENTS
################################################################################################################################################################################
tmp=parents2
tmp_dir=/home/m.miculan/zea_mays_projects/eQTL/tmp_dir/${parents1}

VCF=/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/vcf_tmp/snps.UG.${tmp}.vcf

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
--genotype_likelihoods_model SNP \
--heterozygosity 0.01 \
--reference_sample_name Zea_mays.AGPv4.dna.toplevel \
-L /home/m.miculan/zea_mays_projects/eQTL/contigs.intervals \
-dcov 2000 -nt ${NPROC} >${WORKDIR}/scripts/logs/unifiedgenotyper/${tmp}.out 2>${WORKDIR}/scripts/logs/unifiedgenotyper/${tmp}.err

rm -rf ${tmp_dir}



# Second 4 SAMPLES - PARENTS
################################################################################################################################################################################
tmp=parents3
tmp_dir=/home/m.miculan/zea_mays_projects/eQTL/tmp_dir/${parents1}

VCF=/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/vcf_tmp/snps.UG.${tmp}.vcf

cd /home/m.miculan/zea_mays_projects/eQTL
mkdir -p ${tmp_dir}
cd ${WORKDIR}

${java} -Xmx5g -Djava.io.tmpdir=${tmp_dir} -jar ${GATK} -T UnifiedGenotyper -R ${REF} -o ${VCF} \
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
-L /home/m.miculan/zea_mays_projects/eQTL/contigs.intervals \
-dcov 2000 -nt ${NPROC} >${WORKDIR}/scripts/logs/unifiedgenotyper/${tmp}.out 2>${WORKDIR}/scripts/logs/unifiedgenotyper/${tmp}.err

rm -rf ${tmp_dir}


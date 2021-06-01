#!/bin/bash

#######################################################################################################################################################################################
# MAPPING
#######################################################################################################################################################################################
# a - mapping with STAR (2pass mode) 
# b - Split'N'Trim 
# c - IndelRealignment (GATK)
# d - Indels Left Aligned (GATK)


# variables
STAR=/home/m.miculan/sw/STAR-2.5.3a/bin/Linux_x86_64_static/STAR
java=/home/m.miculan/sw/jre1.8.0_151/bin/java
GATK=/home/m.miculan/sw/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar

REF=/home/m.miculan/reference/zea_mays/Zea_mays.AGPv4.dna.toplevel.fa
ANN=/home/m.miculan/reference/zea_mays/annotation/Zea_mays.AGPv4.37.chr.gtf
WORKDIR=/home/m.miculan/zea_mays_projects/eQTL

cd /home/m.miculan/zea_mays_projects/eQTL



#######################################################################################################################################################################################
# INDEXING REFERENCE
# needed only once! input: reference in fasta file, annotation in gff3
# overlapping splice junction 100bp
# version 4 of the assembly
#######################################################################################################################################################################################
# module load sw/aligners/star/2.5
cd /home/m.miculan/reference/zea_mays/STAR_index_100

#${STAR} --runThreadN 15 --runMode genomeGenerate --genomeDir /home/m.miculan/reference/zea_mays/STAR_index_100 --genomeFastaFiles ${REF} --sjdbGTFfile ${ANN} --sjdbOverhang 100 --sjdbGTFtagExonParentTranscript Parent
${STAR} --runThreadN 15 --runMode genomeGenerate --genomeDir /home/m.miculan/reference/zea_mays/STAR_index_100 --genomeFastaFiles ${REF} --sjdbGTFfile ${ANN} --sjdbOverhang 100



#######################################################################################################################################################################################
# 2a.  Alignment with STAR.
# RGDT= submission date
# before mapping with STAR 2-pass mode (suitable for GATK)
# then, run GATK SplitNCigarReads to fix MAPQ of unique reads from 255 (assigned by STAR) to 60, which is the value for GATK (255 means unknown)
# read groups already inserted in the mapping command 
#######################################################################################################################################################################################

# MAPPING ENA-ERP009123 - PARENTS
#for i in $(ls /home/m.miculan/zea_mays_projects/sequences/ENA-ERP009123_parents/trimmed/*_trimmed_1.fastq.gz); do
for SAMPLE in A632_repl1 A632_repl2 B73_repl1 B73_repl2 B73_repl3 CML91_repl1 CML91_repl2 F7_repl1 F7_repl2 H99_repl1 H99_repl2 HP301_repl1 W153R_repl1 W153R_repl2 HP301_repl2 Mo17_repl1 Mo17_repl2; do
    THR=10
    RANDOMID=$(head -5 /dev/urandom | tr -cd '[:alnum:]' | cut -c -10)
    seqdir=/home/m.miculan/zea_mays_projects/sequences/ENA-ERP009123_parents/trimmed
    INBR=${SAMPLE%_repl*}
    TMPDIR=${WORKDIR}/mappings/${SAMPLE}/tmp

    mkdir -p ${WORKDIR}/mappings/${SAMPLE} ${WORKDIR}/mappings/${SAMPLE}/logs
    cd ${WORKDIR}/mappings/${SAMPLE}
    echo "start ${SAMPLE}" >> /home/m.miculan/zea_mays_projects/eQTL/progress_mapping.txt
    
    ${STAR} --runThreadN ${THR} --genomeDir /home/m.miculan/reference/zea_mays/STAR_index_100 --sjdbGTFfile ${ANN} \
        --readFilesIn ${seqdir}/${SAMPLE}_trimmed_1.fastq.gz ${seqdir}/${SAMPLE}_trimmed_2.fastq.gz \
        --outFileNamePrefix ${WORKDIR}/mappings/${SAMPLE}/${SAMPLE}. --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --outSAMmultNmax 10 --outTmpDir ${TMPDIR} --outBAMsortingThreadN ${THR} \
        --outSAMattrRGline ID:${RANDOMID} PL:ILLUMINA CN:E-MTAB-3173 DS:Illumina_HiSeq2000 DT:2016-05-20 PU:xxx LB:${SAMPLE}_library SM:${INBR} \
        --outWigType bedGraph --outWigStrand Stranded --quantMode GeneCounts --outSAMstrandField intronMotif --readFilesCommand zcat --outWigNorm RPM --twopassMode Basic >${WORKDIR}/mappings/${SAMPLE}/logs/mapping_${SAMPLE}.out 2>${WORKDIR}/mappings/${SAMPLE}/logs/mapping_${SAMPLE}.err &&
    
    mv ${WORKDIR}/mappings/${SAMPLE}/${SAMPLE}.Aligned.sortedByCoord.out.bam ${WORKDIR}/mappings/${SAMPLE}/${SAMPLE}.Aligned.sortedByCoord.out.tmp.bam &&
    
    samtools index ${WORKDIR}/mappings/${SAMPLE}/${SAMPLE}.Aligned.sortedByCoord.out.tmp.bam &&
    
    echo "done ${SAMPLE}" >> /home/m.miculan/zea_mays_projects/eQTL/progress_mapping.txt

done


# MAPPING ENA-ERP011069
for i in $(ls /home/m.miculan/zea_mays_projects/sequences/ENA-ERP011069/trimmed/*_trimmed_1.fastq.gz); do

    THR=18
    RANDOMID=$(head -5 /dev/urandom | tr -cd '[:alnum:]' | cut -c -10)
    FULLSAMPLE=${i/_trimmed_1.fastq.gz/}; SAMPLE=$(basename $FULLSAMPLE)
    TMPDIR=${WORKDIR}/mappings/${SAMPLE}/tmp

    mkdir -p ${WORKDIR}/mappings/${SAMPLE} ${WORKDIR}/mappings/${SAMPLE}/logs
    cd ${WORKDIR}/mappings/${SAMPLE}
    echo "start ${SAMPLE}" >> /home/m.miculan/zea_mays_projects/eQTL/progress_mapping.txt
    
    ${STAR} --runThreadN ${THR} --genomeDir /home/m.miculan/reference/zea_mays/STAR_index_100 --sjdbGTFfile ${ANN} \
        --readFilesIn ${FULLSAMPLE}_trimmed_1.fastq.gz ${FULLSAMPLE}_trimmed_2.fastq.gz \
        --outFileNamePrefix ${WORKDIR}/mappings/${SAMPLE}/${SAMPLE}. --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --outSAMmultNmax 10 --outTmpDir ${TMPDIR} --outBAMsortingThreadN ${THR} \
        --outSAMattrRGline ID:${RANDOMID} PL:ILLUMINA CN:E-MTAB-3758 DS:Illumina_HiSeq2000 DT:2012-09-28 PU:xxx LB:${SAMPLE}_library SM:${SAMPLE} \
        --outWigType bedGraph --outWigStrand Stranded --quantMode GeneCounts --outSAMstrandField intronMotif --readFilesCommand zcat --outWigNorm RPM --twopassMode Basic >${WORKDIR}/mappings/${SAMPLE}/logs/mapping_${SAMPLE}.out 2>${WORKDIR}/mappings/${SAMPLE}/logs/mapping_${SAMPLE}.err &&
    
    mv ${WORKDIR}/mappings/${SAMPLE}/${SAMPLE}.Aligned.sortedByCoord.out.bam ${WORKDIR}/mappings/${SAMPLE}/${SAMPLE}.Aligned.sortedByCoord.out.tmp.bam &&
    
    samtools index ${WORKDIR}/mappings/${SAMPLE}/${SAMPLE}.Aligned.sortedByCoord.out.tmp.bam &&
    
    echo "done ${SAMPLE}" >> /home/m.miculan/zea_mays_projects/eQTL/progress_mapping.txt

done


# MAPPING ENA-ERP012784
for i in NG-6860_R_1 NG-6860_R_10 NG-6860_R_11 NG-6860_R_12 NG-6860_R_13 NG-6860_R_14 NG-6860_R_15 NG-6860_R_16 NG-6860_R_18 \
         NG-6860_R_19 NG-6860_R_2 NG-6860_R_20 NG-6860_R_22 NG-6860_R_23 NG-6860_R_24 NG-6860_R_26 NG-6860_R_27 NG-6860_R_28 \
         NG-6860_R_29 NG-6860_R_3 NG-6860_R_30 NG-6860_R_31 NG-6860_R_32 NG-6860_R_33 NG-6860_R_34 NG-6860_R_35 NG-6860_R_36 \
         NG-6860_R_37 NG-6860_R_38 NG-6860_R_39 NG-6860_R_4 NG-6860_R_40 NG-6860_R_41 NG-6860_R_42 NG-6860_R_44 NG-6860_R_45 \
         NG-6860_R_46 NG-6860_R_47 NG-6860_R_48 NG-6860_R_49 NG-6860_R_5 NG-6860_R_50 NG-6860_R_51 NG-6860_R_52 NG-6860_R_6 \
         NG-6860_R_7 NG-6860_R_8 NG-6860_R_9 NG-7084_R_53 NG-7084_R_54 NG-7084_R_55 NG-7084_R_57 NG-7084_R_58 NG-7084_R_59 NG-7084_R_60 \
         NG-7084_R_61 NG-7084_R_62 NG-7084_R_63 NG-7084_R_64 NG-7084_R_65 NG-7084_R_66 NG-7084_R_67 NG-7084_R_68 NG-7084_R_69 NG-7084_R_70 \
         NG-7084_R_72 NG-7084_R_73 NG-7084_R_74 NG-7084_R_75 NG-7084_R_76 NG-7084_R_77 NG-7084_R_78 NG-7084_R_79 NG-7084_R_80 NG-7084_R_81 \
         NG-7084_R_83 NG-7084_R_84 NG-7084_R_85 NG-7084_R_86 NG-7084_R_87 NG-7084_R_88 NG-7084_R_89 NG-7084_R_90 NG-7084_R_91 NG-7084_R_92 \
         NG-7084_R_94 NG-7084_R_71 NG-7084_R_95 NG-7084_R_96 NG-7084_R_97 NG-7084_R_98 NG-7084_R_99 NG-7084_R_93 NG-7084_R_82 
    do
    echo "start ${i}" >> ${WORKDIR}/mappings/progress_split.txt
    ${java} -Xmx7g -jar ${GATK} -T SplitNCigarReads -R ${REF} -I ${i}/${i}.Aligned.sortedByCoord.out.tmp.bam -o ${i}/${i}.Aligned.sortedByCoord.out.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS >${WORKDIR}/mappings/${i}/logs/split_${i}.out 2>${WORKDIR}/mappings/${i}/logs/split_${i}.err &&
    samtools index ${i}/${i}.Aligned.sortedByCoord.out.split.bam &&
    echo "end ${i}" >> ${WORKDIR}/mappings/progress_split.txt
done




#######################################################################################################################################################################################
# 2b-c. FIXING .bam MAPPINGS for the SNP detection
#    GATK BEST PRACTICES
#######################################################################################################################################################################################
# no duplicates removal since we are not doing DE and we have many samples, the effect should be reduced.
# read this: https://www.nature.com/articles/srep25533#discussion
# no need to run fixmate
GATK=/home/m.miculan/sw/gatk/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
java=/home/m.miculan/sw/jre1.8.0_151/bin/java

REF=/home/m.miculan/reference/zea_mays/Zea_mays.AGPv4.dna.toplevel.fa
WORKDIR=/home/m.miculan/zea_mays_projects/eQTL


# 2b.  Split'N'Trim and reassign mapping qualities (GATK)
##########################################################
for i in $(ls ${WORKDIR}/mappings/ ); do
    echo "start ${i}" >> ${WORKDIR}/mappings/progress_split.txt
    ${java} -Xmx7g -jar ${GATK} -T SplitNCigarReads -R ${REF} -I ${i}/${i}.Aligned.sortedByCoord.out.tmp.bam -o ${i}/${i}.Aligned.sortedByCoord.out.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS >${WORKDIR}/mappings/${i}/logs/split_${i}.out 2>${WORKDIR}/mappings/${i}/logs/split_${i}.err &&
    samtools index ${i}/${i}.Aligned.sortedByCoord.out.split.bam &&
    echo "end ${i}" >> ${WORKDIR}/mappings/progress_split.txt
done

# remove tmp.bam
for i in $(ls ${WORKDIR}/mappings/ ); do rm mappings/${i}/${i}.Aligned.sortedByCoord.out.tmp.ba*; done


# 2c.  IndelRealignment (GATK)
##########################################################
THR=10

# creation of file .intervals
for i in $(ls ${WORKDIR}/mappings/ ); do
    
    echo "start ${i}" >> ${WORKDIR}/progress_realignment.txt
    tmp_dir=${WORKDIR}/tmp_dir/${i}
    mkdir -p ${tmp_dir} ${WORKDIR}/logs/indel_realignment
    
    ${java} -Xmx2g -Djava.io.tmpdir=${tmp_dir} -jar ${GATK} -T RealignerTargetCreator -R ${REF} -I ${WORKDIR}/mappings/${i}/${i}.Aligned.sortedByCoord.out.split.bam -o ${WORKDIR}/mappings/${i}/${i}.intervals --num_threads ${THR} >${WORKDIR}/logs/indel_realignment/${i}.indel.creator.out 2>${WORKDIR}/logs/indel_realignment/${i}.indel.creator.err &&

    echo "end ${i}" >> ${WORKDIR}/progress_realignment.txt
    rm -rf ${tmp_dir}
done

# Indel Realigner
for i in $(ls ${WORKDIR}/mappings/ ); do
    
    echo "start ${i}" >> ${WORKDIR}/progress_realignment.txt
    tmp_dir=${WORKDIR}/tmp_dir/${i}
    mkdir -p ${tmp_dir} ${WORKDIR}/logs/indel_realignment
    
    ${java} -Xmx2g -Djava.io.tmpdir=${tmp_dir} -jar ${GATK} -T IndelRealigner -R ${REF} -I ${WORKDIR}/mappings/${i}/${i}.Aligned.sortedByCoord.out.split.bam -targetIntervals ${WORKDIR}/mappings/${i}/${i}.intervals -o ${WORKDIR}/mappings/${i}/${i}.realigned.bam >${WORKDIR}/logs/indel_realignment/${i}.indel.real.out 2>${WORKDIR}/logs/indel_realignment/${i}.indel.real.err &&

    echo "end ${i}" >> ${WORKDIR}/progress_realignment.txt &&
    rm -rf ${tmp_dir}
done


# remove old bam
for i in $(ls ${WORKDIR}/mappings/ ); do rm ${WORKDIR}/mappings/${i}/${i}.Aligned.sortedByCoord.out.*.ba*; done


# 2d.  Indels Left Aligned (GATK)
##########################################################
GATK=/home/m.miculan/sw/gatk/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
java=/home/m.miculan/sw/jre1.8.0_151/bin/java
REF=/home/m.miculan/reference/zea_mays/Zea_mays.AGPv4.dna.toplevel.fa
WORKDIR=/home/m.miculan/zea_mays_projects/eQTL

for i in $(ls ${WORKDIR}/mappings/ ); do
    
    echo "start ${i}" >> ${WORKDIR}/progress_left.txt
    tmp_dir=${WORKDIR}/tmp_dir/${i}
    mkdir -p ${tmp_dir} ${WORKDIR}/logs/left_indels
    
    ${java} -Xmx2g -Djava.io.tmpdir=${tmp_dir} -jar ${GATK} -T LeftAlignIndels -R ${REF} -I ${WORKDIR}/mappings/${i}/${i}.realigned.bam -o ${WORKDIR}/mappings/${i}/${i}.realigned_L.bam >${WORKDIR}/logs/left_indels/${i}.left.out 2>${WORKDIR}/logs/left_indels/${i}.left.err &&

    echo "end ${i}" >> ${WORKDIR}/progress_left.txt &&
    rm -rf ${tmp_dir}
done

# remove old bam
for i in $(ls ${WORKDIR}/mappings/ ); do rm ${WORKDIR}/mappings/${i}/${i}.realigned.ba* ; done

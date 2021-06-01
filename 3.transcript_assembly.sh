#!/bin/bash

#######################################################################################################################################################################################
# 3. Creation of specific transcript assembly (StringTie)
#######################################################################################################################################################################################
# variables
REF=/home/m.miculan/reference/zea_mays/Zea_mays.AGPv4.dna.toplevel.fa
ANN=/home/m.miculan/reference/zea_mays/annotation/Zea_mays.AGPv4.37.chr.gtf
WORKDIR=/home/m.miculan/zea_mays_projects/eQTL
ANN=/home/m.miculan/reference/zea_mays/annotation/Zea_mays.AGPv4.37.chr.gtf


cd ${WORKDIR}/mappings


for i in NG-6280_RIL138_P1_lib15140_1178_5 NG-6280_RIL138_P2_lib15141_1178_5 NG-6280_RIL83_P1_lib15144_1178_5 NG-6280_RIL83_P2_lib15145_1178_5 NG-6280_RIL94_P1_lib15138_1178_5 NG-6280_RIL94_P2_lib15139_1178_5 NG-6280_RIL96_P1_lib15142_1178_5 NG-6280_RIL96_P2_lib15143_1178_5 NG-6564_RIL100_lib19904_1418_4 NG-6564_RIL103_lib19905_1418_5 NG-6564_RIL104_lib19906_1418_5 NG-6564_RIL105_lib19907_1418_5 NG-6564_RIL106_lib19908_1418_5 NG-6564_RIL107_lib19909_1418_5 NG-6564_RIL109_lib19910_1418_5 NG-6564_RIL111_lib19911_1418_5 NG-6564_RIL113_lib19912_1418_5 NG-6564_RIL114_lib19913_1418_5 NG-6564_RIL115_lib19914_1418_5 NG-6564_RIL116_lib19915_1418_5 NG-6564_RIL117_lib19916_1418_5 NG-6564_RIL118_lib19917_1418_6 NG-6564_RIL119_lib19918_1418_6 NG-6564_RIL11_lib19839_1479_4 NG-6564_RIL120_lib19919_1418_6 NG-6564_RIL121_lib19920_1418_6 NG-6564_RIL122 NG-6564_RIL123_lib19922_1418_6 NG-6564_RIL124_lib19923_1418_6 NG-6564_RIL125_lib19924_1418_6 NG-6564_RIL126_lib19925_1418_6 NG-6564_RIL127_lib19926_1418_6 NG-6564_RIL128_lib19927_1418_6 NG-6564_RIL130_lib19928_1418_6 NG-6564_RIL131_lib19929_1418_7 NG-6564_RIL132_lib19930_1418_7 NG-6564_RIL133_lib19931_1418_7 NG-6564_RIL134_lib19932_1418_7 NG-6564_RIL135_lib19933_1418_7 NG-6564_RIL136_lib19934_1418_7 NG-6564_RIL137_lib19935_1418_7 NG-6564_RIL139_lib19936_1418_7 NG-6564_RIL13_lib19840_1479_4 NG-6564_RIL140_lib19937_1418_7 NG-6564_RIL141_lib19938_1418_7 NG-6564_RIL16_lib19841_1479_4 NG-6564_RIL17_lib19842_1479_4 NG-6564_RIL18_lib19843_1479_4 NG-6564_RIL19_lib19844_1479_4 NG-6564_RIL23_lib19845_1479_4 NG-6564_RIL25_lib19846_1479_2 NG-6564_RIL26_lib19847_1479_4 NG-6564_RIL27_lib19848_1479_2 NG-6564_RIL29_lib19849_1479_2 NG-6564_RIL30_lib19850_1479_2 NG-6564_RIL32_lib19851_1479_2 NG-6564_RIL33_lib19852_1479_3 NG-6564_RIL36_lib19853_1479_3 NG-6564_RIL37_lib19854_1479_2 NG-6564_RIL38_lib19855_1479_2 NG-6564_RIL39_lib19856_1479_2 NG-6564_RIL40_lib19857_1479_3 NG-6564_RIL42_lib19858_1479_2 NG-6564_RIL43_lib19859_1479_2 NG-6564_RIL44_lib19860_1479_3 NG-6564_RIL45_lib19861_1479_2 NG-6564_RIL46_lib19862_1479_2 NG-6564_RIL47_lib19863_1479_3 NG-6564_RIL48_lib19864_1479_3 NG-6564_RIL49_lib19865_1479_3 NG-6564_RIL50_lib19866_1479_3 NG-6564_RIL51_lib19867_1479_6 NG-6564_RIL52_lib19868_1479_1 NG-6564_RIL53_lib19869_1479_1 NG-6564_RIL54_lib19870_1479_1 NG-6564_RIL55_lib19871_1479_1 NG-6564_RIL57_lib19872_1479_1 NG-6564_RIL58_lib19873_1479_1 NG-6564_RIL59_lib19874_1479_1 NG-6564_RIL60_lib19875_1479_1 NG-6564_RIL61_lib19876_1479_1 NG-6564_RIL62_lib19877_1479_1 NG-6564_RIL63_lib19878_1479_1 NG-6564_RIL64_lib19879_1479_1 NG-6564_RIL65_lib19880_1479_6 NG-6564_RIL66_lib19881_1418_3 NG-6564_RIL67_lib19882_1418_3 NG-6564_RIL69_lib19883_1418_3 NG-6564_RIL70_lib19884_1418_3 NG-6564_RIL71_lib19885_1418_3 NG-6564_RIL72_lib19886_1418_3 NG-6564_RIL73_lib19887_1418_3 NG-6564_RIL74_lib19888_1418_3 NG-6564_RIL75 NG-6564_RIL76_lib19890_1418_3 NG-6564_RIL80_lib19891_1418_3 NG-6564_RIL81_lib19892_1418_3 NG-6564_RIL82_lib19893_1418_4 NG-6564_RIL84_lib19894_1418_4 NG-6564_RIL85_lib19895_1418_4 NG-6564_RIL88_lib19896_1418_4 NG-6564_RIL89_lib19897_1418_4 NG-6564_RIL90_lib19898_1418_4 NG-6564_RIL91_lib19899_1418_4 NG-6564_RIL92_lib19900_1418_4 NG-6564_RIL93_lib19901_1418_4 NG-6564_RIL95_lib19902_1418_4 NG-6564_RIL99_lib19903_1418_4 \
	NG-6860_R_1 NG-6860_R_10 NG-6860_R_11 NG-6860_R_12 NG-6860_R_13 NG-6860_R_14 NG-6860_R_15 NG-6860_R_16 NG-6860_R_18 NG-6860_R_19 NG-7084_R_71 NG-7084_R_95 \
	NG-6860_R_2 NG-6860_R_20 NG-6860_R_22 NG-6860_R_23 NG-6860_R_24 NG-6860_R_26 NG-6860_R_27 NG-6860_R_28 NG-6860_R_29 NG-6860_R_3 NG-6860_R_30 NG-6860_R_31 NG-6860_R_32 NG-6860_R_33 NG-6860_R_34 \
	NG-6860_R_35 NG-6860_R_36 NG-6860_R_37 NG-6860_R_38 NG-6860_R_39 NG-6860_R_4 NG-6860_R_40 NG-6860_R_41 NG-6860_R_42 NG-6860_R_44 NG-6860_R_45 NG-6860_R_46 NG-6860_R_47 NG-6860_R_48 NG-6860_R_49 NG-6860_R_5 NG-6860_R_50 NG-6860_R_51 NG-6860_R_52 NG-6860_R_6 NG-6860_R_7 NG-6860_R_8 NG-6860_R_9 NG-7084_R_53 NG-7084_R_54 NG-7084_R_55 NG-7084_R_57 NG-7084_R_58 NG-7084_R_59 NG-7084_R_60 NG-7084_R_61 NG-7084_R_62 NG-7084_R_63 NG-7084_R_64 NG-7084_R_65 NG-7084_R_66 NG-7084_R_67 NG-7084_R_68 NG-7084_R_69 NG-7084_R_70 NG-7084_R_71 NG-7084_R_72 NG-7084_R_73 NG-7084_R_74 NG-7084_R_75 NG-7084_R_76 NG-7084_R_77 NG-7084_R_78 NG-7084_R_79 NG-7084_R_80 NG-7084_R_81 NG-7084_R_82 NG-7084_R_83 NG-7084_R_84 NG-7084_R_85 NG-7084_R_86 NG-7084_R_87 NG-7084_R_88 NG-7084_R_89 NG-7084_R_90 NG-7084_R_91 NG-7084_R_92 NG-7084_R_93 NG-7084_R_94 NG-7084_R_95 NG-7084_R_96 NG-7084_R_97 NG-7084_R_98 NG-7084_R_99
	do
    THR=8
    BAM=${WORKDIR}/mappings/${i}/${i}.Aligned.sortedByCoord.out.split.bam
    OUTDIR=${WORKDIR}/stringtie_output/${i}
    
    mkdir -p ${OUTDIR}
    cd ${WORKDIR}/stringtie_output
    stringtie ${BAM} -o ${OUTDIR}/${i}.StringTie.out.gtf -p ${THR} -G ${ANN} -l ${i} -A ${OUTDIR}/${i}.StringTie.gene_abundance.tab -b ${OUTDIR}/${i}.ballgown.tables -e -v  >${OUTDIR}/${i}.StringTie.out 2>${OUTDIR}/${i}.StringTie.err
done




# put parents together!
#######################################################################################################################################################################################
# since we have replicates, we merge together raw counts, then we normalize
cd /home/m.miculan/zea_mays_projects/eQTL/scripts
Rscript merge_gene_abundance.r

# merge bam replicates in one single bam
cd /home/m.miculan/zea_mays_projects/eQTL/mappings
for i in A632 B73 CML91 F7 H99 HP301 Mo17 W153R; do
    mkdir -p ${i}
    samtools merge ${i}/${i}.realigned_L.bam ${i}_repl1/${i}_repl1.realigned_L.bam ${i}_repl2/${i}_repl2.realigned_L.bam; samtools index ${i}/${i}.realigned_L.bam
done

# STRINGTIE on PARENTS merged and not-merged
for i in A632_repl1 B73_repl1 CML91_repl1 F7_repl1 H99_repl1 HP301_repl1 Mo17_repl1 W153R_repl1 A632_repl2 B73_repl2 CML91_repl2 F7_repl2 H99_repl2 HP301_repl2 Mo17_repl2 W153R_repl2; do
    THR=16
    BAM=${WORKDIR}/mappings/${i}/${i}.realigned_L.bam
    OUTDIR=${WORKDIR}/stringtie_output/${i}; mkdir -p ${OUTDIR}; cd ${WORKDIR}/stringtie_output
    stringtie ${BAM} -o ${OUTDIR}/${i}.StringTie.out.gtf -p ${THR} -G ${ANN} -l ${i} -A ${OUTDIR}/${i}.StringTie.gene_abundance.tab -b ${OUTDIR}/${i}.ballgown.tables -e -v  >${OUTDIR}/${i}.StringTie.out 2>${OUTDIR}/${i}.StringTie.err
done

for i in A632 B73 CML91 F7 H99 HP301 Mo17 W153R ; do
    THR=6
    BAM=${WORKDIR}/mappings/${i}/${i}.realigned_L.bam
    OUTDIR=${WORKDIR}/stringtie_output/${i}; mkdir -p ${OUTDIR}; cd ${WORKDIR}/stringtie_output
    stringtie ${BAM} -o ${OUTDIR}/${i}.StringTie.out.gtf -p ${THR} -G ${ANN} -l ${i} -A ${OUTDIR}/${i}.StringTie.gene_abundance.tab -b ${OUTDIR}/${i}.ballgown.tables -e -v  >${OUTDIR}/${i}.StringTie.out 2>${OUTDIR}/${i}.StringTie.err
done



cd /home/m.miculan/zea_mays_projects/eQTL
awk 'BEGIN{OFS="\t"} { if ($3 == "gene" ) print $10,$1,$4,$5,$7}' ${ANN} | sed 's/"//g' | sed 's/;//' > gene_position.txt

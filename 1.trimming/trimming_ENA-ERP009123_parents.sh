#!/bin/bash

#######################################################################################################################################################################################
#  TRIMMING of each fastq (.gz) using erne-filter.
#  downloaded fastq files from server
# 
#  more cycles to parallelise processes 
#######################################################################################################################################################################################

# exp ENA-ERP009123_parents
cd /home/m.miculan/zea_mays_projects/sequences/ENA-ERP009123_parents

cut -f 28,34 experiment.table | awk 'BEGIN{FS="\t"} {split($2, aa, "/"); print aa[8], $1}'
mv ERR712363_1.fastq.gz A632_repl1_1.fastq.gz
mv ERR712363_2.fastq.gz A632_repl1_2.fastq.gz
mv ERR712361_1.fastq.gz A632_repl2_1.fastq.gz
mv ERR712361_2.fastq.gz A632_repl2_2.fastq.gz
mv ERR712351_1.fastq.gz B73_repl1_1.fastq.gz
mv ERR712351_2.fastq.gz B73_repl1_2.fastq.gz
mv ERR712354_1.fastq.gz B73_repl2_1.fastq.gz
mv ERR712354_2.fastq.gz B73_repl2_2.fastq.gz
mv ERR712365_1.fastq.gz B73_repl3_1.fastq.gz
mv ERR712365_2.fastq.gz B73_repl3_2.fastq.gz
mv ERR712360_1.fastq.gz CML91_repl1_1.fastq.gz
mv ERR712360_2.fastq.gz CML91_repl1_2.fastq.gz
mv ERR712364_1.fastq.gz CML91_repl2_1.fastq.gz
mv ERR712364_2.fastq.gz CML91_repl2_2.fastq.gz
mv ERR712355_1.fastq.gz F7_repl1_1.fastq.gz
mv ERR712355_2.fastq.gz F7_repl1_2.fastq.gz
mv ERR712358_1.fastq.gz F7_repl2_1.fastq.gz
mv ERR712358_2.fastq.gz F7_repl2_2.fastq.gz
mv ERR712362_1.fastq.gz H99_repl1_1.fastq.gz
mv ERR712362_2.fastq.gz H99_repl1_2.fastq.gz
mv ERR712357_1.fastq.gz H99_repl2_1.fastq.gz
mv ERR712357_2.fastq.gz H99_repl2_2.fastq.gz
mv ERR712356_1.fastq.gz H99_repl3_1.fastq.gz
mv ERR712356_2.fastq.gz H99_repl3_2.fastq.gz
mv ERR712353_1.fastq.gz HP301_repl1_1.fastq.gz
mv ERR712353_2.fastq.gz HP301_repl1_2.fastq.gz
mv ERR712350_1.fastq.gz HP301_repl2_1.fastq.gz
mv ERR712350_2.fastq.gz HP301_repl2_2.fastq.gz
mv ERR712359_1.fastq.gz Mo17_repl1_1.fastq.gz
mv ERR712359_2.fastq.gz Mo17_repl1_2.fastq.gz
mv ERR712366_1.fastq.gz Mo17_repl2_1.fastq.gz
mv ERR712366_2.fastq.gz Mo17_repl2_2.fastq.gz
mv ERR712367_1.fastq.gz W153R_repl1_1.fastq.gz
mv ERR712367_2.fastq.gz W153R_repl1_2.fastq.gz
mv ERR712352_1.fastq.gz W153R_repl2_1.fastq.gz
mv ERR712352_2.fastq.gz W153R_repl2_2.fastq.gz

for FULL in $(ls /home/m.miculan/zea_mays_projects/sequences/ENA-ERP009123_parents/*repl1*_1.fastq.gz); do
SAMPLE=$(basename $FULL); SS=${SAMPLE/_1.fastq.gz/}; THR=10; 
OUTDIR=/home/m.miculan/zea_mays_projects/sequences/ENA-ERP009123_parents/trimmed
erne-filter --query1 ${SS}_1.fastq.gz --query2 ${SS}_2.fastq.gz --output ${OUTDIR}/${SS}_trimmed --threads ${THR} --min-size 50 --gzip >logs/log_${SS}.out 2>logs/log_${SS}.err
done

for FULL in $(ls /home/m.miculan/zea_mays_projects/sequences/ENA-ERP009123_parents/*repl2*_1.fastq.gz); do
SAMPLE=$(basename $FULL); SS=${SAMPLE/_1.fastq.gz/}; THR=8;
OUTDIR=/home/m.miculan/zea_mays_projects/sequences/ENA-ERP009123_parents/trimmed
erne-filter --query1 ${SS}_1.fastq.gz --query2 ${SS}_2.fastq.gz --output ${OUTDIR}/${SS}_trimmed --threads ${THR} --min-size 50 --gzip >logs/log_${SS}.out 2>logs/log_${SS}.err
done

for FULL in $(ls /home/m.miculan/zea_mays_projects/sequences/ENA-ERP009123_parents/*repl3*_1.fastq.gz); do
SAMPLE=$(basename $FULL); SS=${SAMPLE/_1.fastq.gz/}; THR=8; 
OUTDIR=/home/m.miculan/zea_mays_projects/sequences/ENA-ERP009123_parents/trimmed
erne-filter --query1 ${SS}_1.fastq.gz --query2 ${SS}_2.fastq.gz --output ${OUTDIR}/${SS}_trimmed --threads ${THR} --min-size 50 --gzip >logs/log_${SS}.out 2>logs/log_${SS}.err
done

for FULL in H99_repl2_1.fastq.gz HP301_repl1_1.fastq.gz; do
SAMPLE=$(basename $FULL); SS=${SAMPLE/_1.fastq.gz/}; THR=6; 
OUTDIR=/home/m.miculan/zea_mays_projects/sequences/ENA-ERP009123_parents/trimmed
erne-filter --query1 ${SS}_1.fastq.gz --query2 ${SS}_2.fastq.gz --output ${OUTDIR}/${SS}_trimmed --threads ${THR} --min-size 50 --gzip >logs/log_${SS}.out 2>logs/log_${SS}.err
done

for FULL in H99_repl3_1.fastq.gz HP301_repl2_1.fastq.gz Mo17_repl1_1.fastq.gz Mo17_repl2_1.fastq.gz; do
SAMPLE=$(basename $FULL); SS=${SAMPLE/_1.fastq.gz/}; THR=6; 
OUTDIR=/home/m.miculan/zea_mays_projects/sequences/ENA-ERP009123_parents/trimmed
erne-filter --query1 ${SS}_1.fastq.gz --query2 ${SS}_2.fastq.gz --output ${OUTDIR}/${SS}_trimmed --threads ${THR} --min-size 50 --gzip >logs/log_${SS}.out 2>logs/log_${SS}.err
done


rm /home/m.miculan/zea_mays_projects/sequences/ENA-ERP009123_parents/*fastq.gz
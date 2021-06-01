#!/bin/bash

#######################################################################################################################################################################################
#  TRIMMING of each fastq (.gz) using erne-filter.
#  downloaded fastq files from server
# 
#  more cycles to parallelise processes 
#######################################################################################################################################################################################

# exp ENA-ERP012784
cd /home/m.miculan/zea_mays_projects/sequences/ENA-ERP012784

for FULL in $(ls /home/m.miculan/zea_mays_projects/sequences/ENA-ERP012784/NG-*1703_1_1.fastq.bz2); do
SAMPLE=$(basename $FULL); SS=${SAMPLE/_1.fastq.bz2/}; #IND=$(cut -d'_' -f5 <<<${SS})
THR=6; OUTDIR=/home/m.miculan/zea_mays_projects/sequences/ENA-ERP012784/trimmed
bunzip2 ${SS}*.fastq.bz2 &&
erne-filter --query1 ${SS}_1.fastq --query2 ${SS}_2.fastq --output ${OUTDIR}/${SS}_trimmed --threads ${THR} --min-size 50 --gzip >logs/log_${SS}.out 2>logs/log_${SS}.err &&
gzip --fast ${SS}*.fastq
done

for FULL in $(ls /home/m.miculan/zea_mays_projects/sequences/ENA-ERP012784/NG-*1703_2_1.fastq.bz2); do
SAMPLE=$(basename $FULL); SS=${SAMPLE/_1.fastq.bz2/}; #IND=$(cut -d'_' -f5 <<<${SS})
THR=6; OUTDIR=/home/m.miculan/zea_mays_projects/sequences/ENA-ERP012784/trimmed
bunzip2 ${SS}*.fastq.bz2 &&
erne-filter --query1 ${SS}_1.fastq --query2 ${SS}_2.fastq --output ${OUTDIR}/${SS}_trimmed --threads ${THR} --min-size 50 --gzip >logs/log_${SS}.out 2>logs/log_${SS}.err &&
gzip --fast ${SS}*.fastq
done

for FULL in $(ls /home/m.miculan/zea_mays_projects/sequences/ENA-ERP012784/NG-*1703_6_1.fastq.bz2); do
SAMPLE=$(basename $FULL); SS=${SAMPLE/_1.fastq.bz2/}; #IND=$(cut -d'_' -f5 <<<${SS})
THR=6; OUTDIR=/home/m.miculan/zea_mays_projects/sequences/ENA-ERP012784/trimmed
bunzip2 ${SS}*.fastq.bz2 &&
erne-filter --query1 ${SS}_1.fastq --query2 ${SS}_2.fastq --output ${OUTDIR}/${SS}_trimmed --threads ${THR} --min-size 50 --gzip >logs/log_${SS}.out 2>logs/log_${SS}.err &&
gzip --fast ${SS}*.fastq; 
done


for FULL in $(ls /home/m.miculan/zea_mays_projects/sequences/ENA-ERP012784/NG-*1746_6_1.fastq.bz2); do
SAMPLE=$(basename $FULL); SS=${SAMPLE/_1.fastq.bz2/}; #IND=$(cut -d'_' -f5 <<<${SS})
THR=6; OUTDIR=/home/m.miculan/zea_mays_projects/sequences/ENA-ERP012784/trimmed
bunzip2 ${SS}*.fastq.bz2 &&
erne-filter --query1 ${SS}_1.fastq --query2 ${SS}_2.fastq --output ${OUTDIR}/${SS}_trimmed --threads ${THR} --min-size 50 --gzip >logs/log_${SS}.out 2>logs/log_${SS}.err &&
gzip --fast ${SS}*.fastq
done

for FULL in $(ls /home/m.miculan/zea_mays_projects/sequences/ENA-ERP012784/NG-*1990_2_1.fastq.bz2); do
SAMPLE=$(basename $FULL); SS=${SAMPLE/_1.fastq.bz2/}; #IND=$(cut -d'_' -f5 <<<${SS})
THR=6; OUTDIR=/home/m.miculan/zea_mays_projects/sequences/ENA-ERP012784/trimmed
bunzip2 ${SS}*.fastq.bz2 &&
erne-filter --query1 ${SS}_1.fastq --query2 ${SS}_2.fastq --output ${OUTDIR}/${SS}_trimmed --threads ${THR} --min-size 50 --gzip >logs/log_${SS}.out 2>logs/log_${SS}.err &&
gzip --fast ${SS}*.fastq
done

for FULL in $(ls /home/m.miculan/zea_mays_projects/sequences/ENA-ERP012784/NG-*1990_7_1.fastq.bz2); do
SAMPLE=$(basename $FULL); SS=${SAMPLE/_1.fastq.bz2/}; #IND=$(cut -d'_' -f5 <<<${SS})
THR=6; OUTDIR=/home/m.miculan/zea_mays_projects/sequences/ENA-ERP012784/trimmed
bunzip2 ${SS}*.fastq.bz2 &&
erne-filter --query1 ${SS}_1.fastq --query2 ${SS}_2.fastq --output ${OUTDIR}/${SS}_trimmed --threads ${THR} --min-size 50 --gzip >logs/log_${SS}.out 2>logs/log_${SS}.err &&
gzip --fast ${SS}*.fastq
done

for FULL in $(ls /home/m.miculan/zea_mays_projects/sequences/ENA-ERP012784/NG-*1992_1_1.fastq.bz2); do
SAMPLE=$(basename $FULL); SS=${SAMPLE/_1.fastq.bz2/}; #IND=$(cut -d'_' -f5 <<<${SS})
THR=6; OUTDIR=/home/m.miculan/zea_mays_projects/sequences/ENA-ERP012784/trimmed
bunzip2 ${SS}*.fastq.bz2 &&
erne-filter --query1 ${SS}_1.fastq --query2 ${SS}_2.fastq --output ${OUTDIR}/${SS}_trimmed --threads ${THR} --min-size 50 --gzip >logs/log_${SS}.out 2>logs/log_${SS}.err &&
gzip --fast ${SS}*.fastq
done

for FULL in $(ls /home/m.miculan/zea_mays_projects/sequences/ENA-ERP012784/NG-*1992_7_1.fastq.bz2); do
SAMPLE=$(basename $FULL); SS=${SAMPLE/_1.fastq.bz2/}; #IND=$(cut -d'_' -f5 <<<${SS})
THR=6; OUTDIR=/home/m.miculan/zea_mays_projects/sequences/ENA-ERP012784/trimmed
bunzip2 ${SS}*.fastq.bz2 &&
erne-filter --query1 ${SS}_1.fastq --query2 ${SS}_2.fastq --output ${OUTDIR}/${SS}_trimmed --threads ${THR} --min-size 50 --gzip >logs/log_${SS}.out 2>logs/log_${SS}.err &&
gzip --fast ${SS}*.fastq
done

for FULL in NG-6860_R_8W__52_lib28139_1731_8 NG-7084_R_8W__70_lib34117_1990_6 NG-7084_R_8W__54_lib38678_2201_2 NG-7084_R_8W__54_lib38678_2201_4; do
SAMPLE=$(basename $FULL); SS=${SAMPLE/_1.fastq.bz2/}; #IND=$(cut -d'_' -f5 <<<${SS})
THR=6; OUTDIR=/home/m.miculan/zea_mays_projects/sequences/ENA-ERP012784/trimmed
bunzip2 ${SS}*.fastq.bz2 &&
erne-filter --query1 ${SS}_1.fastq --query2 ${SS}_2.fastq --output ${OUTDIR}/${SS}_trimmed --threads ${THR} --min-size 50 --gzip >logs/log_${SS}.out 2>logs/log_${SS}.err &&
gzip --fast ${SS}*.fastq
done

rm /home/m.miculan/zea_mays_projects/sequences/ENA-ERP012784/*fastq.gz


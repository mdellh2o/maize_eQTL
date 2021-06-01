#!/bin/bash


#######################################################################################################################################################################################
# TRIMMING of each fastq (.gz) using erne-filter
# downloaded fastq files from server
#
# more cycles to parallelise processes 
#######################################################################################################################################################################################

# exp ENA-ERP011069
cd /home/m.miculan/zea_mays_projects/sequences/ENA-ERP011069

for FULL in $(ls NG-*_1_1.fastq.bz2); do
SAMPLE=$(basename $FULL); SS=${SAMPLE/_1.fastq.bz2/}; #IND=$(cut -d'_' -f2 <<<${SS})
THR=6; OUTDIR=/home/m.miculan/zea_mays_projects/sequences/ENA-ERP011069/trimmed
bunzip2 ${SS}*.fastq.bz2 &&
erne-filter --query1 ${SS}_1.fastq --query2 ${SS}_2.fastq --output ${OUTDIR}/${SS}_trimmed --threads ${THR} --min-size 50 --gzip >logs/trim_${SS}.out 2>logs/trim_${SS}.err &&
gzip ${SS}*.fastq
done

for FULL in $(ls /home/m.miculan/zea_mays_projects/sequences/ENA-ERP011069/NG-*_2_1.fastq.bz2); do
SAMPLE=$(basename $FULL); SS=${SAMPLE/_1.fastq.bz2/}; #IND=$(cut -d'_' -f2 <<<${SS})
THR=6; OUTDIR=/home/m.miculan/zea_mays_projects/sequences/ENA-ERP011069/trimmed
bunzip2 ${SS}*.fastq.bz2 &&
erne-filter --query1 ${SS}_1.fastq --query2 ${SS}_2.fastq --output ${OUTDIR}/${SS}_trimmed --threads ${THR} --min-size 50 --gzip >logs/trim_${SS}.out 2>logs/trim_${SS}.err &&
gzip ${SS}*.fastq
done

for FULL in $(ls /home/m.miculan/zea_mays_projects/sequences/ENA-ERP011069/NG-*_3_1.fastq.bz2); do
SAMPLE=$(basename $FULL); SS=${SAMPLE/_1.fastq.bz2/}; #IND=$(cut -d'_' -f2 <<<${SS})
THR=6; OUTDIR=/home/m.miculan/zea_mays_projects/sequences/ENA-ERP011069/trimmed
bunzip2 ${SS}*.fastq.bz2 &&
erne-filter --query1 ${SS}_1.fastq --query2 ${SS}_2.fastq --output ${OUTDIR}/${SS}_trimmed --threads ${THR} --min-size 50 --gzip >logs/trim_${SS}.out 2>logs/trim_${SS}.err &&
gzip ${SS}*.fastq
done

for SS in NG-6564_RIL113_lib19912_1418_5 NG-6564_RIL122 NG-6564_RIL75; do
THR=6; OUTDIR=/home/m.miculan/zea_mays_projects/sequences/ENA-ERP011069/trimmed
erne-filter --query1 ${SS}_1.fastq --query2 ${SS}_2.fastq --output ${OUTDIR}/${SS}_trimmed --threads ${THR} --min-size 50 --gzip >logs/trim_${SS}.out 2>logs/trim_${SS}.err &&
gzip ${SS}*.fastq
done

for FULL in $(ls /home/m.miculan/zea_mays_projects/sequences/ENA-ERP011069/NG-*_4_1.fastq.bz2); do
SAMPLE=$(basename $FULL); SS=${SAMPLE/_1.fastq.bz2/}; #IND=$(cut -d'_' -f2 <<<${SS})
THR=6; OUTDIR=/home/m.miculan/zea_mays_projects/sequences/ENA-ERP011069/trimmed
bunzip2 ${SS}*.fastq.bz2 &&
erne-filter --query1 ${SS}_1.fastq --query2 ${SS}_2.fastq --output ${OUTDIR}/${SS}_trimmed --threads ${THR} --min-size 50 --gzip >logs/trim_${SS}.out 2>logs/trim_${SS}.err &&
gzip --fast ${SS}*.fastq 
done

for FULL in $(ls /home/m.miculan/zea_mays_projects/sequences/ENA-ERP011069/NG-*_5_1.fastq.bz2); do
SAMPLE=$(basename $FULL); SS=${SAMPLE/_1.fastq.bz2/}; #IND=$(cut -d'_' -f2 <<<${SS})
THR=6; OUTDIR=/home/m.miculan/zea_mays_projects/sequences/ENA-ERP011069/trimmed
bunzip2 ${SS}*.fastq.bz2 &&
erne-filter --query1 ${SS}_1.fastq --query2 ${SS}_2.fastq --output ${OUTDIR}/${SS}_trimmed --threads ${THR} --min-size 50 --gzip >logs/trim_${SS}.out 2>logs/trim_${SS}.err &&
gzip --fast ${SS}*.fastq
done

for FULL in $(ls /home/m.miculan/zea_mays_projects/sequences/ENA-ERP011069/NG-*_6_1.fastq.bz2); do
SAMPLE=$(basename $FULL); SS=${SAMPLE/_1.fastq.bz2/}; #IND=$(cut -d'_' -f2 <<<${SS})
THR=6; OUTDIR=/home/m.miculan/zea_mays_projects/sequences/ENA-ERP011069/trimmed
bunzip2 ${SS}*.fastq.bz2 &&
erne-filter --query1 ${SS}_1.fastq --query2 ${SS}_2.fastq --output ${OUTDIR}/${SS}_trimmed --threads ${THR} --min-size 50 --gzip >logs/trim_${SS}.out 2>logs/trim_${SS}.err &&
gzip --fast ${SS}*.fastq
done

for FULL in $(ls /home/m.miculan/zea_mays_projects/sequences/ENA-ERP011069/NG-*_7_1.fastq.bz2); do
SAMPLE=$(basename $FULL); SS=${SAMPLE/_1.fastq.bz2/}; #IND=$(cut -d'_' -f2 <<<${SS})
THR=6; OUTDIR=/home/m.miculan/zea_mays_projects/sequences/ENA-ERP011069/trimmed
bunzip2 ${SS}*.fastq.bz2 &&
erne-filter --query1 ${SS}_1.fastq --query2 ${SS}_2.fastq --output ${OUTDIR}/${SS}_trimmed --threads ${THR} --min-size 50 --gzip >logs/trim_${SS}.out 2>logs/trim_${SS}.err &&
gzip --fast ${SS}*.fastq
done

rm /home/m.miculan/zea_mays_projects/sequences/ENA-ERP011069/*fastq.gz


#######################################################################################################################################################################################
# this is the README of ANALYSIS eQTL ZEA MAYS TRANSCRIPTOME
#######################################################################################################################################################################################
# steps:
# INDEXING reference .fasta with STAR (overlapping splice junctions 100bp) 
# 1. TRIMMING with erne-filter (done once!)
# 2. mapping with STAR software and fixing Bam 
#     a. mapping
#     b. Split'N'Trim and reassign mapping qualities (GATK)
#     c. IndelRealignment (GATK)
#     d. Left Alignment (GATK)
# 3. Creation of specific transcript assembly (StringTie) (contemporary step 3)
# 4. Extract transcript abundance values and filter samples
# 5. SNP detection from RNAseq data
# 6. Genome Wide Association Analysis
# 7. transcript-trait correlation analysis
# 8. 


# software and reference versions
java --> jre1.8.0_151
STAR --> STAR-2.5.3a
GATK --> GenomeAnalysisTK-3.8-0-ge9d806836

REFERENCE --> Zea_mays.AGPv4.dna.toplevel.fa
ANNOTATION --> Zea_mays.AGPv4.37.chr.gtf


#######################################################################################################################################################################################
# 1.  TRIMMING of each fastq (.gz) using erne-filter.
#     downloaded fastq files from server
#######################################################################################################################################################################################
scripts in --> 1.trimming folder

# exp ENA-ERP011069
trimming_ENA-ERP011069.sh

# exp ENA-ERP012784
trimming_ENA-ERP012784.sh

# exp ENA-ERP009123_parents
trimming_ENA-ERP009123_parents.sh



#######################################################################################################################################################################################
# 2. Alignment with STAR.
# RGDT= submission date
# before mapping with STAR 2-pass mode (suitable for GATK)
# then, run GATK SplitNCigarReads to fix MAPQ of unique reads from 255 (assigned by STAR) to 60, which is the value for GATK (255 means unknown)
# read groups already inserted in the mapping command 
#######################################################################################################################################################################################
2.mapping.sh



#######################################################################################################################################################################################
# 3.  Creation of specific transcript assembly via StringTie. In output table there are normalized expression values as FPKM and/or TPM
#######################################################################################################################################################################################
3.transcript_assembly.sh



#######################################################################################################################################################################################
# 4.  Extract FPKM/TPM from StringTie. Write a table with genes name and unstranded raw counts.
#######################################################################################################################################################################################
4a.extract_transcript_abundance.r

# output:
# /home/m.miculan/zea_mays_projects/eQTL/stringtie_output/all.tpm.chr.txt



# filter Expression
###############################################################################################################################################################
# FILTER genes that are expressed with a variance at least 0.5 and expressed in 60% of samples
4b.filter_expression_data.r

# output: "/home/m.miculan/zea_mays_projects/eQTL/stringtie_output/all.tpm.chr.0.5var.txt"



#######################################################################################################################################################################################
# 5.  SNPs detection
# for samples with 2 libraries we considered only P1
#######################################################################################################################################################################################
all scripts in the folder --> 5.snp_detection
start with 5.workflow_snps_UG.sh



###############################################################################################################################################
# 6 - Genome Wide Association Analysis
###############################################################################################################################################
6.GWAS.anlaysis.R



###############################################################################################################################################
# 7 - Correlation analysis
###############################################################################################################################################
7.RNA.pheno.cor.splitted.pops.r

# gene ontology on genes coming from trait associations with transcript levels  analysis 



###############################################################################################################################################
# 8.  eQTL analysis
###############################################################################################################################################
#scripts in the folder --> 6.eqtl_analysis

# expression_file_all --> all.tpm.chr.0.05var.txt
# SNP_file_name_all --> genotypes.chr.all.filtered.biallelic.txt
# covariates_file_name --> PCA_covariates_top10_and_pop.all.txt


# CALCULATE PCA on expression data and on genotypes - ALLPOPs
#########################################################################
8a.calculate_pca.allpop.r #script run with different number of PCAs, selected the best with 10


# CALCULATE eQTL - ALLPOPs 
#########################################################################
8b.eQTL_calculation_all.r 


# select independent SNPs (not in the same LD block)
##########################################
8c.select_snps_by_LD.all.r

# results:
# gwa.eqtl.red.txt, pcc.eqtl.red
# OBJECTS: me, eqtl_maf, pcc.eqtl, gwa.eqtl, int.eqtl, pcc.gwa.eqtl --> me.all.pca10.all.reduced.2.Rdata



# 9.GRAPHS and TABLEs
##########################################
folder 9.graphs_and_tables

# scripts:
plot.GWA.MVP.r
fig2.graph.correlation.r

prepare_data.r 
get.tables.from.eQTL
fig4.scatterplots.eQTL.r
figS5.histogram.eQTL.r

fig5.boxplot.loop.r



######################################################################################################
######################################################################################################
# 10. NETWORK ANALYSIS 
# pops splitted & analysis considering all genes! work to find the consensus between both populations
######################################################################################################
######################################################################################################
# scripts in the folder --> 10.wgcna

# clean data and create consensus network
10a.wgcna_onestep.consensus.r

# relate phenotypic data with consensus modules
10b.wgcna_relating_traits.consensus.r

# prepare data to visualization the network
10c.wgcna_annotation_visualization.consensus.r

# plots 
10d.wgcna_plots_consensus.r
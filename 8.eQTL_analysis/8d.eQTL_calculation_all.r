###############################################################################################################################################
# ANALYSIS eQTL
###############################################################################################################################################
# STEPs:
# 1. calculated PCA o enpression data and on genotypes to use as covariates
# 2. subsampled sNPs and gene expression with p-value 1
# 3. eQTL with smaller p-value

#install.packages("MatrixEQTL")
#install.packages("GenABEL")

library(MatrixEQTL)
library(data.table)
# library(GenABEL)
library(stringr) # to manipulate strings


# VARIABLES 
###############################################################################################################################################
workdir="/home/maize.vib/eqtl_results_2"
base.dir = find.package("MatrixEQTL") # to find the location of the package

useModel = modelLINEAR; # modelANOVA or modelLINEAR or modelLINEAR_CROSS 

SNP_file_name = "/home/maize.vib/genotypic_data_2/genotypes.chr.all.filtered.biallelic.txt"
snps_location_file_name = "/home/maize.vib/genotypic_data_2/snp.pos.filtered.biallelic.txt"

expression_file_name = "/home/maize.vib/expression_data/all.tpm.chr.0.05var.txt"
gene_location_file_name = "/home/maize.vib/annotation/gene_positions_V4.txt"

covariates_file_name = "/home/maize.vib/eqtl_covariates/PCA_covariates_top10_all.txt"

errorCovariance = numeric()
pvOutputThreshold = 1e-7 #to set the threshold, p-value threshold determines which gene-SNP associations are saved in the output file output_file_name

# errorCovariance = read.table("Sample_Data/errorCovariance.txt");


setwd(workdir)


# Variables RUN1 
####################################################
# object: me.all.pca10

# LOADING files with genotipes, gene expression, and covariates
# SlicedData is a class for storing large matrices (for fast and memory efficient manipulations with large datasets presented in matrix form)
############################################################################################################################################################

snps = SlicedData$new() 
snps$fileDelimiter = "\t"      # the TAB character
snps$fileOmitCharacters = "NA" # denote missing values;
snps$fileSkipRows = 1          # one row of column labels
snps$fileSkipColumns = 1       # one column of row labels
snps$fileSliceSize = 2000      # read file in pieces of 2,000 rows
snps$LoadFile( SNP_file_name )


gene = SlicedData$new()
gene$fileDelimiter = "\t"      # the TAB character
gene$fileOmitCharacters = "NA" # denote missing values;
gene$fileSkipRows = 1          # one row of column labels
gene$fileSkipColumns = 1       # one column of row labels
gene$fileSliceSize = 2000      # read file in pieces of 2,000 rows
gene$LoadFile( expression_file_name )



# CHECK OUTLIERS IN GENOTYPE
#######################################################
maf.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
  slice = snps[[sl]];
  maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
  maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
}
maf = unlist(maf.list)

## Look at the distribution of MAF
pdf("/home/maize.vib/eqtl_results_2/eqtl_graphs/hist_maf.pdf")
hist(maf[maf<0.1],seq(0,0.1,length.out=100))
dev.off()

cat('SNPs before filtering:',nrow(snps))
# snps$RowReorderSimple(maf>0.05);
snps$RowReorder(maf>0.05);
cat('SNPs before filtering:',nrow(snps))



#  COVARIATES FILE - PCA TOP10
############################################################################################################################################################
output_file_name = paste0(workdir, "/", "eqtl_all_pca10.txt")
#covariates_file_name = PCA_file_top10

cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"      # the TAB character
cvrt$fileOmitCharacters = "NA" # denote missing values;
cvrt$fileSkipRows = 1          # one row of column labels
cvrt$fileSkipColumns = 1       # one column of row labels
cvrt$fileSliceSize = 2000      # read file in pieces of 2,000 rows
if(length(covariates_file_name)>0) { cvrt$LoadFile(covariates_file_name) }



#  RUN1 - MAIN function - PCA TOP10
############################################################################################################################################################
# me.all.pca10 = Matrix_eQTL_engine(
#   snps = snps,
#   gene = gene,
#   cvrt = cvrt,
#   output_file_name = output_file_name,
#   pvOutputThreshold = pvOutputThreshold,
#   useModel = useModel,
#   errorCovariance = errorCovariance,
#   verbose = TRUE,
#   pvalue.hist = "qqplot",
#   min.pv.by.genesnp = FALSE,
#   noFDRsaveMemory = FALSE)


# # t-statistic is equal to the effect size estimate divided by its standard error. 
# # Thus, the standard deviation of the estimate can be easily calculated as the ration of the effect size estimates and their t-statistics
# me = me.all.pca10
# me.all.pca10$all$eqtls$beta_se = me$all$eqtls$beta / me$all$eqtls$statistic

# # correlations from t-statistics:
# tstat = me$all$eqtls$statistic
# r = tstat / sqrt( dfFull + tstat^2 )
# R2 = r^2
# me.all.pca10$all$eqtls$r2 = R2


# save(me.all.pca10, SNP_file_name, expression_file_name, covariates_file_name, file=paste0(workdir, "/", "me.all.pca10.Rdata"))



#  RUN2 - MAIN function - PCA TOP10 - MAIN
############################################################################################################################################################
output_file_name = paste0(workdir, "/", "eqtl_all_pca10_mainf.2.txt")
pvOutputThreshold = 1e-7 #to set the threshold, p-value threshold determines which gene-SNP associations are saved in the output file output_file_name

me.all.pca10.2 = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = TRUE,
  noFDRsaveMemory = FALSE);

# me.all.pca10.2 = me.all.pca10


# t-statistic is equal to the effect size estimate divided by its standard error. 
# Thus, the standard deviation of the estimate can be easily calculated as the ration of the effect size estimates and their t-statistics
me = me.all.pca10.2
me$all$eqtls$beta_se = me$all$eqtls$beta / me$all$eqtls$statistic


# correlations from t-statistics:
################################
dfFull = me$param$dfFull # degrees of freedom
tstat = me$all$eqtls$statistic
r = tstat / sqrt( dfFull + tstat^2 )
R2 = r^2
me$all$eqtls$r2 = R2


# not shown in paper
pdf('pvalues_distribution.all.pdf')
hist(me$all$eqtls$"pvalue", nclass = 50, main = "histogram pvalues all")
dev.off()

pdf('fdr_distribution.all.pdf')
hist(me$all$eqtls$"FDR", nclass = 50, main = "histogram FDR all")
dev.off()

pdf('fdr-r2.all.pdf')
plot(me$all$eqtls$"FDR", me$all$eqtls$"r2", main = "FDR vs R2 all")
dev.off()

pdf('min.pv.snps.all.pdf')
hist(me$all$min.pv.snps[me$all$min.pv.snps <= 1e-15] , nclass = 30, main = "histogram min pvalues snps all")
dev.off()

pdf('min.pv.gene.all.pdf')
hist(me$all$min.pv.snps[me$all$min.pv.snps <= 1e-15], nclass = 30, main = "histogram min pvalues gene all")
dev.off()


# FILTER eQTL
################################
me.df <- as.data.frame(me$all$eqtls) #not ordered!

library(plyr)
library(LSD)
library(stringr)
library(dplyr)

# eqtl <- me.df %>% filter(me.df$"FDR" < 0.05)
maf <- unlist(lapply(seq_len(length(snps)), function(sl) {
  slice <- snps[[sl]]
  af <- rowMeans(slice, na.rm = TRUE) / 2
  pmin(af, 1 - af)
}))
maf_df <- data.frame(snp = snps$GetAllRowNames(), maf = maf, stringsAsFactors = FALSE) %>% tbl_df
dim(maf_df)


save(me, maf_df, file=paste0(workdir, "/", "me.all.pca10.all.Rdata"))




# pvalue threshold:
pval.threshold <- 0.01/nrow(geno)

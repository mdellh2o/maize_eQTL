###############################################################################################################################################
# ANALYSIS eQTL of SUB-SAMPLED DATA
###############################################################################################################################################
#install.packages("MatrixEQTL")
#install.packages("GenABEL")

library(MatrixEQTL)
library(data.table)
library(GenABEL)


###############################################################################################################################################
# 1. variables (these variables are to be kept! ALWAYS during this analysis till the symbol ***
###############################################################################################################################################
base.dir = find.package("MatrixEQTL") # to find the location of the package
useModel = modelLINEAR; # modelANOVA or modelLINEAR or modelLINEAR_CROSS 

SNP_file_name = "/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/genotypes.chr.all.filtered.biallelic.txt"
expression_file_name = "/home/m.miculan/zea_mays_projects/eQTL/stringtie_output/all.tpm.chr.0.5var.txt"
#population_file_name= "/home/m.miculan/zea_mays_projects/eQTL/covariates/population_covariates.txt"
errorCovariance = numeric()

subSNPFILE = "/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/genotypes_sub_1000.all.txt"
subexpfile = "/home/m.miculan/zea_mays_projects/eQTL/stringtie_output/all.tpm.sub1000.txt"


# Create SUBSAMPLE to easily generate qqplot for initial evaluation of type 1 error
###############################################################################################################################################
ngenes=1000
ngeno=1000

SNPdata<-fread(SNP_file_name,data.table=F)
expdata<-fread(expression_file_name,data.table=F)

subSNPdata<-SNPdata[sample(nrow(SNPdata),ngeno),]
subexpdata<-expdata[sample(nrow(expdata),ngenes),]

write.table(subSNPdata,subSNPFILE,row.names=F,quote=F,sep="\t")
write.table(subexpdata,subexpfile,row.names=F,quote=F,sep="\t")


 
###############################################################################################################################################
# RUN matrix_eQTL on sampled data
###############################################################################################################################################
# these files remain the same
SNP_file_name= subSNPFILE
expression_file_name= subexpfile


# LOADING files with genotipes, gene expression, and covariates
################################################################################################################################################
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

# cvrt = SlicedData$new()
# cvrt$fileDelimiter = "\t"      # the TAB character
# cvrt$fileOmitCharacters = "NA" # denote missing values;
# cvrt$fileSkipRows = 1          # one row of column labels
# cvrt$fileSkipColumns = 1       # one column of row labels
# cvrt$fileSliceSize = 2000      # read file in pieces of 2,000 rows
# if(length(covariates_file_name)>0) { cvrt$LoadFile(covariates_file_name) }


################################################################################################################################################
# *** END of fixed variables ***
################################################################################################################################################


# from here at each run we have to:
# 1. load variables names (section VARIABLES)
# 2. load the covariates (section load COVARIATES)
# 3. run the eQTL function (MAIN function)
# 4. run the analysis part (GRAPHS)


# Section: VARIABLES of 1000-SUBSAMPLE run WITHOUT COVARIATES AT ALL 
###############################################################
covariates_file_name=character()
output_file_name="/home/m.miculan/zea_mays_projects/eQTL/eQTL_results/subsampled/eQTL_sub1000.all.txt"
pvOutputThreshold = 1 


# Variables of 1000-SUBSAMPLE run with PCA - TOP1
covariates_file_name="/home/m.miculan/zea_mays_projects/eQTL/covariates/PCA_covariates_top1_all.txt"
output_file_name="/home/m.miculan/zea_mays_projects/eQTL/eQTL_results/subsampled/eQTL_sub1000.PCA-TOP1.all.txt"
pvOutputThreshold = 1 


# Variables of 1000-SUBSAMPLE run with PCA - TOP3
covariates_file_name="/home/m.miculan/zea_mays_projects/eQTL/covariates/PCA_covariates_top3_all.txt"
output_file_name="/home/m.miculan/zea_mays_projects/eQTL/eQTL_results/subsampled/eQTL_sub1000.PCA-TOP3.all.txt"
pvOutputThreshold = 1 


# Variables of 1000-SUBSAMPLE run with PCA - TOP5
covariates_file_name="/home/m.miculan/zea_mays_projects/eQTL/covariates/PCA_covariates_top5_all.txt"
output_file_name="/home/m.miculan/zea_mays_projects/eQTL/eQTL_results/subsampled/eQTL_sub1000.PCA-TOP5.all.txt"
pvOutputThreshold = 1 


# Variables of 1000-SUBSAMPLE run with PCA - TOP10
covariates_file_name="/home/m.miculan/zea_mays_projects/eQTL/covariates/PCA_covariates_top10_all.txt"
output_file_name="/home/m.miculan/zea_mays_projects/eQTL/eQTL_results/subsampled/eQTL_sub1000.PCA-TOP10.all.txt"
pvOutputThreshold = 1 




#  Section: Load COVARIATES
############################################################################################################################################################
cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"      # the TAB character
cvrt$fileOmitCharacters = "NA" # denote missing values;
cvrt$fileSkipRows = 1          # one row of column labels
cvrt$fileSkipColumns = 1       # one column of row labels
cvrt$fileSliceSize = 2000      # read file in pieces of 2,000 rows
if(length(covariates_file_name)>0) { cvrt$LoadFile(covariates_file_name) }


#  Section: MAIN function
############################################################################################################################################################
me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)



# Section: GRAPHS - ANALYSIS - FDR
######################################################################################################################################################

# qqplot subsample WITHOUT covariates at all
eQTL_file_name = "/home/m.miculan/zea_mays_projects/eQTL/eQTL_results/subsampled/eQTL_sub1000.all.txt"
pdf_file_name= "/home/m.miculan/zea_mays_projects/eQTL/graphs/qqplot_subsampled/all.qqplot_without_cov.pdf"
pp<-fread(eQTL_file_name,data.table=F)
datafdr<-sort(pp$FDR)
theofdr<-ppoints(datafdr)
datafdr<--log10(datafdr)
theofdr<--log10(theofdr)
#Cosa vietata che si fa solo ora per vedere!!!
#datafdr[datafdr>50]<-50
#Fine della cosa vietata!!
sdatafdr<-datafdr[seq(1,length(datafdr),by=10)]
stheofdr<-theofdr[seq(1,length(theofdr),by=10)]
pdf(pdf_file_name)
plot(stheofdr,sdatafdr,pch=19,cex=0.5)
abline(0,1,col="red")
dev.off()


# qqplot subsample with PCAtop1 as covariates
eQTL_file_name="/home/m.miculan/zea_mays_projects/eQTL/eQTL_results/subsampled/eQTL_sub1000.PCA-TOP1.all.txt"
pdf_file_name= "/home/m.miculan/zea_mays_projects/eQTL/graphs/qqplot_subsampled/all.qqplot.PCA-TOP1.pdf"
pp<-fread(eQTL_file_name,data.table=F)
datafdr<-sort(pp$FDR)
theofdr<-ppoints(datafdr)
datafdr<--log10(datafdr)
theofdr<--log10(theofdr)
#Cosa vietata che si fa solo ora per vedere!!!
#datafdr[datafdr>50]<-50
#Fine della cosa vietata!!
sdatafdr<-datafdr[seq(1,length(datafdr),by=10)]
stheofdr<-theofdr[seq(1,length(theofdr),by=10)]
pdf(pdf_file_name)
plot(stheofdr,sdatafdr,pch=19,cex=0.5)
abline(0,1,col="red")
dev.off()


# qqplot subsample with PCAtop3 as covariates
eQTL_file_name="/home/m.miculan/zea_mays_projects/eQTL/eQTL_results/subsampled/eQTL_sub1000.PCA-TOP3.all.txt"
pdf_file_name= "/home/m.miculan/zea_mays_projects/eQTL/graphs/qqplot_subsampled/all.qqplot.PCA-TOP3.pdf"
pp<-fread(eQTL_file_name,data.table=F)
datafdr<-sort(pp$FDR)
theofdr<-ppoints(datafdr)
datafdr<--log10(datafdr)
theofdr<--log10(theofdr)
#Cosa vietata che si fa solo ora per vedere!!!
#datafdr[datafdr>50]<-50
#Fine della cosa vietata!!
sdatafdr<-datafdr[seq(1,length(datafdr),by=10)]
stheofdr<-theofdr[seq(1,length(theofdr),by=10)]
pdf(pdf_file_name)
plot(stheofdr,sdatafdr,pch=19,cex=0.5)
abline(0,1,col="red")
dev.off()


# qqplot subsample with PCAtop5 as covariates
eQTL_file_name="/home/m.miculan/zea_mays_projects/eQTL/eQTL_results/subsampled/eQTL_sub1000.PCA-TOP5.all.txt"
pdf_file_name= "/home/m.miculan/zea_mays_projects/eQTL/graphs/qqplot_subsampled/all.qqplot.PCA-TOP5.pdf"
pp<-fread(eQTL_file_name,data.table=F)
datafdr<-sort(pp$FDR)
theofdr<-ppoints(datafdr)
datafdr<--log10(datafdr)
theofdr<--log10(theofdr)
#Cosa vietata che si fa solo ora per vedere!!!
#datafdr[datafdr>50]<-50
#Fine della cosa vietata!!
sdatafdr<-datafdr[seq(1,length(datafdr),by=10)]
stheofdr<-theofdr[seq(1,length(theofdr),by=10)]
pdf(pdf_file_name)
plot(stheofdr,sdatafdr,pch=19,cex=0.5)
abline(0,1,col="red")
dev.off()



# qqplot subsample with PCAtop10 as covariates
eQTL_file_name="/home/m.miculan/zea_mays_projects/eQTL/eQTL_results/subsampled/eQTL_sub1000.PCA-TOP10.all.txt"
pdf_file_name= "/home/m.miculan/zea_mays_projects/eQTL/graphs/qqplot_subsampled/all.qqplot.PCA-TOP10.pdf"
pp<-fread(eQTL_file_name,data.table=F)
datafdr<-sort(pp$FDR)
theofdr<-ppoints(datafdr)
datafdr<--log10(datafdr)
theofdr<--log10(theofdr)
#Cosa vietata che si fa solo ora per vedere!!!
#datafdr[datafdr>50]<-50
#Fine della cosa vietata!!
sdatafdr<-datafdr[seq(1,length(datafdr),by=10)]
stheofdr<-theofdr[seq(1,length(theofdr),by=10)]
pdf(pdf_file_name)
plot(stheofdr,sdatafdr,pch=19,cex=0.5)
abline(0,1,col="red")
dev.off()



# me is the eQTL without any covariate
summary(me$all$eqtls$pvalue)





# Section: GRAPHS - ANALYSIS - p-value
######################################################################################################################################################

# qqplot subsample WITHOUT covariates at all
eQTL_file_name = "/home/m.miculan/zea_mays_projects/eQTL/eQTL_results/subsampled/eQTL_sub1000.all.txt"
pdf_file_name= "/home/m.miculan/zea_mays_projects/eQTL/graphs/qqplot_subsampled/all.qqplot_without_cov_pvalue.pdf"
pp<-fread(eQTL_file_name,data.table=F)
datapv<-sort(pp$"p-value")
theopv<-ppoints(datapv)
datapv<--log10(datapv)
theopv<--log10(theopv)
sdatapv<-datapv[seq(1,length(datapv),by=10)]
stheopv<-theopv[seq(1,length(theopv),by=10)]
pdf(pdf_file_name)
plot(stheopv,sdatapv,pch=19,cex=0.5)
abline(0,1,col="red")
dev.off()


# qqplot subsample with PCAtop1 as covariates
eQTL_file_name="/home/m.miculan/zea_mays_projects/eQTL/eQTL_results/subsampled/eQTL_sub1000.PCA-TOP1.all.txt"
pdf_file_name= "/home/m.miculan/zea_mays_projects/eQTL/graphs/qqplot_subsampled/all.qqplot.PCA-TOP1_pvalue.pdf"
pp<-fread(eQTL_file_name,data.table=F)
datapv<-sort(pp$"p-value")
theopv<-ppoints(datapv)
datapv<--log10(datapv)
theopv<--log10(theopv)
sdatapv<-datapv[seq(1,length(datapv),by=10)]
stheopv<-theopv[seq(1,length(theopv),by=10)]
pdf(pdf_file_name)
plot(stheopv,sdatapv,pch=19,cex=0.5)
abline(0,1,col="red")
dev.off()


# qqplot subsample with PCAtop3 as covariates
eQTL_file_name="/home/m.miculan/zea_mays_projects/eQTL/eQTL_results/subsampled/eQTL_sub1000.PCA-TOP3.all.txt"
pdf_file_name= "/home/m.miculan/zea_mays_projects/eQTL/graphs/qqplot_subsampled/all.qqplot.PCA-TOP3_pvalue.pdf"
pp<-fread(eQTL_file_name,data.table=F)
datapv<-sort(pp$"p-value")
theopv<-ppoints(datapv)
datapv<--log10(datapv)
theopv<--log10(theopv)
sdatapv<-datapv[seq(1,length(datapv),by=10)]
stheopv<-theopv[seq(1,length(theopv),by=10)]
pdf(pdf_file_name)
plot(stheopv,sdatapv,pch=19,cex=0.5)
abline(0,1,col="red")
dev.off()


# qqplot subsample with PCAtop5 as covariates
eQTL_file_name="/home/m.miculan/zea_mays_projects/eQTL/eQTL_results/subsampled/eQTL_sub1000.PCA-TOP5.all.txt"
pdf_file_name= "/home/m.miculan/zea_mays_projects/eQTL/graphs/qqplot_subsampled/all.qqplot.PCA-TOP5_pvalue.pdf"
pp<-fread(eQTL_file_name,data.table=F)
datapv<-sort(pp$"p-value")
theopv<-ppoints(datapv)
datapv<--log10(datapv)
theopv<--log10(theopv)
sdatapv<-datapv[seq(1,length(datapv),by=10)]
stheopv<-theopv[seq(1,length(theopv),by=10)]
pdf(pdf_file_name)
plot(stheopv,sdatapv,pch=19,cex=0.5)
abline(0,1,col="red")
dev.off()



# qqplot subsample with PCAtop10 as covariates
eQTL_file_name="/home/m.miculan/zea_mays_projects/eQTL/eQTL_results/subsampled/eQTL_sub1000.PCA-TOP10.all.txt"
pdf_file_name= "/home/m.miculan/zea_mays_projects/eQTL/graphs/qqplot_subsampled/all.qqplot.PCA-TOP10_pvalue.pdf"
pp<-fread(eQTL_file_name,data.table=F)
datapv<-sort(pp$"p-value")
theopv<-ppoints(datapv)
datapv<--log10(datapv)
theopv<--log10(theopv)
sdatapv<-datapv[seq(1,length(datapv),by=10)]
stheopv<-theopv[seq(1,length(theopv),by=10)]
pdf(pdf_file_name)
plot(stheopv,sdatapv,pch=19,cex=0.5)
abline(0,1,col="red")
dev.off()



# me is the eQTL without any covariate
summary(me$all$eqtls$pvalue)




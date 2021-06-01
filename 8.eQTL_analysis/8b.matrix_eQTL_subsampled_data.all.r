###############################################################################################################################################
# ANALYSIS eQTL of SUB-SAMPLED DATA to SELECT the NUMBER of COVARIATES
###############################################################################################################################################
library(MatrixEQTL)
library(data.table)
library(dplyr)
#library(GenABEL)


###############################################################################################################################################
# 1. variables (these variables are to be kept! ALWAYS during this analysis till the symbol ***
###############################################################################################################################################
base.dir = find.package("MatrixEQTL") # to find the location of the package
useModel = modelLINEAR; # modelANOVA or modelLINEAR or modelLINEAR_CROSS 

workdir= "/glab/projects/eQTL_maize"
SNP_file_name = "/glab/projects/eQTL_maize/input_data/genotypes.chr.all.filtered.biallelic.txt"
expression_file_name = "/glab/projects/eQTL_maize/input_data/all.tpm.chr.0.05var.txt"

subSNPFILE = "/glab/projects/eQTL_maize/input_data/genotypes_sub_1000.all.txt"
subexpfile = "/glab/projects/eQTL_maize/input_data/all.tpm.sub1000.txt"

outdir="/glab/projects/eQTL_maize/comments_reviewers/covariates/eQTL_subsample_PCA"

#load("/glab/projects/eQTL_maize/comments_reviewers/covariates/pc_analysis.RData") #covariate data


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


################################################################################################################################################
# *** END of fixed variables ***
################################################################################################################################################
# from here at each run we have to:
# 1. load variables names (section VARIABLES)
# 2. load the covariates (section load COVARIATES)
# 3. run the eQTL function (MAIN function)
# 4. run the analysis part (GRAPHS)



# Section: VARIABLES of 1000-SUBSAMPLE run WITHOUT COVARIATES AT ALL 
################################################################################################################################################
covariates_file_name=character()
errorCovariance = numeric()

output_file_name= paste0(outdir, "/eQTL_sub1000_PCA-TOP0_all.txt")
pvOutputThreshold = 1 

#  Section: Load COVARIATES
##############################################################################
cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"      # the TAB character
cvrt$fileOmitCharacters = "NA" # denote missing values;
cvrt$fileSkipRows = 1          # one row of column labels
cvrt$fileSkipColumns = 1       # one column of row labels
cvrt$fileSliceSize = 2000      # read file in pieces of 2,000 rows
if(length(covariates_file_name)>0) { cvrt$LoadFile(covariates_file_name) }

#  Section: MAIN function
###############################################################################
me = Matrix_eQTL_main(
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



# Section: VARIABLES of 1000-SUBSAMPLE run WITH DIFFERENT NUMBER OF PCs (0 to 20)
################################################################################################################################################
pvOutputThreshold = 1 
errorCovariance = numeric()


for (i in 1:20) {
    covariates_file_name = paste0("/glab/projects/eQTL_maize/comments_reviewers/covariates/cov_files/PCA_top", i, "_all.txt")
    output_file_name = paste0("/glab/projects/eQTL_maize/comments_reviewers/covariates/eQTL_subsample_PCA/eQTL_sub1000_PCA-TOP", i, "_all.txt")

    #  Section: Load COVARIATES
    ###########################################
    cvrt = SlicedData$new()
    cvrt$fileDelimiter = "\t"      # the TAB character
    cvrt$fileOmitCharacters = "NA" # denote missing values;
    cvrt$fileSkipRows = 1          # one row of column labels
    cvrt$fileSkipColumns = 1       # one column of row labels
    cvrt$fileSliceSize = 2000      # read file in pieces of 2,000 rows
    if(length(covariates_file_name)>0) { cvrt$LoadFile(covariates_file_name) }


    #  Section: MAIN function
    ###########################################
    # use engine instead off main to produce the graph
    me = Matrix_eQTL_engine(
        snps = snps,
        gene = gene,
        cvrt = cvrt,
        output_file_name = output_file_name,
        pvOutputThreshold = pvOutputThreshold,
        useModel = useModel,
        errorCovariance = errorCovariance,
        verbose = TRUE,
        pvalue.hist = "qqplot",
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE)
    
    pdf(paste0("/glab/projects/eQTL_maize/comments_reviewers/covariates/eQTL_subsample_PCA/pval_engine_PC", i, ".pdf")
    plot(me, pch = 16, cex = 0.7)
    dev.off()

} #for i 



# Section: GRAPHS - Q-Q plot 
######################################################################################################################################################
# http://www.sthda.com/english/wiki/qq-plots-quantile-quantile-plots-r-base-graphs
library("car")
library(data.table)

eQTLdir = "/glab/projects/eQTL_maize/comments_reviewers/covariates/eQTL_subsample_PCA"


# #FDR
# ###########################################################################
# for (i in 0:20) {
#   eQTL_file_name = paste0(eQTLdir, "/eQTL_sub1000_PCA-TOP", i, "_all.txt")
#   pdf_file_name = paste0(eQTLdir, "/qqplot_TOP", i, "_FDR.pdf")
#   eqtl<-fread(eQTL_file_name, data.table=F)

#   datafdr<-sort(eqtl$FDR)
#   theofdr<-ppoints(datafdr)
#   datafdr<--log10(datafdr)
#   theofdr<--log10(theofdr)
#   sdatafdr<-datafdr[seq(1,length(datafdr),by=10)]
#   stheofdr<-theofdr[seq(1,length(theofdr),by=10)]

  
#   pdf(pdf_file_name)
#     #qqnorm(eqtl$"FDR", pch = 1, frame=FALSE)
#     #qqline(eqtl$"FDR", col = "red", lwd = 2)
#     #qqPlot(eqtl$"FDR")
#     plot(stheofdr,sdatafdr,pch=19,cex=0.5)
#     abline(0,1,col="red")
#   dev.off()
# }


# ############################################################
# # qqplot subsample WITHOUT covariates at all
# eQTL_file_name = paste0(outdir, "/eQTL_sub1000.all_noCov.txt")
# pdf_file_name= paste0("/glab/projects/eQTL_maize/comments_reviewers/covariates/eQTL_subsample_PCA/all.qqplot_without_cov.pdf")
# eqtl<-fread(eQTL_file_name, data.table=F)
# datafdr<-sort(eqtl$FDR)
# theofdr<-ppoints(datafdr)
# datafdr<--log10(datafdr)
# theofdr<--log10(theofdr)
# #Cosa vietata che si fa solo ora per vedere!!!
# #datafdr[datafdr>50]<-50
# #Fine della cosa vietata!!
# sdatafdr<-datafdr[seq(1,length(datafdr),by=10)]
# stheofdr<-theofdr[seq(1,length(theofdr),by=10)]
# pdf(pdf_file_name)
# plot(stheofdr,sdatafdr,pch=19,cex=0.5)
# abline(0,1,col="red")
# dev.off()



# num of eQTL with FDR < 0.05
#################################
library(magicfor)
library(data.table)
library(dplyr)
eQTLdir = "/glab/projects/eQTL_maize/comments_reviewers/covariates/eQTL_subsample_PCA"

magic_for(silent = TRUE)

for (i in 0:20) {
  eQTL_file_name = paste0(eQTLdir, "/eQTL_sub1000_PCA-TOP", i, "_all.txt")
  eqtl<-fread(eQTL_file_name, data.table=F)
  fdr05 <- nrow(eqtl %>% filter(FDR < 0.05))
  fdr01 <- nrow(eqtl %>% filter(FDR < 0.01))

  put(fdr05, fdr01)
}

magic_result_as_dataframe()








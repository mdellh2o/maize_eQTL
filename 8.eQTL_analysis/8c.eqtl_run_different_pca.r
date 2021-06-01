###############################################################################################################################################
# RUN eQTL with different PCs as COVARIATES --> to choose the right number
###############################################################################################################################################
library(MatrixEQTL)
library(data.table)
library(dplyr)
library(magicfor)



###############################################################################################################################################
# 1. variables (these variables are to be kept! ALWAYS during this analysis till the symbol ***
###############################################################################################################################################
base.dir = find.package("MatrixEQTL") # to find the location of the package
useModel = modelLINEAR; # modelANOVA or modelLINEAR or modelLINEAR_CROSS 

workdir= "/glab/projects/eQTL_maize"
SNP_file_name = "/glab/projects/eQTL_maize/input_data/genotypes.chr.all.filtered.biallelic.txt"
expression_file_name = "/glab/projects/eQTL_maize/input_data/all.tpm.chr.0.05var.txt"
errorCovariance = numeric()

#outdir="/glab/projects/eQTL_maize/"
#load("/glab/projects/eQTL_maize/comments_reviewers/covariates/pc_analysis.RData") #covariate dat

 
###############################################################################################################################################
# RUN matrix_eQTL on sampled data
##############################################################################################################################################

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
# 2. load the covariates in the loop  (section load COVARIATES)
# 3. run the eQTL function (MAIN function)
# 4. run the analysis part (GRAPHS)



# Section: VARIABLES of 1000-SUBSAMPLE run WITHOUT COVARIATES AT ALL 
################################################################################################################################################
covariates_file_name=character()
output_file_name= NULL
pvOutputThreshold = 1e-7



# Section: VARIABLES of 1000-SUBSAMPLE run WITH DIFFERENT NUMBER OF PCs (0 to 20)
################################################################################################################################################
magic_for(silent = TRUE)

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

    # num of eQTL with FDR < 0.05
    fdr05 <- nrow(me$all$eqtls %>% filter(FDR < 0.05))
    fdr01 <- nrow(me$all$eqtls %>% filter(FDR < 0.01))
    pv05 <- nrow(me$all$eqtls %>% filter(pvalue < 0.05))
    pv01 <- nrow(me$all$eqtls %>% filter(pvalue < 0.01))
    put(fdr05, fdr01, pv05, pv01)

} #for i 


df <- magic_result_as_dataframe()








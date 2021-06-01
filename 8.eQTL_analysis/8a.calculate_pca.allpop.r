#  CALCULATE PCA on expression data and on genotypes (BOTH POPULATION)
#########################################################################
# conda activate eQTL
R
library(data.table)
library(factoextra)

expression_file_name = "/glab/projects/eQTL_maize/input_data/all.tpm.chr.0.05var.txt"
SNP_file_name = "/glab/projects/eQTL_maize/input_data/genotypes.chr.all.filtered.biallelic.txt"
cov_file_name = "/glab/projects/eQTL_maize/input_data/population_covariates.txt"

setwd("/glab/projects/eQTL_maize/comments_reviewers/covariates")


# covariates file of population
cov <- fread(cov_file_name, sep="\t", data.table=F, header = T)
cov <- as.data.frame(t(cov[, -1]))
names(cov) <- c("cov")
cov$"population" <- "NA"
cov[which(cov$cov == 2), ]$"population" <- "parent"
cov[which(cov$cov == 1), ]$"population" <- "multiparental"
cov[which(cov$cov == 0), ]$"population" <- "biparental"



##################################
# PRINCIPAL COMPONENTS ANALYSIS
##################################

# PCA of Expression Data
rr<-fread(expression_file_name, sep="\t", data.table=F)
pca.rr<-prcomp(t(rr[,-1]), scale=FALSE)
#pcr<-pca.rr$x

# PCA of Genotypes
gg<-fread(SNP_file_name, sep="\t", data.table=F)
gg=na.omit(gg)
pca.gg<-prcomp(t(gg[,-1]), scale=FALSE)
#pcp<-mypcag$x



# SELECTION of PCs NUMBER
########################################################################
# prepare file with n PCs
dir.create("cov_files")

for (i in 1:20) {
    top.r <- t(pca.rr$x[,c(1:i)]) # take PC-i
    rownames(top.r)<-paste0("Exp_",colnames(pca.rr$x)[1:i])

    top.g <- t(pca.gg$x[,c(1:i)]) # take PC-i
    rownames(top.g)<-paste0("Geno_",colnames(pca.gg$x)[1:i])
    
    top<-rbind(top.r,top.g)
    top<-as.data.frame(top)
    top<-cbind(row.names(top),top)
    names(top)[1]<-"covariate"

    PCA_file = paste0("cov_files/PCA_top", i, "_all.txt")
    write.table(top,PCA_file,quote=F,sep="\t",row.names=FALSE)
}



# PLOT VARIANCE EXPLAINED 
########################################################################
library(factoextra)


pdf("variance_explained_expr.pdf")
fviz_eig(pca.rr,addlabels = TRUE, ylim = c(0, 15))
dev.off()


pdf("variance_explained_geno.pdf")
fviz_eig(pca.gg, addlabels = TRUE, ylim = c(0, 15))
dev.off()


#extract some infos
eig.val.r <- get_eigenvalue(pca.rr) # to get variance % and cumulative vaiance %
var.r <- get_pca_var(pca.rr) # list of matrices containing all the results for the active variables (coordinates, correlation between variables and axes, squared cosine and contributions)

eig.val.g <- get_eigenvalue(pca.gg) # to get variance % and cumulative vaiance %
var.g <- get_pca_var(pca.gg) # list of matrices containing all the results for the active variables (coordinates, correlation between variables and axes, squared cosine and contributions)

# 20 Pcs explain 63.38809% of variance in expresison data and 44.57511% in genotype data
# we should have taken num. of PCs to explain either 70% or till eigenvalue >1 (Kaiser 1961)
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/


# graphs of individuals (DO IT REMOVING LABELS)
ind.r <-  get_pca_ind(pca.rr)
ind.g <-  get_pca_ind(pca.gg)

pdf("pca_ind_expression.pdf")
fviz_pca_ind(pca.rr, label = "none", col.ind = cov$"population",
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
            legend.title = "population"
)
dev.off()

pdf("pca_ind_geno.pdf")
fviz_pca_ind(pca.gg, label = "none", col.ind = cov$"population",
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
)
dev.off()


# save results
#########################################
save(expression_file_name, cov_file_name, SNP_file_name, pca.rr, pca.gg, eig.val.r, eig.val.g,var.g, var.r, 
    file = "pc_analysis.RData")


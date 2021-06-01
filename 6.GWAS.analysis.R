maindir<-"home/output"
setwd(maindir)

library(MVP)

# curhmp<-"../input/genotypic.data.hmp"
# curphe<-"../input/allpheno.txt"
# MVP.Data(fileHMP=curhmp,
		# filePhe=curphe,
		# sep.hmp="\t",
        # sep.phe="\t",
        # SNP.effect="Add",
        # fileKin=FALSE, 
		# filePC=FALSE, 
        # out="mvp.hmp",
        # priority="speed"
        # )

#load geno, map, and pheno file generated here above
genotype <- attach.big.matrix("mvp.hmp.geno.desc")
phenotype <- read.table("mvp.hmp.phe",head=TRUE)
phenotype<-phenotype[,-c(2:4)]
#select only relevant phenos
phenotype<-phenotype[,c(1,3,7,2,15)]
map <- read.table("mvp.hmp.map" , head = TRUE)

#assign best models
#bestmod<-read.table("best.models.MD.txt", header=T)

#get in best model directory
setwd("best.model.out")

#run MVP
for(i in 2:ncol(phenotype)){
  curname<-colnames(phenotype)[i]
  imMVP <- MVP(
    phe=phenotype[, c(1, i)],
    geno=genotype,
    map=map,
 #   K=Kinship,
 #   CV.GLM=Covariates,
 #   CV.MLM=Covariates,
 #   CV.FarmCPU=Covariates,
    nPC.FarmCPU=bestmod[grep(curname, bestmod[,1]),2],
    perc=1,
    priority="speed",
    ncpus=6,
    vc.method="EMMA",
    maxLoop=10,
    method.bin="FaST-LMM",#"FaST-LMM","EMMA", "static"
    permutation.threshold=F,
   # permutation.rep=199,
    threshold=0.1,
    method="FarmCPU"
  )
  save(imMVP, file=paste(curname, "MVP.object.Rdata", sep="."))
}


maindir<-"/home/maize.vib/phenotypes/splitted.pops"
setwd(maindir)
library(data.table)
library(string)

#load in phenotypes
pheno<-read.delim("/home/maize.vib/PUBBLICATION/7.RNA.pheno.cor/allpheno.txt", header=T)
rownames(pheno)<-pheno[,1]
pheno<-pheno[,-1]

#load in expression data already filtered by variance (5%)
expr<-fread("/home/maize.vib/PUBBLICATION/7.RNA.pheno.cor/all.tpm.chr.0.05var.txt", data.table = F)
rownames(expr)<-expr[,1]
expr<-expr[,-1]

#fix colnames in expression data
sampnames<-colnames(expr)
sampnames<-sub("^.*\\_R", "R",sampnames)
sampnames<-sub("R_", "RIL_8W_",sampnames)
sampnames<-sub("\\_lib.*$", "",sampnames)
sampnames<-sub("\\_P.*$", "",sampnames)

colnames(expr)<-sampnames

#fix dataframes for correlations
texpr_backup<-t(expr)

#divide by population 
bh<-pheno[grepl("BxH|AVERAGE", pheno$"family"),]
bh<-bh[,-c(1:3)]

magic<-pheno[grep(("MAGIC|AVERAGE"), pheno$"family"),] 
magic<-magic[,-c(1:3)]


#make it a list and go in a loop
poplist<-list(bh, magic)
names(poplist)<-c("bh", "magic")

#make an out list to store results
outlist<-list()

for(i in 1:length(poplist)){
	tmp<-poplist[[i]]
	#subset them to contain the same samples
	texpr<-texpr_backup[order(rownames(texpr_backup)),] #BUG was present
	tmp<-tmp[order(rownames(tmp)),]

	texpr<-texpr[which(rownames(texpr) %in% rownames(tmp)),]
	tmp<-tmp[which(rownames(tmp) %in% rownames(texpr)),]

	stopifnot(all(rownames(texpr) == rownames(tmp)))

	#calculate correlations
	crp<-cor(tmp, texpr, use = "pairwise.complete.obs",  method ="pearson")
	#crs<-cor(tmp, texpr, use = "pairwise.complete.obs",  method ="spearman")

	#store them in the outlist
	outlist[[i]]<-crp
	
}#for i

names(outlist)<-names(poplist)
save(outlist, file="/home/maize.vib/PUBBLICATION/7.RNA.pheno.cor/RNA.pheno.correlations.Rdata")
# load("/home/maize.vib/PUBBLICATION/7.RNA.pheno.cor/RNA.pheno.correlations.Rdata")

#################################
#calculate permutations and extract quantiles

#set n of permutation and list to store results
nperm<-7000
permoutlist<-list()

for(p in 1:length(poplist)){
	tmp<-poplist[[p]]
	tmp=tmp[ , c("LER", "leaf.length", "DZ.size", "LED")] # MM
	
	#set dataframe to store results
	qout<-matrix(nrow=6, ncol=ncol(tmp))
	colnames(qout)<-colnames(tmp)
	rownames(qout)<-c("q001","q005","q01","q99", "q995","q999")

	#subset them to contain the same samples
	texpr<-texpr_backup[order(rownames(texpr_backup)),]
	tmp<-tmp[order(rownames(tmp)),]

	texpr<-texpr[which(rownames(texpr) %in% rownames(tmp)),]
	tmp<-tmp[which(rownames(tmp) %in% rownames(texpr)),]

	stopifnot(all(rownames(texpr) == rownames(tmp)))

	#initialize permutation obj
	crpermpheno<-list()
	
	for (i in 1:ncol(tmp)){
			print(colnames(tmp)[i])

			tmpp<-data.frame(tmp[,i])
			rownames(tmpp)<-rownames(tmp)
			names(tmpp) <- names(tmp)[i] # MM 

			#initialize permutation obj
			crperm<-list()

			for (j in 1:nperm){
					if(j%%100==0){
							print(j)
					}
					tmpp[,1]<-sample(tmpp[,1], size = nrow(tmpp))
					crperm[[j]]<-cor(tmpp[,1], texpr,  use = "pairwise.complete.obs", method ="pearson")
			}#for j

			crpermpheno[[i]]<-unlist(crperm)
			names(crpermpheno)[i]<-colnames(tmp)[i]

			#extract quantiles
			qout[1,i]<-quantile(crpermpheno[[i]], 0.001, na.rm=T)
			qout[2,i]<-quantile(crpermpheno[[i]], 0.005, na.rm=T)
			qout[3,i]<-quantile(crpermpheno[[i]], 0.01, na.rm=T)
			qout[4,i]<-quantile(crpermpheno[[i]], 0.99, na.rm=T)
			qout[5,i]<-quantile(crpermpheno[[i]], 0.995, na.rm=T)
			qout[6,i]<-quantile(crpermpheno[[i]], 0.999, na.rm=T)

	}#for i

	#store in list object and save separatedly
	permoutlist[[j]]<-qout
	
	write.table(qout, file=paste0("permutation.quantiles.", names(poplist)[p],".txt"), quote=F, sep="\t")
#	write.table(qout, file=paste0("permutation.quantiles.", names(poplist)[p],".2red.txt"), quote=F, sep="\t")
#	save(qout, file=paste0("permutation.quantiles.", names(poplist)[p],".Rdata"))
	
	
	save(crpermpheno, file=paste0("permutation.", names(poplist)[p],".Rdata"))
#	save(crpermpheno, file=paste0("permutations.", names(poplist)[p],".only4.Rdata"))

}#for p
 
 

######################
#now get significant correlations for each population

#load in correlations
load("/home/maize.vib/PUBBLICATION/7.RNA.pheno.cor/RNA.pheno.correlations.Rdata")

#set list to store results
popsig<-list()

# outlist is a list with 2 elements BH and MULTI, each of them has the correlations between pheno and genes
for (i in 1:length(outlist)){
	tmp<-outlist[[i]] #load correlation
	tmp <- tmp[c("LER", "leaf.length", "DZ.size", "LED"), ]
	curpop<-names(outlist)[i]
	
	quants<-read.table(paste0("/home/maize.vib/phenotypes/splitted.pops/permutation.quantiles.", curpop,".txt"))
	
	#check that everything is as expected
	stopifnot(all(colnames(quants)==rownames(tmp)))
	
	#set list to store results
	outsig<-list()
	for (j in 1:nrow(tmp)){
		#get thresholds
		thrlow<-quants[2,j] # !!! ATTENTION, THE QUANTILE! 
		thrhi<-quants[5,j] # !!! ATTENTION, THE QUANTILE! 
		
		#set significant genes
		sigs<-which(tmp[j,]<thrlow | tmp[j,]>thrhi)
		outsig[[j]]<-colnames(tmp)[sigs]
		print(paste(curpop, "-", rownames(tmp)[j], "- genes corr:", length(sigs)))
		print(paste(curpop, "-", rownames(tmp)[j], "- negat corr:", length(which(tmp[j,sigs]<0))))
		print(paste(curpop, "-", rownames(tmp)[j], "- posit corr:", length(which(tmp[j,sigs]>0))))
	}
	names(outsig)<-rownames(tmp)
	popsig[[i]]<-outsig
}

#now get overlapping genes
overlap<-list()
for (i in 1:length(popsig[[1]])){
	overlap[[i]]<-intersect(popsig[[1]][[i]], popsig[[2]][[i]])
}#for i
names(overlap)<-names(outsig) # then it will be replaced


# how many genes overlap in each phenotype
overlap.corr <- list()
tmp <- data.frame()

for (i in 1:length(overlap)) {
  tmp <- outlist[[1]][i,(overlap)[[i]]] # create a small data frame for each pheno with the correlation values for both pop each with
  tmp <- rbind(tmp, outlist[[2]][i,(overlap)[[i]]])
  row.names(tmp) <- c("biparental", "multiparental")
  overlap.corr[[i]] <- tmp
  names(overlap.corr)[[i]] <- names(overlap)[[i]]
  
  print(paste(names(overlap)[[i]], "- tot genes", length(overlap[[i]]))) # how many overlapping genes each phenotype
  print(paste(names(overlap)[[i]], "- negat corr:", length(which(tmp[1, ]<0 & tmp[2,]<0))))
  print(paste(names(overlap)[[i]], "- posit corr:", length(which(tmp[1, ]>0 & tmp[2,]>0))))
  overlap[[i]] <- names(which( (tmp[1, ]<0 & tmp[2,]<0) | (tmp[1, ]>0 & tmp[2,]>0) ))
} #for i 



#save dataset
save(overlap, overlap.corr, file="/home/maize.vib/PUBBLICATION/7.RNA.pheno.cor/overlapping.genes.Rdata")


# CREATE TABLE FOR PUBLICATION
library(dplyr)
#traits<-c("leaf.length", "DZ.size", "LER", "LED")
pcc.table <- t(outlist[[1]])
pcc.table <- as_tibble(pcc.table[, traits])
pcc.table <- rename(pcc.table, leaf.length.bi = leaf.length)
pcc.table <- rename(pcc.table, DZ.size.bi = DZ.size)
pcc.table <- rename(pcc.table, LER.bi = LER)
pcc.table <- rename(pcc.table, LED.bi = LED)

pcc.table <- as_tibble(cbind(pcc.table, t(outlist[[2]][traits, ])))
pcc.table <- rename(pcc.table, leaf.length.multi = leaf.length)
pcc.table <- rename(pcc.table, DZ.size.multi = DZ.size)
pcc.table <- rename(pcc.table, LER.multi = LER)
pcc.table <- rename(pcc.table, LED.multi = LED)
pcc.table$"gene_id" <- colnames(outlist[[1]])

write.table(pcc.table,"/home/maize.vib/PUBBLICATION/7.RNA.pheno.cor/correlations.txt", row.names=F,quote=F,sep="\t")

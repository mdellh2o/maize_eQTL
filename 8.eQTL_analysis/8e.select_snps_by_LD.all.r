###############################################################################################################################################
# PARSE PCC by LD
###############################################################################################################################################
R

library(data.table)
library(stringr) # to manipulate strings
library(dplyr)
library (plyr)

# variables
workdir <- "/home/maize.vib/eqtl_results_2"
setwd(workdir)

load("me.all.pca10.all.Rdata")

phenotypes <- c("LER", "LED", "leaf.length", "DZ.size") # define which phenotypes

# define a function
'%!in%' <- function(x,y)!('%in%'(x,y))

#get LD info
LDdecay<-read.delim("/home/maize.vib/LD/subsetbothpop/subsampling/LD.decay.info.txt")
#move it to bp if necessary
if(max(LDdecay[3,])<1000){
  LDdecay[3,]<-LDdecay[3,]*1e6
}



# get eQTL output and prepare file
################################################ 
#work gene by gene
#eqtl <- fread(eqtl_file_name, data.table=F)

eqtl <- me$all$eqtls

eqtl$"snps" <- as.character(eqtl$"snps")
eqtl$"snp.chr" <- as.numeric(do.call(rbind, strsplit(eqtl$"snps", "_", fixed=TRUE))[,1])
eqtl$"snp.pos" <- as.numeric(do.call(rbind, strsplit(eqtl$"snps", "_", fixed=TRUE))[,2])
eqtl <- eqtl[, c("snps", "snp.chr", "snp.pos", "gene", "statistic", "pvalue", "FDR", "beta")]
gene_to_check <- unique(eqtl$"gene")
snp_to_check <- unique(eqtl$"snps")


# LOAD GWAS RESULTs
################################################ 
# add columns to eqtl table, one for each GWA.phenotype results, to record if the SNP was found even in the GWAs analysis

snps.notok <- list() #create list of sign snps from GWAS that are not present in the eqtl analysis (if they are not ;) )

for(pheno in phenotypes){
  gwas.perm.LD.reduced.file = paste0("/home/maize.vib/results_dataframes/GWAS/LD_reduced/LD.reduced.", pheno, ".GWAS.pvalue.MVP.perm.hits.Rdata") # object --> LDredGWA
  load(gwas.perm.LD.reduced.file) 
  LDredGWA.perm = LDredGWA # significant SNP_ID obtained by GWAS permutations and selected with LD decay 
  LDredGWA.perm$"SNP" = gsub("S","", LDredGWA.perm$"SNP")
  
  gwa.perm = LDredGWA.perm$"SNP" # ID of sign SNPs
  
  eqtl[ ,paste0("gwa.", pheno)] <- 0 # add column of the phenotype
  #eqtl[which(eqtl$"SNP" %in% gwa.perm), paste0("gwa.", pheno)] <- 1 # write 1 if the SNPs from GWA analysis is present in the eqtl results
  eqtl[is.element(eqtl$"snps", gwa.perm), paste0("gwa.", pheno)] <- 1 # write 1 if the SNPs from GWA analysis is present in the eqtl results
  
  # store in the object snps.notok the snps that are not found in eqtl results
  if(length(which(gwa.perm %!in% snp_to_check == TRUE)) > 0) {
    #snps.notok[[pheno]] <- which(gwa.perm %!in% snp_to_check)
    snps.notok[[pheno]] <- gwa.perm[which(gwa.perm %!in% snp_to_check)]
  } else { snps.notok[[pheno]] <- 0 } 
} # for pheno

snps.notok <- ldply(snps.notok, data.frame)
snps.notok


# SELECT SNPs indipendent in eQTL table
################################################ 
#pval_thr <- as.data.frame(quantile(eqtl$pvalue))[3, ] # to reduce computation I will take only the eqtl within the 75th quantile
pval_thr <- 0.01/373769 # total number of SNPs
eqtl <- eqtl[which(eqtl$"pvalue" <= pval_thr), ]  # look at the distribution of the min.pvalue! 
dim(eqtl)
eqtl <- eqtl[order(as.numeric(eqtl[,"snp.chr"]), as.numeric(eqtl[,"snp.pos"])), ]
bychr <- split(eqtl, eqtl$"snp.chr") # create a list where each element come from eqtl results for one chr


#set output list
outbygene <- list() #this will go inside outchr list and it will have one element each gene
outchr <- list()


# START!
for(ch in 1:length(bychr)){
  
  curchr<-bychr[[ch]] # working one chr each time
  
  table_gene <- ddply(curchr, .(gene), nrow) # return number of SNPs each gene

    for (g in 1:nrow(table_gene)) {
    
    curgene <- curchr[which(curchr$"gene" == as.character(table_gene[ g, "gene"])), ] #there was a problem, this might solve it
    curgene <- curgene[order(as.numeric(curgene[,"snp.pos"])),]
    
    #add LD peak information
    LDpeakidx<-1
    curgene$"LDpeak"<-LDpeakidx
    
    curLD<-LDdecay[ 3, ch] #get LD specific to this chromosome in Mb
    
    idx<-curchr[1,"snp.pos"] #get the position of the first snp
    
    if(nrow(curgene)>1) {
      
      #look at all the other SNP and take note whether they belong to the same LD peak
      for(j in 2:nrow(curgene)) {
        #if distance is above LD decay, assign subsequent LD peak
        if(abs(curgene[j,"snp.pos"]-idx) > 2*curLD){
          LDpeakidx<-LDpeakidx+1 #change LD peak number
          curgene$"LDpeak"[j]<-LDpeakidx #update the LDpeak number
          
          idx<-curgene[j,"snp.pos"] #update position in the index
          #if distance is below LD decay, then just update SNP position
        } else {
          #update position
          idx<-curgene[j,"snp.pos"]
          curgene$LDpeak[j]<-LDpeakidx
        } # # if abs position
      } # j nrow(curgene)
      
      #now extract from each LD peak JUST the most significant SNP
      #llply(curgene, .(LDpeak), min(curpeak[,"pvalue"]))
      bypeak<-split(curgene, curgene$"LDpeak")
      
      curgene.reduced<-list() #set list to output # 1 snp per ogni picco
      
      for (p in 1:length(bypeak)){
        curpeak<-bypeak[[p]]
        check.gwas <- as.vector(rowSums(curpeak[ , c("gwa.LER", "gwa.LED", "gwa.leaf.length", "gwa.DZ.size")]))
        if(any(check.gwas > 0 )) {
          curgene.reduced[[p]]<-curpeak[rowSums(curpeak[ , c("gwa.LER", "gwa.LED", "gwa.leaf.length", "gwa.DZ.size")]) > 0 , ]
        } else { 
          # if the SNP was found in the GWAS considered that one otherwise take the one with the lowest pvalue
          curgene.reduced[[p]]<-curpeak[which(curpeak[,"pvalue"]==min(curpeak[,"pvalue"], na.rm=T)),]
        } #if check.gwas
      } #for p peak
          
      # transfor the list curgene.reduced[[p]] into a data frame
      curgene.reduced <- do.call(rbind, curgene.reduced) # each gene
    } #if length curgene>1
    
    outbygene[[g]] <- curgene.reduced # store sig snps to each gene in the elemnt
    
  } #for g in table_gene 
  
  # transform the outgred.df[[g]] list into a data frame
  #outbygene.df <- do.call(rbind, outbygene) 
  #outchr[[ch]] <- outbygene.df
  outchr[[ch]] <- do.call(rbind, outbygene) 

} #for ch



# transform the outgred.df[[g]] list into a data frame
outchr.df <- ldply(outchr, data.frame)
outchr.df <- outchr.df[order(as.numeric(outchr.df[,"snp.chr"]), as.numeric(outchr.df[,"snp.pos"])),]
outchr.df <- unique(outchr.df)


# SHARED GENES
# load genes from PCC results
load("/home/maize.vib/PUBBLICATION/7.RNA.pheno.cor/overlapping.genes.Rdata")

# create a dataframe with for each gene --> which phenotype has a correlation with
genebyph <- ldply (overlap, data.frame)

names(genebyph) <- c("phenotype", "gene")
genebyph$"pcc.LER" <- 0
genebyph$"pcc.LED" <- 0
genebyph$"pcc.leaf.length" <- 0
genebyph$"pcc.DZ.size" <- 0

for(ph in c("LER", "LED", "DZ.size", "leaf.length")){
  genebyph[which(genebyph[ , "gene"] %in% overlap[[ph]]), paste0("pcc.", ph)] <- 1
}

genebyph <- genebyph[ , -1]
genebyph <- genebyph[unique(genebyph$"gene"), ]
genebyph$"tot.pheno" <- rowSums(genebyph[ , 2:length(names(genebyph))])

write.table(genebyph, file = "/home/maize.vib/PUBBLICATION/7.RNA.pheno.cor/sig.genes.all.2.txt")


# merge with GWAS data 
outchr.df$"pcc.LER" <- 0
outchr.df$"pcc.LED" <- 0
outchr.df$"pcc.leaf.length" <- 0
outchr.df$"pcc.DZ.size" <- 0
for(ph in c("LER", "LED", "DZ.size", "leaf.length")){
  outchr.df[which(outchr.df[ , "gene"] %in% overlap[[ph]]), paste0("pcc.", ph)] <- 1
}


#sharedgLD = gwa.pos[which(gwa.pos[,"gene"] %in% pcc.pos.reduced[,"gene"]), "gene"]
#sharedgLD = unique(sharedgLD) 


# save min.pv.snps and min.pv.gene in the eqtl result
#min.pv.snsp <- as.data.frame(me.all.pca10.2$all$min.pv.snps)
min.pv.snsp <- as.data.frame(me$all$min.pv.snps)
names(min.pv.snsp) <- "min.pv.snps"
min.pv.snsp$"snps" <- rownames(min.pv.snsp)
outchr.df <- merge.data.frame(outchr.df, min.pv.snsp, by = "snps")

#min.pv.gene <- as.data.frame(me.all.pca10.2$all$min.pv.gene)
min.pv.gene <- as.data.frame(me$all$min.pv.gene)
names(min.pv.gene) <- "min.pv.gene"
min.pv.gene$"gene" <- rownames(min.pv.gene)
outchr.df <- merge.data.frame(outchr.df, min.pv.gene, by = "gene")

outchr.df <- outchr.df[order(as.numeric(outchr.df[,"snp.chr"]), as.numeric(outchr.df[,"snp.pos"])), ]


# calculate MAF
##########################
# load original geno file, to apply the MAF 
SNP_file_name = "/home/maize.vib/genotypic_data_2/genotypes.chr.all.filtered.biallelic.txt"
snps = SlicedData$new() 
snps$fileDelimiter = "\t"      # the TAB character
snps$fileOmitCharacters = "NA" # denote missing values;
snps$fileSkipRows = 1          # one row of column labels
snps$fileSkipColumns = 1       # one column of row labels
snps$fileSliceSize = 2000      # read file in pieces of 2,000 rows
snps$LoadFile( SNP_file_name )

eqtl.backup <- outchr.df # just in case, but the FDR is ok
outchr.df <- outchr.df %>% filter(FDR < 0.01)
maf <- unlist(lapply(seq_len(length(snps)), function(sl) {
  slice <- snps[[sl]]
  af <- rowMeans(slice, na.rm = TRUE) / 2
  pmin(af, 1 - af)
}))

maf_df <- data.frame(snp = snps$GetAllRowNames(), maf = maf, stringsAsFactors = FALSE) %>% tbl_df
head(maf_df)
#colnames(eqtl)[2] <- "snps"
eqtl_maf <- outchr.df %>% left_join(maf_df, by = c(snps = "snp"))
orig_cor <- with(eqtl_maf, cor.test(maf, abs(beta), method = "spearman", exact = FALSE))


# SELECT GWAS AND PCC eQTLs
eqtl <- eqtl_maf
gwa.eqtl <- eqtl[which(apply(eqtl[ , 9:12], 1, sum) > 0), ]
pcc.eqtl <- eqtl[which(apply(eqtl[ , 14:17], 1, sum) > 0), ]
int.eqtl <- eqtl[which(apply(eqtl[ , 9:12], 1, sum) > 0 & apply(eqtl[ , 14:17], 1, sum) > 0), ]

int.eqtl$"tot.gwa" <- NA
int.eqtl$"tot.gwa" <- apply(int.eqtl[ , 9:12], 1, sum)
int.eqtl$"tot.pcc" <- NA
int.eqtl$"tot.pcc" <- apply(int.eqtl[ , 14:17], 1, sum)

#print how many genes and SNPs each phenotype and blablabla
print("tot eQTLs for PCC-eQTLs"); for (i in phenotypes) {print(c(i, nrow(unique(pcc.eqtl[which(pcc.eqtl[ , paste0("pcc.", i)] > 0) , ]))))}
print("unique genes for PCC-eQTLs"); for (i in phenotypes) {print(c(i, length(unique(pcc.eqtl[which(pcc.eqtl[ , paste0("pcc.", i)] > 0) , ]$gene))))}
print("unique SNPs for PCC-eQTLs"); for (i in phenotypes) {print(c(i, length(unique(pcc.eqtl[which(pcc.eqtl[ , paste0("pcc.", i)] > 0) , ]$snps))))}

print("tot eQTLs for GWA-eQTLs"); for (i in phenotypes) {print(c(i, nrow(unique(gwa.eqtl[which(gwa.eqtl[ , paste0("gwa.", i)] > 0) , ]))))}
print("unique genes for GWA-eQTLs"); for (i in phenotypes) {print(c(i, length(unique(gwa.eqtl[which(gwa.eqtl[ , paste0("gwa.", i)] > 0) , ]$gene))))}
print("unique SNPs for GWA-eQTLs"); for (i in phenotypes) {print(c(i,length(unique(gwa.eqtl[which(gwa.eqtl[ , paste0("gwa.", i)] > 0) , ]$snps))))}

print("tot eQTLs for PCC-GWA-eQTLs"); for (i in phenotypes) {print(c(i, nrow(unique(pcc.gwa.eqtl[which(pcc.gwa.eqtl[ , paste0("pcc.", i)] > 0 | pcc.gwa.eqtl[ , paste0("gwa.", i)] > 0) , ]))))}
print("unique genes for PCC-GWA-eQTLs"); for (i in phenotypes) {print(c(i, length(unique(pcc.gwa.eqtl[which(pcc.gwa.eqtl[ , paste0("pcc.", i)] > 0 | pcc.gwa.eqtl[ , paste0("gwa.", i)] > 0) , ]$gene))))}
print("unique SNPs for PCC-GWA-eQTLs"); for (i in phenotypes) {print(c(i, length(unique(pcc.gwa.eqtl[which(pcc.gwa.eqtl[ , paste0("pcc.", i)] > 0 | pcc.gwa.eqtl[ , paste0("gwa.", i)] > 0) , ]$snps))))}

print("tot eQTLs INTERSECTION PCC-GWA-eQTLs"); for (i in phenotypes) {print(c(i, nrow(unique(int.eqtl[which(int.eqtl[ , paste0("pcc.", i)] > 0 & int.eqtl[ , paste0("gwa.", i)] > 0) , ]))))}
print("unique genes INTERSECTION PCC-GWA-eQTLs"); for (i in phenotypes) {print(c(i, length(unique(int.eqtl[which(int.eqtl[ , paste0("pcc.", i)] > 0 & int.eqtl[ , paste0("gwa.", i)] > 0) , ]$gene))))}
print("unique SNPs INTERSECTION PCC-GWA-eQTLs"); for (i in phenotypes) {print(c(i, length(unique(int.eqtl[which(int.eqtl[ , paste0("pcc.", i)] > 0 & int.eqtl[ , paste0("gwa.", i)] > 0) , ]$snps))))}





# SAVE RESULTS  
write.table(outchr.df, file = paste0(workdir, "/eqtl.all.red.2.txt"), quote = F, sep = "\t", row.names = F)
#write.table(eqtl_maf, file = paste0(workdir, "/eqtl.all.red.2.maf.txt"), quote = F, sep = "\t", row.names = F)
save(me, eqtl_maf, pcc.eqtl, gwa.eqtl, int.eqtl, pcc.gwa.eqtl, file=paste0(workdir, "/", "me.all.pca10.all.reduced.2.Rdata"))

pcc.eqtl <- pcc.eqtl[ , -c(9:13)] # remove gwa columns and LDpeak column
gwa.eqtl <- gwa.eqtl[ , -c(13:17)] # remove gwa columns and LDpeak column

write.table(pcc.eqtl, file = paste0(workdir, "/pcc.eqtl.red.txt"), quote = F, sep = "\t", row.names = F)
write.table(gwa.eqtl, file = paste0(workdir, "/gwa.eqtl.red.txt"), quote = F, sep = "\t", row.names = F)

  
# PLOT results
# not shown in main text

## Plot the FDR 
###############
#eqtl reduced
pdf("/home/maize.vib/eqtl_results_2/eqtl_graphs/FDR.all.red.2.pdf")
hist(eqtl_maf$FDR)
dev.off()
## pcc.eqtl
pdf("/home/maize.vib/eqtl_results_2/eqtl_graphs/FDR.pcc.eqtl.pdf")
hist(pcc.eqtl$FDR, nclass = 20)
dev.off()
## gwa.eqtl
pdf("/home/maize.vib/eqtl_results_2/eqtl_graphs/FDR.gwa.eqtl.pdf")
hist(gwa.eqtl$FDR, nclass = 20)
dev.off()


## Plot the pvalue 
###############
#eqtl reduced
pdf("/home/maize.vib/eqtl_results_2/eqtl_graphs/pvalue.all.red.2.pdf")
hist(eqtl_maf$pvalue)
dev.off()
## pcc.eqtl
pdf("/home/maize.vib/eqtl_results_2/eqtl_graphs/pvalue.pcc.eqtl.pdf")
hist(pcc.eqtl$pvalue, nclass = 20)
dev.off()
## gwa.eqtl
pdf("/home/maize.vib/eqtl_results_2/eqtl_graphs/pvalue.gwa.eqtl.pdf")
hist(gwa.eqtl$pvalue, nclass = 20)
dev.off()


# already done !!
thr.cis = 1000000 # 1Mb
int.eqtl$"cistrans" <- "NA"
int.eqtl[which(int.eqtl$"snp.chr" != int.eqtl$"gene.chr"), ]$"cistrans" <- "distant_chr"
int.eqtl[which(int.eqtl$"snp.chr" == int.eqtl$"gene.chr" & abs(int.eqtl$"snp.pos" - int.eqtl$"gene.start") > thr.cis), ]$"cistrans" <- "distant"
int.eqtl[which(int.eqtl$"snp.chr" == int.eqtl$"gene.chr" & abs(int.eqtl$"snp.pos" - int.eqtl$"gene.start") <= thr.cis), ]$"cistrans" <- "local"
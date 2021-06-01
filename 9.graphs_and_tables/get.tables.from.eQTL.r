#######################################
#STATISTICS eQTL_file_name
#######################################
library(plyr)
library(ggplot2)
library(dplyr)
library(data.table)


#maindir<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/MAIZE_ghent/analyses/11.manage.eQTL.out"
workdir<-"/home/maize.vib/eqtl_results_2/eqtl_graphs"
eQTL_file_name = "/home/maize.vib/eqtl_results_2/me.all.pca10.all.reduced.3.Rdata"
traits<-c("leaf.length", "DZ.size", "LER", "LED")

biparental <- c("RIL100","RIL103","RIL104","RIL105","RIL106","RIL107","RIL109","RIL11","RIL111","RIL113","RIL114","RIL115","RIL116","RIL117","RIL118","RIL119","RIL120","RIL121","RIL122","RIL123","RIL124","RIL125",
                "RIL126","RIL127","RIL128","RIL13","RIL130","RIL131","RIL132","RIL133","RIL134","RIL135","RIL136","RIL137","RIL138","RIL139","RIL140","RIL141","RIL16","RIL17","RIL18","RIL19","RIL23","RIL25","RIL26","RIL27",
                "RIL29","RIL30","RIL32","RIL33","RIL36","RIL37","RIL38","RIL39","RIL40","RIL42","RIL43","RIL44","RIL45","RIL46","RIL47","RIL48","RIL49","RIL50","RIL51","RIL52","RIL53","RIL54","RIL55","RIL57","RIL58","RIL59","RIL60","RIL61","RIL62","RIL64","RIL65","RIL66","RIL67","RIL69","RIL70","RIL71","RIL72","RIL73","RIL74","RIL75","RIL76","RIL80","RIL81","RIL82","RIL83","RIL84","RIL85","RIL88","RIL89","RIL90","RIL91","RIL92","RIL93","RIL94","RIL95","RIL96","RIL99")

magic <- c("RIL_8W_1", "RIL_8W_10", "RIL_8W_11","RIL_8W_12","RIL_8W_13","RIL_8W_14","RIL_8W_15","RIL_8W_16","RIL_8W_18","RIL_8W_19","RIL_8W_2","RIL_8W_20","RIL_8W_22","RIL_8W_23","RIL_8W_24","RIL_8W_26","RIL_8W_27","RIL_8W_28","RIL_8W_29","RIL_8W_3","RIL_8W_30","RIL_8W_31","RIL_8W_32","RIL_8W_33","RIL_8W_34","RIL_8W_35","RIL_8W_36","RIL_8W_37","RIL_8W_38","RIL_8W_39","RIL_8W_4","RIL_8W_40","RIL_8W_41","RIL_8W_42","RIL_8W_44","RIL_8W_45","RIL_8W_46",
           "RIL_8W_47","RIL_8W_48","RIL_8W_49","RIL_8W_5","RIL_8W_50","RIL_8W_51","RIL_8W_52","RIL_8W_53","RIL_8W_54","RIL_8W_55","RIL_8W_57","RIL_8W_58","RIL_8W_59","RIL_8W_6","RIL_8W_60","RIL_8W_61","RIL_8W_62","RIL_8W_63","RIL_8W_64","RIL_8W_65","RIL_8W_66","RIL_8W_67","RIL_8W_68","RIL_8W_69","RIL_8W_7","RIL_8W_70","RIL_8W_71","RIL_8W_72","RIL_8W_73","RIL_8W_74","RIL_8W_75","RIL_8W_76","RIL_8W_77",
           "RIL_8W_78","RIL_8W_79","RIL_8W_8","RIL_8W_80","RIL_8W_81","RIL_8W_82","RIL_8W_83","RIL_8W_84","RIL_8W_85","RIL_8W_86","RIL_8W_87","RIL_8W_88","RIL_8W_89","RIL_8W_9","RIL_8W_90","RIL_8W_91","RIL_8W_92","RIL_8W_93","RIL_8W_94","RIL_8W_95","RIL_8W_96","RIL_8W_97","RIL_8W_98","RIL_8W_99")

parents <- c("B73", "H99", "A632", "CML91", "F7", "HP301", "Mo17","W153R")

traits<-c("leaf.length", "DZ.size", "LER", "LED")

# fake variables
workdir = "/Users/mara/work_analysis"
setwd(workdir)
SNP_file_name = "/Users/mara/work_analysis/maize.vib/inputs/genotypes.chr.all.filtered.biallelic.txt"
expression_file_name = "/Users/mara/work_analysis/maize.vib/inputs/all.tpm.chr.0.05var.txt"
pheno_file_name = "/Users/mara/work_analysis/maize.vib/inputs/allpheno.txt" #"/home/maize.vib/phenotypes/allpheno.txt"
eQTL_file_name = "/Users/mara/work_analysis/maize.vib/me.all.pca10.all.reduced.3.Rdata"

setwd(workdir)


load(eQTL_file_name)
# there is the object eqtl and 
# eqtl.pheno with the counts of genes each SNP



setwd("/home/maize.vib/eqtl_results_2/tables_for_graphs")
####################################
# tot eQTL get the TABLE STATISTICS
####################################
#start getting stats from GWA-eQTL
#neqtl<-nrow(gwa.pos)
neqtl<-nrow(gwa.eqtl)
uniq.s<-length(unique(gwa.eqtl[,"snps"]))
uniq.g<-length(unique(gwa.eqtl[,"gene"]))

tinfo<-data.frame(table(gwa.eqtl$cistrans),(table(gwa.eqtl$cistrans)/neqtl)*100)
tinfo<-tinfo[,-3]
colnames(tinfo)<-c("class", "n", "%")

gwa.g <- table(gwa.eqtl$gene)
gwa.s <- table(gwa.eqtl$snps)

#get number of genes by snps and the opposite
genebys <-	max(gwa.g)
snpsbyg <- max(gwa.s)

#get number of genes by snps and the opposite
#genebys<-max(gwa.s[,2])
#snpsbyg<-max(gwa.g[,2])

#prepare a txt file for output
out = file(paste("GWA.eqtl.output.txt", sep=""), 'a')

#write stats on file
write(paste(neqtl,"eQTL"), out, append=TRUE)
write(paste(uniq.s,"SNPs"), out, append=TRUE)
write(paste(uniq.g,"genes"), out, append=TRUE)
write(paste(genebys,"maximum genes by SNP"), out, append=TRUE)
write(paste(snpsbyg,"maximum SNPs per gene"), out, append=TRUE)

write("cis/trans information below", out, append=TRUE)
write(tinfo[,1], out, append=TRUE)
write(tinfo[,2], out, append=TRUE)
write(tinfo[,3], out, append=TRUE)
close(out)

#get out info useful for plots
gwa.eqtl.s<-as.data.frame(gwa.s)[,2]
gwa.eqtl.g<-as.data.frame(gwa.g)[,2]

#write out supplemental tables
#write.table(gwa.eqtl, file=paste("GWA.eQTL.table.txt", sep=""), sep="\t", quote=F, row.names=F)


####################################
# get the TABLE STATISTICS
####################################
# for each trait
################
for (f in 1:length(traits)){
  tmptrait<-traits[f]
  print(tmptrait)
  #get everything in memory
  tmp <- paste0("gwa.", tmptrait)
  gwa.pheno <- gwa.eqtl[which(gwa.eqtl[ , tmp] > 0) , ]
  #	load(filen[grep(tmptrait, filen)][1])
  #	load(filen[grep(tmptrait, filen)][2])
  
  ####################################
  #start getting stats from GWA-eQTL
  neqtl<-nrow(gwa.pheno)
  uniq.s<-length(unique(gwa.pheno[,"snps"]))
  uniq.g<-length(unique(gwa.pheno[,"gene"]))
  
  #prova <- count(gwa.leaf.length, "gene")
  
  tinfo<-data.frame(table(gwa.pheno$cistrans),(table(gwa.pheno$cistrans)/neqtl)*100)
  tinfo<-tinfo[,-3]
  colnames(tinfo)<-c("class", "n", "%")
  
  gwa.g <- count(gwa.pheno, "gene")
  gwa.s <- count(gwa.pheno, "snps")
  
  #get number of genes by snps and the opposite
  genebys<-max(gwa.s[,2])
  snpsbyg<-max(gwa.g[,2])
  
  #prepare a txt file for output
  out = file(paste(tmptrait, ".gwa.pheno.output.txt", sep=""), 'a')
  
  #write stats on file
  write(paste(tmptrait, "- GWA eQTL"), out, append=TRUE)
  write(paste(neqtl,"eQTL"), out, append=TRUE)
  write(paste(uniq.s,"SNPs"), out, append=TRUE)
  write(paste(uniq.g,"genes"), out, append=TRUE)
  write(paste(genebys,"maximum genes by SNP"), out, append=TRUE)
  write(paste(snpsbyg,"maximum SNPs per gene"), out, append=TRUE)
  
  write("cis/trans information below - cis/trans/tran_chr", out, append=TRUE)
  write(tinfo[,1], out, append=TRUE)
  write(tinfo[,2], out, append=TRUE)
  write(tinfo[,3], out, append=TRUE)
  close(out)
  
  #get out info useful for plots
  #gwa.gbyslist[[f]]<-gwa.s[,2]
  #gwa.sbyglist[[f]]<-gwa.g[,2]
  gwa.gbyslist[[f]]<-as.data.frame(gwa.s)[,2]
  gwa.sbyglist[[f]]<-as.data.frame(gwa.g)[,2]
  
  #write out supplemental tables
  #write.table(gwa.pos, file=paste(tmptrait, ".gwa.pheno.table.txt", sep=""), sep="\t", quote=F, row.names=F)
  
  ####################################
  #get stats from PCC-eQTL
  tmp <- paste0("pcc.", tmptrait)
  pcc.pheno <- pcc.eqtl[which(pcc.eqtl[ , tmp] > 0) , ]
  
  neqtl<-nrow(pcc.pheno)
  uniq.s<-length(unique(pcc.pheno[,"snps"]))
  uniq.g<-length(unique(pcc.pheno[,"gene"]))
  
  tinfo<-data.frame(table(pcc.pheno$cistrans),(table(pcc.pheno$cistrans)/neqtl)*100)
  tinfo<-tinfo[,-3]
  colnames(tinfo)<-c("class", "n", "%")
  
  
  #get number of genes by snps and the opposite
  byg <- count(pcc.pheno, "gene")
  bys<-count(pcc.pheno, "snps")
  
  #byg<-split(pcc.pheno, pcc.pheno[,"gene"])
  #bys<-split(pcc.pheno, pcc.pheno[,"snps"])
  #byg<-sort(unlist(lapply(byg, nrow)))
  #bys<-sort(unlist(lapply(bys, nrow)))
  
  #prepare a txt file for output
  out = file(paste(tmptrait, ".PCC.eqtl.output.txt", sep=""), 'a')
  
  #write stats on file
  write(paste(tmptrait), out, append=TRUE)
  write(paste(neqtl,"eQTL"), out, append=TRUE)
  write(paste(uniq.s,"SNPs"), out, append=TRUE)
  write(paste(uniq.g,"genes"), out, append=TRUE)
  write(paste(max(bys[,2]),"maximum genes by SNP"), out, append=TRUE)
  write(paste(max(byg[,2]),"maximum SNPs per gene"), out, append=TRUE)
  
  write("cis/trans information below - cis/trans/tran_chr", out, append=TRUE)
  write(tinfo[,1], out, append=TRUE)
  write(tinfo[,2], out, append=TRUE)
  write(tinfo[,3], out, append=TRUE)
  close(out)
  
  
  #get out info useful for plots
  pcc.gbyslist[[f]]<-bys[ ,2]
  pcc.sbyglist[[f]]<-byg[ ,2]
  
  #write out supplemental tables
  #write.table(pcc.pheno, file=paste(tmptrait, ".PCC.eQTL.table.txt", sep=""), sep="\t", quote=F, row.names=F)
  
  #get shared genes
  #	sharedg<-pcc.pos.reduced[which(unique(pcc.pos.reduced$gene) %in% unique(gwa.pos$gene)), 4:7]
  sharedg<-pcc.pheno[which(unique(pcc.pheno$gene) %in% unique(gwa.pheno$gene)), c("gene", "gene.chr", "gene.start", "gene.end")]
  print(paste(tmptrait, "features", nrow(sharedg), "unique genes"))
  shareglist[[f]]<-sharedg
  write.table(sharedg, file=paste(tmptrait, ".shared.genes.txt", sep=""), sep="\t", quote=F, row.names=F)
} # for f
####################################

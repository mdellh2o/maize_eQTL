#######################################
# CREATE FILE for PLOTs
#######################################
library(plyr)
library(dplyr)
library(data.table)


workdir<-"/home/maize.vib/eqtl_results_2/eqtl_graphs"
eQTL_file_name = "/home/maize.vib/eqtl_results_2/me.all.pca10.all.reduced.2.Rdata"
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
eQTL_file_name = "/Users/mara/work_analysis/maize.vib/me.all.pca10.all.reduced.2.Rdata"


setwd(workdir)
load(eQTL_file_name)
# there is the object eqtl_maf 
# and add genetic positions!

gene_location_rdata = "/home/maize.vib/annotation/gene_positions_V4.Rdata" #the same of the .txt, object called genepos
load(gene_location_rdata) #object in bed file --> genepos

genepos <- fread("/Users/mara/work_analysis/gene_position.txt", data.table=F)
genepos <- genepos[, 1:4]

# add gene positions at eqtl file
names(genepos) <- c("gene", "gene.chr", "gene.start", "gene.end")
eqtl <- merge(eqtl_maf, genepos, by = "gene")


# calculate CIS-TRANS
######################
thr.cis = 1000000 # 1Mb
eqtl$"cistrans" <- "NA"
eqtl[which(eqtl$"snp.chr" != eqtl$"gene.chr"), ]$"cistrans" <- "trans_chr"
eqtl[which(eqtl$"snp.chr" == eqtl$"gene.chr" & abs(eqtl$"snp.pos" - eqtl$"gene.start") > thr.cis), ]$"cistrans" <- "trans"
eqtl[which(eqtl$"snp.chr" == eqtl$"gene.chr" & abs(eqtl$"snp.pos" - eqtl$"gene.start") <= thr.cis), ]$"cistrans" <- "cis"

eqtl$"tot.gwa" <- NA
eqtl$"tot.gwa" <- apply(eqtl[ , 9:12], 1, sum)
eqtl$"tot.pcc" <- NA
eqtl$"tot.pcc" <- apply(eqtl[ , 14:17], 1, sum)

eqtl<- eqtl[ , -13] # remove column "LDpeak"

#gwa.eqtl <- eqtl[which(apply(eqtl[ , 9:12], 1, sum) > 0), ]
#pcc.eqtl <- eqtl[which(apply(eqtl[ , 14:17], 1, sum) > 0), ]
#int.eqtl <- eqtl[which(apply(eqtl[ , 9:12], 1, sum) > 0 & apply(eqtl[ , 14:17], 1, sum) > 0), ]
gwa.eqtl <- eqtl[which(eqtl$"tot.gwa" > 0), ]
pcc.eqtl <- eqtl[which(eqtl$"tot.pcc" > 0), ]
int.eqtl <- eqtl[which(eqtl$"tot.gwa" > 0 & eqtl$"tot.pcc" > 0), ]

gwa.eqtl$gene <- droplevels(gwa.eqtl$gene)
#gwa.eqtl$snps <- droplevels(gwa.eqtl$snps)
pcc.eqtl$gene <- droplevels(pcc.eqtl$gene)
#pcc.eqtl$snps <- droplevels(pcc.eqtl$snps)
int.eqtl$gene <- droplevels(int.eqtl$gene)
#int.eqtl$snps <- droplevels(int.eqtl$snps)



# COUNT number of transcripts for EACH TRAIT! 
#############################################
gwa.leaf.length <- gwa.eqtl[which(gwa.eqtl$gwa.leaf.length > 0), ]
gwa.snp.leaf.length <- plyr::count(gwa.leaf.length, "snps")
gwa.snp.leaf.length <- unique(merge(gwa.snp.leaf.length, gwa.eqtl[ , c("snps", "maf")]))
gwa.snp.leaf.length$"trait" <- "leaf.length"
gwa.snp.leaf.length$"eqtl" <- "GWA.eQTL"

gwa.LER <- gwa.eqtl[which(gwa.eqtl$gwa.LER > 0), ]
gwa.snp.LER <- plyr::count(gwa.LER, "snps")
gwa.snp.LER <- unique(merge(gwa.snp.LER, gwa.eqtl[ , c("snps", "maf")]))
gwa.snp.LER$"trait" <- "LER"
gwa.snp.LER$"eqtl" <- "GWA.eQTL"

gwa.LED <- gwa.eqtl[which(gwa.eqtl$gwa.LED > 0), ]
gwa.snp.LED <- plyr::count(gwa.LED, "snps")
gwa.snp.LED <- unique(merge(gwa.snp.LED, gwa.eqtl[ , c("snps", "maf")]))
gwa.snp.LED$"trait" <- "LED"
gwa.snp.LED$"eqtl" <- "GWA.eQTL"

gwa.DZ.size <- gwa.eqtl[which(gwa.eqtl$gwa.DZ.size > 0), ]
gwa.snp.DZ.size <- plyr::count(gwa.DZ.size, "snps")
gwa.snp.DZ.size <- unique(merge(gwa.snp.DZ.size, gwa.eqtl[ , c("snps", "maf")]))
gwa.snp.DZ.size$"trait" <- "DZ.size"
gwa.snp.DZ.size$"eqtl" <- "GWA.eQTL"


# COUNT number of transcripts for EACH TRAIT! 
#############################################
pcc.leaf.length <- pcc.eqtl[which(pcc.eqtl$pcc.leaf.length > 0), ]
pcc.snp.leaf.length <- plyr::count(pcc.leaf.length, "snps")
pcc.snp.leaf.length <- unique(merge(pcc.snp.leaf.length, pcc.eqtl[ , c("snps", "maf")]))
pcc.snp.leaf.length$"trait" <- "leaf.length"
pcc.snp.leaf.length$"eqtl" <- "pcc.eQTL"

pcc.LER <- pcc.eqtl[which(pcc.eqtl$pcc.LER > 0), ]
pcc.snp.LER <- plyr::count(pcc.LER, "snps")
pcc.snp.LER <- unique(merge(pcc.snp.LER, pcc.eqtl[ , c("snps", "maf")]))
pcc.snp.LER$"trait" <- "LER"
pcc.snp.LER$"eqtl" <- "pcc.eQTL"

pcc.LED <- pcc.eqtl[which(pcc.eqtl$pcc.LED > 0), ]
pcc.snp.LED <- plyr::count(pcc.LED, "snps")
pcc.snp.LED <- unique(merge(pcc.snp.LED, pcc.eqtl[ , c("snps", "maf")]))
pcc.snp.LED$"trait" <- "LED"
pcc.snp.LED$"eqtl" <- "pcc.eQTL"

pcc.DZ.size <- pcc.eqtl[which(pcc.eqtl$pcc.DZ.size > 0), ]
pcc.snp.DZ.size <- plyr::count(pcc.DZ.size, "snps")
pcc.snp.DZ.size <- unique(merge(pcc.snp.DZ.size, pcc.eqtl[ , c("snps", "maf")]))
pcc.snp.DZ.size$"trait" <- "DZ.size"
pcc.snp.DZ.size$"eqtl" <- "pcc.eQTL"




# COUNT TRANSCRIPTS/SNP all gwa.eqtl and pcc.eqtl 
eqtl.leaf.length <- eqtl[which(eqtl$gwa.leaf.length > 0 | eqtl$pcc.leaf.length > 0), ]
eqtl.leaf.length <- plyr::count(eqtl.leaf.length, "snps")
eqtl.leaf.length <- unique(merge(eqtl.leaf.length, eqtl[ , c("snps", "maf")]))
eqtl.leaf.length$"trait" <- "leaf.length"
eqtl.leaf.length$"eqtl" <- "NA"
eqtl.leaf.length[which(eqtl.leaf.length$"snps" %in% eqtl[eqtl$gwa.leaf.length > 0, "snps"]), ]$"eqtl" <- "GWA.eQTL"
eqtl.leaf.length[which(eqtl.leaf.length$"snps" %in% eqtl[eqtl$pcc.leaf.length > 0, "snps"]), ]$"eqtl" <- "PCC.eQTL"
eqtl.leaf.length[which(eqtl.leaf.length$"snps" %in% eqtl[(eqtl$gwa.leaf.length > 0 & eqtl$pcc.leaf.length > 0), "snps"]), ]$"eqtl" <- "GWA.PCC.eQTL"

eqtl.LER <- eqtl[which(eqtl$gwa.LER > 0 | eqtl$pcc.LER > 0), ]
eqtl.LER <- plyr::count(eqtl.LER, "snps")
eqtl.LER <- unique(merge(eqtl.LER, eqtl[ , c("snps", "maf")]))
eqtl.LER$"trait" <- "LER"
eqtl.LER$"eqtl" <- "NA"
eqtl.LER[which(eqtl.LER$"snps" %in% eqtl[eqtl$gwa.LER > 0, "snps"]), ]$"eqtl" <- "GWA.eQTL"
eqtl.LER[which(eqtl.LER$"snps" %in% eqtl[eqtl$pcc.LER > 0, "snps"]), ]$"eqtl" <- "PCC.eQTL"
eqtl.LER[which(eqtl.LER$"snps" %in% eqtl[(eqtl$gwa.LER > 0 & eqtl$pcc.LER > 0), "snps"]), ]$"eqtl" <- "GWA.PCC.eQTL"

eqtl.LED <- eqtl[which(eqtl$gwa.LED > 0 | eqtl$pcc.LED > 0), ]
eqtl.LED <- plyr::count(eqtl.LED, "snps")
eqtl.LED <- unique(merge(eqtl.LED, eqtl[ , c("snps", "maf")]))
eqtl.LED$"trait" <- "LED"
eqtl.LED$"eqtl" <- "NA"
eqtl.LED[which(eqtl.LED$"snps" %in% eqtl[eqtl$gwa.LED > 0, "snps"]), ]$"eqtl" <- "GWA.eQTL"
eqtl.LED[which(eqtl.LED$"snps" %in% eqtl[eqtl$pcc.LED > 0, "snps"]), ]$"eqtl" <- "PCC.eQTL"
eqtl.LED[which(eqtl.LED$"snps" %in% eqtl[(eqtl$gwa.LED > 0 & eqtl$pcc.LED > 0), "snps"]), ]$"eqtl" <- "GWA.PCC.eQTL"

eqtl.DZ.size <- eqtl[which(eqtl$gwa.DZ.size > 0 | eqtl$pcc.DZ.size > 0), ]
eqtl.DZ.size <- plyr::count(eqtl.DZ.size, "snps")
eqtl.DZ.size <- unique(merge(eqtl.DZ.size, eqtl[ , c("snps", "maf")]))
eqtl.DZ.size$"trait" <- "DZ.size"
eqtl.DZ.size$"eqtl" <- "NA"
eqtl.DZ.size[which(eqtl.DZ.size$"snps" %in% eqtl[eqtl$gwa.DZ.size > 0, "snps"]), ]$"eqtl" <- "GWA.eQTL"
eqtl.DZ.size[which(eqtl.DZ.size$"snps" %in% eqtl[eqtl$pcc.DZ.size > 0, "snps"]), ]$"eqtl" <- "PCC.eQTL"
eqtl.DZ.size[which(eqtl.DZ.size$"snps" %in% eqtl[(eqtl$gwa.DZ.size > 0 & eqtl$pcc.DZ.size > 0), "snps"]), ]$"eqtl" <- "GWA.PCC.eQTL"

eqtl.pheno <- rbind(eqtl.leaf.length, eqtl.LER, eqtl.LED, eqtl.DZ.size)
eqtl.pheno[which(eqtl.pheno$"snps" %in% eqtl[(eqtl$tot.gwa > 1 | eqtl$tot.pcc > 1), "snps"]), ]$"trait" <- "more"

save(me, eqtl, eqtl.pheno, gwa.eqtl, pcc.eqtl, int.eqtl, eqtl.DZ.size, eqtl.leaf.length, eqtl.LER, eqtl.LED, 
     file=paste0(workdir, "/me.all.pca10.all.reduced.3.Rdata"))




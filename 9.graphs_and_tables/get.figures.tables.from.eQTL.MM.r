#######################################
# MAKE PLOTS and TABLEs
#######################################
library(plyr)
library(ggplot2)
library(dplyr)
library(data.table)


#maindir<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/MAIZE_ghent/analyses/11.manage.eQTL.out"
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
# and add gentic positions!

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

gwa.eqtl <- eqtl[which(apply(eqtl[ , 9:12], 1, sum) > 0), ]
pcc.eqtl <- eqtl[which(apply(eqtl[ , 14:17], 1, sum) > 0), ]
int.eqtl <- eqtl[which(apply(eqtl[ , 9:12], 1, sum) > 0 & apply(eqtl[ , 14:17], 1, sum) > 0), ]

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


# all gwa.eqtl and pcc.eqtl 
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


######################################################
# GRAPHS
######################################################
#load colors
#load("/home/maize.vib/PUBBLICATION/colors.by.phenotype.Rdata")
load("/Users/mara/work_analysis/colors.by.phenotype.Rdata")
phenocol<-phenocol[c(2,1,4,3),]
stopifnot(all(phenocol[,1]==traits))
phenocol <- rbind(phenocol, c("more", "#000000"))
traits<-c("leaf.length", "DZ.size", "LER", "LED", "more")

# GRAPH FREQ/MAF by phenotype
#datasnp <- rbind(eqtl.DZ.size, eqtl.leaf.length, eqtl.LER, eqtl.LED)
datasnp <- eqtl.pheno
datasnp <- datasnp[order(datasnp$maf), ]

plot.maf <- ggplot(datasnp, aes(x=maf, y=freq)) +
  geom_point(aes(color=trait, shape = eqtl)) + scale_color_manual(breaks = traits, values = phenocol[ ,"color"]) +
  scale_shape_manual(breaks=c("GWA.eQTL", "GWA.PCC.eQTL", "PCC.eQTL"), values=c(16, 17, 18)) +
  labs(color = "phenotype")
plot.maf
ggsave("maf_by_pheno.tiff", plot = plot.maf)



# GRAPH FREQ/MAF by LOCATION
snp.freq <- plyr::count(eqtl, "snps")
names(snp.freq)[2]<- "num.genes"
snp.freq <- unique(merge(snp.freq, eqtl[ , c("snps", "maf", "cistrans")], by = "snps"))
snp.freq$"eqtl" <- "eqtl"
snp.freq[which(snp.freq$"snps" %in% pcc.eqtl$snps), ]$eqtl <- "PCC.eQTL"
snp.freq[which(snp.freq$"snps" %in% gwa.eqtl$snps), ]$eqtl <- "GWA.eQTL"
snp.freq[which(snp.freq$"snps" %in% pcc.eqtl$snps & snp.freq$"snps" %in% gwa.eqtl$snps), ]$eqtl <- "GWA.eQTL and PCC.eQTL"


plot.eqtl.freq = ggplot(snp.freq, aes(x=maf, y=num_genes, group=cistrans, color=cistrans)) +
  geom_point(size = 0.5) +
  labs(x = "MAF", y = "num transcript/snp", fill = "eQTL location")
plot.snp.freq
ggsave("n.genes_cistrans.tiff", plot = plot.snp.freq)


# GRAPH FREQ/MAF by eQTL (+ GWA and PCC)
library(RColorBrewer) #Create a custom color scale
myColors <- as.character(c("eqtl"="black", "GWA.eQTL"="red", "PCC.eQTL"="darkgray", "GWA.eQTL and PCC.eQTL" = "green"))
#mySize <- as.character(c("eqtl"=0.5, "GWA.eQTL"=0.8, "PCC.eQTL"=0.7, "LED"=0.8))

target <- c("PCC.eQTL", "GWA.eQTL and PCC.eQTL", "GWA.eQTL")
gwa.pcc.snp <- snp.freq[which(snp.freq$eqtl != "eqtl"), ]
gwa.pcc.snp <- gwa.pcc.snp %>% arrange(factor(eqtl, levels = target))

myColors <- as.character(c("brown1", "blue", "antiquewhite4"))
plot.snp.eqtl = ggplot(gwa.pcc.snp, aes(x=maf, y=num.genes, group=eqtl, color=eqtl)) +
  geom_point(size = 0.6) +
  scale_color_manual(values=myColors) + 
  labs(x = "MAF", y = "num transcript/snp")
plot.snp.eqtl 

ggsave("n.genes.pcc.gwa.tiff", plot = plot.snp.freq)



######################################
#get variance explained by eQTL - BETA
# VIOLIN PLOTS
######################################
#library(vioplot)
setwd(maindir)

traits<-c("leaf.length", "DZ.size", "LER", "LED")

#load colors
load("/home/maize.vib/PUBBLICATION/colors.by.phenotype.Rdata")
phenocol<-phenocol[c(2,1,4,3),]
stopifnot(all(phenocol[,1]==traits))


gwa.violin <- list()
pcc.violin <- list()

for (f in 1:length(traits)){
  tmptrait<-traits[f]
  print(tmptrait)
  
  tmp <- paste0("gwa.", tmptrait)
  gwa.pheno <- eqtl[which(eqtl[ , tmp] > 0) , ]
  gwa.pheno<-gwa.pheno[order(gwa.pheno[,"cistrans"]), ]
  gwa.pheno$"cistrans" = factor(gwa.pheno$"cistrans")
  
  gwa.violin[[f]] <- ggplot(gwa.pheno, aes(x=cistrans, y=beta)) +
    geom_violin(fill = phenocol[f, 2], colour = "black") +
    labs(title = tmptrait, x = NULL, y = NULL, fill = "Tipologia")
  
  tmp <- paste0("pcc.", tmptrait)
  pcc.pheno <- eqtl[which(eqtl[ , tmp] > 0) , ]
  pcc.pheno<-pcc.pheno[order(pcc.pheno[,"cistrans"]), ]
  pcc.pheno$"cistrans" = factor(pcc.pheno$"cistrans")
  
  pcc.violin[[f]] <- ggplot(pcc.pheno, aes(x=cistrans, y=beta)) +
    geom_violin(fill = phenocol[f, 2], colour = "black") +
    labs(title = tmptrait, x = NULL, y = NULL, fill = "Tipologia")
}
library(patchwork)
plot.gwa.violin <- gwa.violin[[1]] | gwa.violin[[2]] | gwa.violin[[3]] | gwa.violin[[4]]
plot.gwa.violin
plot.pcc.violin <- pcc.violin[[1]] | pcc.violin[[2]] | pcc.violin[[3]] | pcc.violin[[4]]
plot.pcc.violin
plot.gwa.violin / plot.pcc.violin



############################  
# BOXPLOT (instead of violin) 
############################
#library(vioplot)
setwd(maindir)
traits<-c("leaf.length", "DZ.size", "LER", "LED")

#load colors
load("/home/maize.vib/PUBBLICATION/colors.by.phenotype.Rdata")
phenocol<-phenocol[c(2,1,4,3),]
stopifnot(all(phenocol[,1]==traits))

gwa.boxplot <- list()
pcc.boxplot <- list()

for (f in 1:length(traits)){
  tmptrait<-traits[f]
  print(tmptrait)
  
  tmp <- paste0("gwa.", tmptrait)
  gwa.pheno <- eqtl[which(eqtl[ , tmp] > 0) , ]
  gwa.pheno<-gwa.pheno[order(gwa.pheno[,"cistrans"]), ]
  gwa.pheno$"cistrans" = factor(gwa.pheno$"cistrans")
  
  gwa.boxplot[[f]] <- ggplot(gwa.pheno, aes(x=cistrans, y=beta)) +
    geom_boxplot(fill = phenocol[f, 2], colour = "black") +
    labs(title = tmptrait, x = NULL, y = NULL, fill = "Tipologia")
  
  tmp <- paste0("pcc.", tmptrait)
  pcc.pheno <- eqtl[which(eqtl[ , tmp] > 0) , ]
  pcc.pheno<-pcc.pheno[order(pcc.pheno[,"cistrans"]), ]
  pcc.pheno$"cistrans" = factor(pcc.pheno$"cistrans")
  
  pcc.boxplot[[f]] <- ggplot(pcc.pheno, aes(x=cistrans, y=beta)) +
    geom_boxplot(fill = phenocol[f, 2], colour = "black") +
    labs(title = tmptrait, x = NULL, y = NULL, fill = "Tipologia")
}

library(patchwork)
plot.gwa.boxplot <- gwa.boxplot[[1]] | gwa.boxplot[[2]] | gwa.boxplot[[3]] | gwa.boxplot[[4]]
plot.gwa.boxplot
plot.pcc.boxplot <- pcc.boxplot[[1]] | pcc.boxplot[[2]] | pcc.boxplot[[3]] | pcc.boxplot[[4]]
plot.pcc.boxplot





############################  
# BARPLOT 
############################
##### GWA ####
gwa.cistrans <- list()
for (f in 1:length(traits)){
  tmptrait<-traits[f]
  print(tmptrait)
  
  tmp <- paste0("gwa.", tmptrait)
  gwa.pheno <- eqtl[which(eqtl[ , tmp] > 0) , ]
  gwa.pheno$"cistrans" = factor(gwa.pheno$"cistrans")
  
  gwa.pheno.cistrans <- plyr::count(gwa.pheno, "cistrans")
  gwa.pheno.cistrans$"phenotype" <- tmptrait
  gwa.pheno.cistrans <- gwa.pheno.cistrans %>%
    mutate(perc = gwa.pheno.cistrans$"freq" / sum(gwa.pheno.cistrans$"freq"))
  gwa.cistrans[[f]] <- gwa.pheno.cistrans
}

gwa.cistrans = dplyr::bind_rows(gwa.cistrans) # transform into a dataframe
gwa.cistrans$"cistrans" = factor(gwa.cistrans$"cistrans")
gwa.cistrans$"phenotype" = factor(gwa.cistrans$"phenotype")

# change colors 
library(RColorBrewer) #Create a custom color scale
myColors <- as.character(c("leaf.length"="#c2b21d", "DZ.size"="#26cf74", "LER"="#c52d2d", "LED"="#006322"))

gwa.bar <- ggplot(gwa.cistrans, aes(x=cistrans, y=perc, fill=phenotype, group = phenotype)) +
  geom_bar(position = "dodge", stat="identity")  + 
  scale_fill_manual(values = myColors) +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.1)) +
  labs(title = "GWA.eQTL localization", x = NULL, y = "frequency")
gwa.bar  



##### PCC ####
pcc.cistrans <- list()
for (f in 1:length(traits)){
  tmptrait<-traits[f]
  print(tmptrait)
  
  tmp <- paste0("pcc.", tmptrait)
  pcc.pheno <- eqtl[which(eqtl[ , tmp] > 0) , ]
  pcc.pheno$"cistrans" = factor(pcc.pheno$"cistrans")
  
  pcc.pheno.cistrans <- plyr::count(pcc.pheno, "cistrans")
  pcc.pheno.cistrans$"phenotype" <- tmptrait
  pcc.pheno.cistrans <- pcc.pheno.cistrans %>%
    mutate(perc = pcc.pheno.cistrans$"freq" / sum(pcc.pheno.cistrans$"freq"))
  pcc.cistrans[[f]] <- pcc.pheno.cistrans
}

pcc.cistrans = dplyr::bind_rows(pcc.cistrans) # transform into a dataframe
pcc.cistrans$"cistrans" = factor(pcc.cistrans$"cistrans")
pcc.cistrans$"phenotype" = factor(pcc.cistrans$"phenotype")

library(RColorBrewer) #Create a custom color scale
myColors <- as.character(c("leaf.length"="#c2b21d", "DZ.size"="#26cf74", "LER"="#c52d2d", "LED"="#006322"))

pcc.bar <- ggplot(pcc.cistrans, aes(x=cistrans, y=perc, fill=phenotype, group = phenotype)) +
  geom_bar(position = "dodge", stat="identity")  + 
  scale_fill_manual(values = myColors) +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.1)) +
  labs(title = "PCC.eQTL localization", x = NULL, y = "frequency")
pcc.bar  
plot.bar <- gwa.bar + pcc.bar + plot_layout(guides = 'collect')
ggsave("cistrans.plot.tiff", plot=plot.bar)


# MERGE GRAPHS
library(patchwork)
plot.gwa.boxplot / plot.pcc.boxplot


fig5 <- (plot.maf | ((gwa.bar | pcc.bar) + plot_layout(guides = 'collect'))) / (plot.gwa.violin / plot.pcc.violin)
fig5

#  theme(panel.background = element_rect(fill = "white"), 
#        panel.grid.major.y = element_line(color = "grey90"), panel.grid.minor.y = element_line(color = "grey90"), panel.grid.major.x = element_line(color = "grey90"), legend.position = "null") +

ggsave("Fig5.tiff", plot = fig5)



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


####################################
# GET THE GRAPH
####################################
setwd(maindir)

#load colors
load("/home/maize.vib/PUBBLICATION/colors.by.phenotype.Rdata")
#load("/Users/mara/work_analysis/colors.by.phenotype.Rdata")
phenocol<-phenocol[c(2,1,4,3),]
stopifnot(all(phenocol[,1]==traits))


gwa.freq.snps <- plyr::count(gwa.eqtl, "snps")
gwa.freq.genes <- plyr::count(gwa.eqtl, "gene")
gwa.freq.genes <- gwa.freq.genes[order(gwa.freq.genes$freq), ]

pcc.freq.snps <- plyr::count(pcc.eqtl, "snps")
pcc.freq.genes <- plyr::count(pcc.eqtl, "gene")
pcc.freq.genes <- pcc.freq.genes[order(pcc.freq.genes$freq), ]
hist(pcc.freq.genes$freq)


pdf("eQTL.by.transcript.pdf", width=8, height=6)
par(mfrow=c(1,2))
#make a plot for GWA-eQTL putting in relation SNPs and number of transcripts per snp
plot(gwa.gbyslist[[1]], log="x", type="l", xlab="eQTL (log scale)", 
     #plot(gwa.gbyslist[[1]], type="l", xlab="eQTL", 
     ylab="Transcripts", lwd=1.5, main="GWA-eQTL", col=phenocol[1,2], xlim=c(1,20), ylim=c(1,55))
for(i in 2:length(gwa.gbyslist)){
  lines(gwa.gbyslist[[i]], type="l", lwd=1.5,  col=phenocol[i,2])
}
legend("topright", legend=traits, col=phenocol[,2], lty=1, cex=0.9, bty = "n")


#make a plot for PCC-eQTL putting in relation SNPs and number of transcripts per snp
plot(rev(pcc.gbyslist[[1]]), log="x", type="l", main="PCC-eQTL", xlab="eQTL",
     #plot(rev(pcc.gbyslist[[1]]), type="l", main="PCC-eQTL", xlab="eQTL",
     ylab="Transcripts", lwd=1.5, col=phenocol[1,2])
for(i in 2:length(pcc.gbyslist)){
  lines(rev(pcc.gbyslist[[i]]), type="l", lwd=1.5,  col=phenocol[i,2])
}
legend("topright", legend=traits, col=phenocol[,2], lty=1, cex=0.9, bty = "n")
dev.off()

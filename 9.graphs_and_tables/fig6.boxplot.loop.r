# MAKE BOXPLOT OF THE MOST IMPORTANT GENES
# # PLOT EXPRESSION and PHENOTYPE - SNP of INTEREST
####################################################
# performing a linear regression of sample genotypes on sample gene expression levels
R

library(data.table)
library(plyr)
library(tidyr)
library(ggplot2)

# variables
workdir = "/home/maize.vib/eqtl_results_2/eqtl_graphs"
SNP_file_name = "/home/maize.vib/genotypic_data_2/genotypes.chr.all.filtered.biallelic.txt"
snps_location_file_name = "/home/maize.vib/genotypic_data_2/snp.pos.filtered.biallelic.txt"
expression_file_name = "/home/maize.vib/expression_data/all.tpm.chr.0.05var.txt"
gene_location_file_name = "/home/maize.vib/annotation/gene_positions_V4.txt"
pheno_file_name = "/home/maize.vib/PUBBLICATION/7.RNA.pheno.cor/allpheno.txt"
interest_traits = c("leaf.length", "LER", "LED", "DZ.size") # then remove some traits we do not need (we can do with phenotypes of interests)


biparental <- c("RIL100","RIL103","RIL104","RIL105","RIL106","RIL107","RIL109","RIL11","RIL111","RIL113","RIL114","RIL115","RIL116","RIL117","RIL118","RIL119","RIL120","RIL121","RIL122","RIL123","RIL124","RIL125",
                "RIL126","RIL127","RIL128","RIL13","RIL130","RIL131","RIL132","RIL133","RIL134","RIL135","RIL136","RIL137","RIL138","RIL139","RIL140","RIL141","RIL16","RIL17","RIL18","RIL19","RIL23","RIL25","RIL26","RIL27",
                "RIL29","RIL30","RIL32","RIL33","RIL36","RIL37","RIL38","RIL39","RIL40","RIL42","RIL43","RIL44","RIL45","RIL46","RIL47","RIL48","RIL49","RIL50","RIL51","RIL52","RIL53","RIL54","RIL55","RIL57","RIL58","RIL59","RIL60","RIL61","RIL62","RIL64","RIL65","RIL66","RIL67","RIL69","RIL70","RIL71","RIL72","RIL73","RIL74","RIL75","RIL76","RIL80","RIL81","RIL82","RIL83","RIL84","RIL85","RIL88","RIL89","RIL90","RIL91","RIL92","RIL93","RIL94","RIL95","RIL96","RIL99")

magic <- c("RIL_8W_1", "RIL_8W_10", "RIL_8W_11","RIL_8W_12","RIL_8W_13","RIL_8W_14","RIL_8W_15","RIL_8W_16","RIL_8W_18","RIL_8W_19","RIL_8W_2","RIL_8W_20","RIL_8W_22","RIL_8W_23","RIL_8W_24","RIL_8W_26","RIL_8W_27","RIL_8W_28","RIL_8W_29","RIL_8W_3","RIL_8W_30","RIL_8W_31","RIL_8W_32","RIL_8W_33","RIL_8W_34","RIL_8W_35","RIL_8W_36","RIL_8W_37","RIL_8W_38","RIL_8W_39","RIL_8W_4","RIL_8W_40","RIL_8W_41","RIL_8W_42","RIL_8W_44","RIL_8W_45","RIL_8W_46",
           "RIL_8W_47","RIL_8W_48","RIL_8W_49","RIL_8W_5","RIL_8W_50","RIL_8W_51","RIL_8W_52","RIL_8W_53","RIL_8W_54","RIL_8W_55","RIL_8W_57","RIL_8W_58","RIL_8W_59","RIL_8W_6","RIL_8W_60","RIL_8W_61","RIL_8W_62","RIL_8W_63","RIL_8W_64","RIL_8W_65","RIL_8W_66","RIL_8W_67","RIL_8W_68","RIL_8W_69","RIL_8W_7","RIL_8W_70","RIL_8W_71","RIL_8W_72","RIL_8W_73","RIL_8W_74","RIL_8W_75","RIL_8W_76","RIL_8W_77",
           "RIL_8W_78","RIL_8W_79","RIL_8W_8","RIL_8W_80","RIL_8W_81","RIL_8W_82","RIL_8W_83","RIL_8W_84","RIL_8W_85","RIL_8W_86","RIL_8W_87","RIL_8W_88","RIL_8W_89","RIL_8W_9","RIL_8W_90","RIL_8W_91","RIL_8W_92","RIL_8W_93","RIL_8W_94","RIL_8W_95","RIL_8W_96","RIL_8W_97","RIL_8W_98","RIL_8W_99")

parents <- c("B73", "H99", "A632", "CML91", "F7", "HP301", "Mo17","W153R")


# LOAD INPUT DATA (and clean)
#############################
pheno = fread(pheno_file_name, data.table=F)
names(pheno) <- gsub(" ", ".", names(pheno))
names(pheno) = gsub("code", "samples", names(pheno))
pheno <- pheno[order(pheno$"samples"), ]
pheno <- pheno[which(pheno$"samples" %in% c(biparental, magic)), ]
samples <- pheno$"samples"


gt = fread(SNP_file_name, data.table=F)
names(gt) = gsub("NG-6280_", "", names(gt))
names(gt) = gsub("NG-6564_", "", names(gt))
names(gt) = gsub("NG-6860_", "", names(gt))
names(gt) = gsub("NG-7084_", "", names(gt))
names(gt) = gsub("P1_", "", names(gt))
names(gt) = gsub("R_", "RIL_8W_", names(gt))
new_names <- sapply(strsplit(names(gt), "_lib"), "[[", 1) # delete "_lib*", so the names correspond with pheno file
names(gt) = new_names
gt <- gt[ , c("SNP", biparental, magic)]


expr = fread(expression_file_name, data.table=F)
names(expr) = gsub("NG-6280_", "", names(expr))
names(expr) = gsub("NG-6564_", "", names(expr))
names(expr) = gsub("NG-6860_", "", names(expr))
names(expr) = gsub("NG-7084_", "", names(expr))
names(expr) = gsub("P1_", "", names(expr))
names(expr) = gsub("R_", "RIL_8W_", names(expr))
new_names <- sapply(strsplit(names(expr), "_lib"), "[[", 1) # delete "_lib*", so the names correspond with pheno file
names(expr) = new_names
expr <- expr[ , c("gene_id", biparental, magic)] # already ordered


library(RColorBrewer) #Create a custom color scale
myColors <- as.character(c("leaf.length"="#c2b21d", "DZ.size"="#26cf74", "LER"="#c52d2d", "LED"="#006322"))

load("/Users/mara/work_analysis/colors.by.phenotype.Rdata")
traits <- c("DZ.size", "leaf.length", "LED", "LER")
phenocol<-phenocol[c(2,1,4,3),]
phenocol <- as.data.frame(phenocol, stringsAsFactors = FALSE)
stopifnot(all(phenocol[,1]==traits))


######################################################
# PLOT EXPRESSION and PHENOTYPE - SNP of INTEREST
######################################################
snp_to_plot = "4_238222947"
snp.gt <- as.data.frame(gt[which(gt$SNP == snp_to_plot), ])
snp.gt <- t(snp.gt[, -1])
colnames(snp.gt)[[1]] <- snp_to_plot
snp.gt <- snp.gt[order(row.names(snp.gt)), ]
genoLong<-as.data.frame(snp.gt)
genoLong$"snp" = snp_to_plot
colnames(genoLong)[[1]] <- "genotype"
genoLong <- genoLong[order(row.names(genoLong)), ]


#################################
gene_list <- c("Zm00001d053442", "Zm00001d053633")
genes_plotted <-  list()

for (g in 1:length(gene_list)) {
  gene_to_plot = gene_list[[g]]
  gene.expr <- as.data.frame(expr[which(expr$gene_id == gene_to_plot), ])
  gene.expr <- t(gene.expr[, -1])
  colnames(gene.expr) <- gene_to_plot
  gene.expr <- gene.expr[order(rownames(gene.expr)), ]
  exprLong<-as.data.frame(gene.expr)
  exprLong$"gene" = gene_to_plot
  colnames(exprLong)[[1]] <- "expression"
  exprLong <- exprLong[order(rownames(exprLong)), ]
  
  dataLong = cbind(genoLong[,c("snp", "genotype")], exprLong[,c("gene", "expression")])
  #dataLong$"comparison" = paste(dataLong$snp, "vs", dataLong$gene)
  dataLong$"genotype" = factor(dataLong$"genotype")
  dataLong <- dataLong[!is.na(dataLong$genotype), ]
  dataLong$"population" <- NA
  dataLong[which(rownames(dataLong) %in% biparental), ]$"population" <- "bip"
  dataLong[which(rownames(dataLong) %in% magic), ]$"population" <- "magic"
  
  pl = ggplot(dataLong, aes(x=genotype, y=expression)) +
    geom_jitter(colour= "blue", position=position_jitter(width=0.25), na.rm = TRUE) +
    geom_boxplot(outlier.size=0, alpha=0.6, fill="grey") +
    facet_wrap(~gene) + theme_bw()
  genes_plotted[[g]] <- pl
} #for g

genes_plotted[[1]] | genes_plotted[[2]]


### phenotype
##################################
pheno_list <- c("LED", "leaf.length")
pheno_plotted <- list()

for (ph in 1:length(pheno_list)) {
  pheno_to_plot = pheno_list[[ph]]
  pheno_color <-  phenocol[which(phenocol$phenotype == pheno_to_plot),2]
  ph.val <- as.data.frame(pheno[ , pheno_to_plot])
  rownames(ph.val) = pheno$samples
  phenoLong<-ph.val
  phenoLong$"phenotype" = pheno_to_plot
  colnames(phenoLong)[[1]] <- "phenotype.value"
  phenoLong <- phenoLong[order(rownames(phenoLong)), ]
  
  #genoLong <- genoLong[-which(rownames(genoLong) == "RIL63"), ]
  
  #genoLong = tidyr::gather(snp.gt, snp, genotype, snp_to_plot)
  #exprLong = tidyr::gather(expr, gene, expression, gene_to_plot)
  dataLong = cbind(genoLong[,c("snp", "genotype")], phenoLong[,c("phenotype", "phenotype.value")], by = rownames(phenoLong))
  dataLong$"genotype" = factor(dataLong$"genotype")
  dataLong <- dataLong[!is.na(dataLong$genotype), ]
  dataLong$"population" <- NA
  dataLong[which(rownames(dataLong) %in% biparental), ]$"population" <- "bip"
  dataLong[which(rownames(dataLong) %in% magic), ]$"population" <- "magic"
  
  pl = ggplot(dataLong, aes(genotype, phenotype.value)) +
    geom_jitter(colour=pheno_color, position=position_jitter(width=0.25), na.rm = TRUE) +
    geom_boxplot(outlier.size=0, alpha=0.6, fill="grey", na.rm = TRUE) +
    facet_wrap(~phenotype) + theme_bw()
  pheno_plotted[[ph]] <- pl

} #for ph
  
pheno_plotted[[1]] | pheno_plotted[[2]]


# PUT ALL TOGETHER
########################################
library(patchwork)

p.tot = ((genes_plotted[[1]] | genes_plotted[[2]]) / (pheno_plotted[[1]] | pheno_plotted[[2]])) + 
  plot_annotation(title = "SNP in Chr 4 23,822,2947")
p.tot

ggsave(paste0(workdir, "/boxplot.4_238222947.tiff"), plot = p.tot)


##################################################################################################






# TEMPLATE TO PLOT WHATEVER GENE AND WHATEVER SNP:
######################################################
# PLOT EXPRESSION and PHENOTYPE - SNP of INTEREST
######################################################
snp_to_plot = "8_89590962"
snp.gt <- as.data.frame(gt[which(gt$SNP == snp_to_plot), ])
snp.gt <- t(snp.gt[, -1])
colnames(snp.gt)[[1]] <- snp_to_plot
snp.gt <- snp.gt[order(row.names(snp.gt)), ]
genoLong<-as.data.frame(snp.gt)
genoLong$"snp" = snp_to_plot
colnames(genoLong)[[1]] <- "genotype"
genoLong <- genoLong[order(row.names(genoLong)), ]


#################################
gene_list <- c("Zm00001d008742")
genes_plotted <-  list()

for (g in 1:length(gene_list)) {
  gene_to_plot = gene_list[[g]]
  gene.expr <- as.data.frame(expr[which(expr$gene_id == gene_to_plot), ])
  gene.expr <- t(gene.expr[, -1])
  colnames(gene.expr) <- gene_to_plot
  gene.expr <- gene.expr[order(rownames(gene.expr)), ]
  exprLong<-as.data.frame(gene.expr)
  exprLong$"gene" = gene_to_plot
  colnames(exprLong)[[1]] <- "expression"
  exprLong <- exprLong[order(rownames(exprLong)), ]
  
  dataLong = cbind(genoLong[,c("snp", "genotype")], exprLong[,c("gene", "expression")])
  #dataLong$"comparison" = paste(dataLong$snp, "vs", dataLong$gene)
  dataLong$"genotype" = factor(dataLong$"genotype")
  dataLong <- dataLong[!is.na(dataLong$genotype), ]
  dataLong$"population" <- NA
  dataLong[which(rownames(dataLong) %in% biparental), ]$"population" <- "bip"
  dataLong[which(rownames(dataLong) %in% magic), ]$"population" <- "magic"
  
  pl = ggplot(dataLong, aes(x=genotype, y=expression)) +
    geom_jitter(colour= "blue", position=position_jitter(width=0.25), na.rm = TRUE) +
    geom_boxplot(outlier.size=0, alpha=0.6, fill="grey") +
    facet_wrap(~gene) + theme_bw()
  genes_plotted[[g]] <- pl
} #for g

genes_plotted[[1]]


### phenotype
##################################
pheno_list <- c("leaf.length")
pheno_plotted <- list()

for (ph in 1:length(pheno_list)) {
  pheno_to_plot = pheno_list[[ph]]
  pheno_color <-  phenocol[which(phenocol$phenotype == pheno_to_plot),2]
  ph.val <- as.data.frame(pheno[ , pheno_to_plot])
  rownames(ph.val) = pheno$samples
  phenoLong<-ph.val
  phenoLong$"phenotype" = pheno_to_plot
  colnames(phenoLong)[[1]] <- "phenotype.value"
  phenoLong <- phenoLong[order(rownames(phenoLong)), ]
  
  #genoLong <- genoLong[-which(rownames(genoLong) == "RIL63"), ]
  
  #genoLong = tidyr::gather(snp.gt, snp, genotype, snp_to_plot)
  #exprLong = tidyr::gather(expr, gene, expression, gene_to_plot)
  dataLong = cbind(genoLong[,c("snp", "genotype")], phenoLong[,c("phenotype", "phenotype.value")], by = rownames(phenoLong))
  dataLong$"genotype" = factor(dataLong$"genotype")
  dataLong <- dataLong[!is.na(dataLong$genotype), ]
  dataLong$"population" <- NA
  dataLong[which(rownames(dataLong) %in% biparental), ]$"population" <- "bip"
  dataLong[which(rownames(dataLong) %in% magic), ]$"population" <- "magic"
  
  pl = ggplot(dataLong, aes(genotype, phenotype.value)) +
    geom_jitter(colour=pheno_color, position=position_jitter(width=0.25), na.rm = TRUE) +
    geom_boxplot(outlier.size=0, alpha=0.6, fill="grey", na.rm = TRUE) +
    facet_wrap(~phenotype) + theme_bw()
  pheno_plotted[[ph]] <- pl
  
} #for ph

pl <- pheno_plotted[[1]] | genes_plotted[[1]]
ggsave(file="Fig.S14.tiff", plot=pl)

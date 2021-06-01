#######################################
# MAKE CORRELATION PLOTS
# FIG-2 in the article
#######################################
#check correlation among traits, pops separated
#library(corrplot)
library(ggupset)
library(ggplot2)
library(data.table)
library(reshape2)


# variables son server
workdir<-"/home/maize.vib/eqtl_results_2/eqtl_graphs"


# fake variables on mara's mac
workdir = "/Users/mara/work_analysis"
SNP_file_name = "/Users/mara/work_analysis/maize.vib/inputs/genotypes.chr.all.filtered.biallelic.txt"
expression_file_name = "/Users/mara/work_analysis/maize.vib/inputs/all.tpm.chr.0.05var.txt"
pheno_file_name = "/Users/mara/work_analysis/maize.vib/inputs/allpheno.txt" #"/home/maize.vib/phenotypes/allpheno.txt"


biparental <- c("RIL100","RIL103","RIL104","RIL105","RIL106","RIL107","RIL109","RIL11","RIL111","RIL113","RIL114","RIL115","RIL116","RIL117","RIL118","RIL119","RIL120","RIL121","RIL122","RIL123","RIL124","RIL125",
                "RIL126","RIL127","RIL128","RIL13","RIL130","RIL131","RIL132","RIL133","RIL134","RIL135","RIL136","RIL137","RIL138","RIL139","RIL140","RIL141","RIL16","RIL17","RIL18","RIL19","RIL23","RIL25","RIL26","RIL27",
                "RIL29","RIL30","RIL32","RIL33","RIL36","RIL37","RIL38","RIL39","RIL40","RIL42","RIL43","RIL44","RIL45","RIL46","RIL47","RIL48","RIL49","RIL50","RIL51","RIL52","RIL53","RIL54","RIL55","RIL57","RIL58","RIL59","RIL60","RIL61","RIL62","RIL64","RIL65","RIL66","RIL67","RIL69","RIL70","RIL71","RIL72","RIL73","RIL74","RIL75","RIL76","RIL80","RIL81","RIL82","RIL83","RIL84","RIL85","RIL88","RIL89","RIL90","RIL91","RIL92","RIL93","RIL94","RIL95","RIL96","RIL99")
magic <- c("RIL_8W_1", "RIL_8W_10", "RIL_8W_11","RIL_8W_12","RIL_8W_13","RIL_8W_14","RIL_8W_15","RIL_8W_16","RIL_8W_18","RIL_8W_19","RIL_8W_2","RIL_8W_20","RIL_8W_22","RIL_8W_23","RIL_8W_24","RIL_8W_26","RIL_8W_27","RIL_8W_28","RIL_8W_29","RIL_8W_3","RIL_8W_30","RIL_8W_31","RIL_8W_32","RIL_8W_33","RIL_8W_34","RIL_8W_35","RIL_8W_36","RIL_8W_37","RIL_8W_38","RIL_8W_39","RIL_8W_4","RIL_8W_40","RIL_8W_41","RIL_8W_42","RIL_8W_44","RIL_8W_45","RIL_8W_46",
           "RIL_8W_47","RIL_8W_48","RIL_8W_49","RIL_8W_5","RIL_8W_50","RIL_8W_51","RIL_8W_52","RIL_8W_53","RIL_8W_54","RIL_8W_55","RIL_8W_57","RIL_8W_58","RIL_8W_59","RIL_8W_6","RIL_8W_60","RIL_8W_61","RIL_8W_62","RIL_8W_63","RIL_8W_64","RIL_8W_65","RIL_8W_66","RIL_8W_67","RIL_8W_68","RIL_8W_69","RIL_8W_7","RIL_8W_70","RIL_8W_71","RIL_8W_72","RIL_8W_73","RIL_8W_74","RIL_8W_75","RIL_8W_76","RIL_8W_77",
           "RIL_8W_78","RIL_8W_79","RIL_8W_8","RIL_8W_80","RIL_8W_81","RIL_8W_82","RIL_8W_83","RIL_8W_84","RIL_8W_85","RIL_8W_86","RIL_8W_87","RIL_8W_88","RIL_8W_89","RIL_8W_9","RIL_8W_90","RIL_8W_91","RIL_8W_92","RIL_8W_93","RIL_8W_94","RIL_8W_95","RIL_8W_96","RIL_8W_97","RIL_8W_98","RIL_8W_99")
parents <- c("B73", "H99", "A632", "CML91", "F7", "HP301", "Mo17","W153R")

traits<-c("leaf.length", "DZ.size", "LER", "LED")


setwd(workdir)


pheno<-fread(pheno_file_name, data.table=F)
rownames(pheno)<-pheno[,1]
pheno<-pheno[,-1]

#divide by population
bh<-pheno[grepl("BxH|AVERAGE", pheno$"family"),]
magic<-pheno[grep("MAGIC", pheno$"family"),]


# function to Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

# BIPARENTAL
####################################
#subset to interesting traits BIPARENTAL
bh.pheno<-bh[,c("leaf.length", "DZ.size", "LER", "LED")]

bh.cor<-cor(bh.pheno, use="complete.obs")
bh.cor <- get_upper_tri(bh.cor)
bh.cor <- melt(bh.cor, na.rm = T)
bh.cor$value <- as.numeric(format(round(bh.cor$value, 2), nsmall = 2))

cor.plot.bh <- ggplot(data = bh.cor, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+
  labs(y = NULL, x = NULL) + coord_fixed() + geom_text(aes(label = value))
cor.plot.bh



# MAGIC
####################################
#subset to interesting traits BIPARENTAL
magic.pheno<-magic[,c("leaf.length", "DZ.size", "LER", "LED")]

magic.cor<-cor(magic.pheno, use="complete.obs")
magic.cor <- get_upper_tri(magic.cor)
magic.cor <- melt(magic.cor, na.rm = T)
magic.cor$value <- as.numeric(format(round(magic.cor$value, 2), nsmall = 2))

cor.plot.magic <- ggplot(data = magic.cor, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+
  labs(y = NULL, x = NULL) + coord_fixed() + geom_text(aes(label = value))
cor.plot.magic



# SIMILAR TO VENN DIAGRAM PLOT
#################################
load("/Users/mara/work_analysis/colors.by.phenotype.Rdata")
phenocol<-phenocol[c(2,1,4,3),]
stopifnot(all(phenocol[,1]==traits))

library(UpSetR)
load("/Users/mara/work_analysis/maize.vib/overlapping.genes.Rdata")
listInput <- overlap # all overlapping genes each phenotype
bar.color=c(leaf.length = "#c2b21d", DZ.size = "#26cf74", LER = "#c52d2d", LED = "#006322")

pdf(NULL)
dev.control(displaylist="enable")
upset(fromList(listInput), order.by = "freq", 
      sets.bar.color=bar.color);
int.plot <- recordPlot()
invisible(dev.off())

# Display the saved plot
#grid::grid.newpage()
tiff("Fig.2C.venn.tiff")
int.plot
dev.off()



# CREATE FIG.2 FOR THE PAPER
###########################
library(gridExtra)
library(grid)
library(lattice)
library(gtable)
library(egg)
library(RColorBrewer)
library(ggplotify)
library(cowplot)

library(corrplot)
tiff("Fig.2A.biparental.tiff")
correlation_matrix.bh <-cor(bh.pheno, use="complete.obs")
cor.plot.bh <- corrplot.mixed(correlation_matrix.bh, upper = "circle", tl.col = "black", lower.col = "black")
dev.off()

tiff("Fig.2B.magic.tiff")
correlation_matrix.magic <-cor(magic.pheno, use="complete.obs")
cor.plot.magic <- corrplot.mixed(correlation_matrix.magic, upper = "circle", tl.col = "black", lower.col = "black")
dev.off()




# One figure in row 1 and two figures in row 2
par(mfcol=c(2,2))
corrplot.mixed(correlation_matrix.bh, upper = "circle", tl.col = "black", lower.col = "black")
corrplot.mixed(correlation_matrix.magic, upper = "circle", tl.col = "black", lower.col = "black")
upset(fromList(listInput), order.by = "freq", sets.bar.color=bar.color)

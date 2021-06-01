################################################################################################################
# MAKE histograms FROM eQTL RESULTS
# figure SUPPLEMENTARY 7/8/
################################################################################################################
R
library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)



workdir<-"/home/maize.vib/eqtl_results_2/eqtl_graphs"

# input files
eqtl_df = "/home/maize.vib/eqtl_results_2/eqtl.all.red.2.maf.txt" # dataframe eqtl reduced by LD, with all statistics
eqtl_file_name = "/home/maize.vib/eqtl_results_2/me.all.pca10.all.reduced.3.Rdata"
#gene_location_file_name = "/home/maize.vib/annotation/gene_positions_V4.txt"
gene_location_rdata = "/home/maize.vib/annotation/gene_positions_V4.Rdata" #the same of the .txt, object called genepos

load(gene_location_rdata) #object in bed file --> genepos


# fake variables
workdir = "/Users/mara/work_analysis"
SNP_file_name = "/Users/mara/work_analysis/maize.vib/inputs/genotypes.chr.all.filtered.biallelic.txt"
expression_file_name = "/Users/mara/work_analysis/maize.vib/inputs/all.tpm.chr.0.05var.txt"
pheno_file_name = "/Users/mara/work_analysis/maize.vib/inputs/allpheno.txt" #"/home/maize.vib/phenotypes/allpheno.txt"
eqtl_file_name = "/Users/mara/work_analysis/me.all.pca10.all.reduced.3.Rdata"
gene_position_file = "/Users/mara/work_analysis/gene_position.txt"


# fake variables 2
workdir = "C:/Users/Mara/OneDrive/eQTL_maize"
eqtl_file_name = "C:/Users/Mara/OneDrive/eQTL_maize/data/me.all.pca10.all.reduced.3.Rdata"
gene_position_file = "C:/Users/Mara/OneDrive/eQTL_maize/data/gene_position.txt"


traits<-c("leaf.length", "DZ.size", "LER", "LED")
setwd(workdir)

load(eqtl_file_name) # objects: me, eqtl_maf, pcc.eqtl, gwa.eqtl (Bonferroni threshold)

# change names trans cis 
eqtl$cistrans<-str_replace(eqtl$cistrans, "trans_chr", "distant")
eqtl$cistrans<-str_replace(eqtl$cistrans, "trans", "distant")
eqtl$cistrans<-str_replace(eqtl$cistrans, "cis", "local")

gwa.eqtl$cistrans<-str_replace(gwa.eqtl$cistrans, "trans_chr", "distant")
gwa.eqtl$cistrans<-str_replace(gwa.eqtl$cistrans, "trans", "distant")
gwa.eqtl$cistrans<-str_replace(gwa.eqtl$cistrans, "cis", "local")

pcc.eqtl$cistrans<-str_replace(pcc.eqtl$cistrans, "trans_chr", "distant")
pcc.eqtl$cistrans<-str_replace(pcc.eqtl$cistrans, "trans", "distant")
pcc.eqtl$cistrans<-str_replace(pcc.eqtl$cistrans, "cis", "local")




###################
# FIGURE SUPPL 7
###################
p<-eqtl[, c("gene", "snps", "cistrans")]

both <-p %>% 
  group_by(gene,cistrans) %>% summarise(count=n()) %>% 
  group_by(gene) %>% summarise(count=n()) %>% filter(count >1 ) # select genes that have both distant and local regulation

p[which(p$gene %in% both$gene), ]$cistrans = "both"

fig<-p %>% group_by(gene, cistrans) %>% summarise(count=n())

#fix naming
colnames(fig)<-c("Gene ID", "Class", "Count")
fig$Class <- factor(fig$Class, levels=c("local", "distant", "both"), labels=c("Local", "Distant", "Both"))


#ggplot(fig, aes(x=count, fill=cistrans))+geom_bar()+ 
#  scale_y_continuous(trans = "log10")+theme_bw()

g1<-ggplot(fig, aes(y=Count, x= Class, fill=Class))+geom_bar(stat="identity")+
  theme_bw() + ggtitle("Number of eQTL [genes]")
# +ggtitle("b")
g1


g2<-ggplot(fig, aes(y=Count, x=Class, fill=Class))+geom_boxplot()+ 
  #scale_y_continuous(trans = "log10")+
  theme_bw() + ggtitle("Number of SNPs controlling eQTL genes")
g2


fig.s7 <- g1 + (g2 + theme(legend.position = "none")) + plot_annotation(tag_levels = c('a','b')) + plot_layout(guides="collect")
ggsave(fig.s7, file="Fig.S7.pdf", width = 8, height=6)





############################
# SUPPLEMENTARY FIGURE 8 - B
# BOXPLOT 
############################
setwd(maindir)
traits<-c("leaf.length", "DZ.size", "LER", "LED")

#load colors
load("C:/Users/Mara/OneDrive/eQTL_maize/data/colors.by.phenotype.Rdata")
phenocol<-phenocol[c(2,1,4,3),]
stopifnot(all(phenocol[,1]==traits))

gwa.boxplot <- list()
pcc.boxplot <- list()
s8.boxplot <- list()

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
  
  all.pheno <- rbind(gwa.pheno, pcc.pheno)
  s8.boxplot[[f]] <- ggplot(all.pheno, aes(x=cistrans, y=beta)) +
    geom_boxplot(fill = phenocol[f, 2], colour = "black") +
    labs(title = tmptrait, x = NULL, y = NULL, fill = "Tipologia")
  
}

library(patchwork)
plot.s8.boxplot <- s8.boxplot[[1]] | s8.boxplot[[2]] | s8.boxplot[[3]] | s8.boxplot[[4]]
plot.s8.boxplot





# FIGURE SUPPL 8
###################
LED <- eqtl %>% filter(gwa.LED > 0 | pcc.LED >0)
LER <- eqtl %>% filter(gwa.LER > 0 | pcc.LER >0) 
LL <- eqtl %>% filter(gwa.leaf.length > 0 | pcc.leaf.length >0) 
DZ <- eqtl %>% filter(gwa.DZ.size > 0 | pcc.DZ.size >0)

LED$trait <- "LED"
LER$trait <- "LER"
LL$trait <- "leaf.length"
DZ$trait <- "DZ.size"

cp <- bind_rows(list(LED, LER, LL, DZ))
cp <- cp[ , c(1,2,8,23,26)]


# calculate percentage
freq.trait <- cp %>%
  group_by(cistrans, trait) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

freq.eqtl <- eqtl %>%
  group_by(cistrans) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
freq.eqtl$trait <- "tot"


data <- rbind(freq.eqtl, freq.trait)






# da qui
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











# FIGURE SUPPLEMENTARY (not used)
#################################

# calculate percentage
freq.eqtl <- eqtl %>%
  group_by(cistrans) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
freq.eqtl$"type" <- "tot"

freq.gwa <- gwa.eqtl %>%
  group_by(cistrans) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
freq.gwa$"type" <- "gwa.eqtl"

freq.pcc <- pcc.eqtl %>%
  group_by(cistrans) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
freq.pcc$"type" <- "pcc.eqtl"


data <- rbind(freq.eqtl, freq.gwa, freq.pcc)


hist.freq <- ggplot(data) + geom_bar(aes(type, freq, fill=cistrans), stat="identity") + 
  labs(title = "frequency eQTL localization", x= NULL) + theme(legend.title = element_blank())
hist.freq
#ggsave(filename="Fig.S5.tiff", plot =hist.freq)



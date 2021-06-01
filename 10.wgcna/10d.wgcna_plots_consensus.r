########################################################################################################################
# NETWORK ANALYSIS - BOXPLOT
# PLOTTING GENE MODULE MEMBERSHIP and GENE TRAIT SIGNIFICANCE
########################################################################################################################
R

library(WGCNA)
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)


workdir = "/home/maize.vib/network_splitted_pop/consensus"
workdir = "C:/Users/Mara/OneDrive/eQTL_maize"
setwd(workdir)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

#pop <-  c("bi", "multi")
setLabels = c("bi", "multi")
phenotypes = c("leaf.length", "LER", "LED", "DZ.size")


# Load the data saved in the first part
load(file = "/home/maize.vib/network_splitted_pop/consensus/Consensus-dataInput.RData")
load(file = "/home/maize.vib/network_splitted_pop/consensus/Consensus-NetworkConstruction-auto.RData")
load(file = "/home/maize.vib/network_splitted_pop/consensus/consensus_relating_traits.RData")


# Load the data saved in the first part FAKE
load(file = "C:/Users/Mara/OneDrive/eQTL_maize/data/Consensus-dataInput.RData")
load(file = "C:/Users/Mara/OneDrive/eQTL_maize/data/Consensus-NetworkConstruction-auto.RData")
load(file = "C:/Users/Mara/OneDrive/eQTL_maize/data/consensus_relating_traits.RData")



# FIGURE 5A 
############################
library(reshape2)
consensusCor
ok.modules <- c("MEmagenta", "MEtan", "MEroyalblue", "MEdarkorange", "MEgreen", "MEcyan", "MEorange", "MEturquoise", "MEpurple", "MElightcyan", "MEyellowgreen")
cons <- consensusCor[which(rownames(consensusCor) %in% ok.modules), ]
cons <- melt(cons)

subtype_colors= rownames(cons) %>% str_replace("ME", "")

# options(digits=2)
fig5a <- ggplot(data = cons, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab") +
  geom_text(aes(Var2, Var1, label = signif(value, digits=2)), color = "black", size = 4) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1,  size = 12, hjust = 1))+
  coord_fixed() +
  ggtitle("Consensus module-trait relationship") +
  xlab("") + ylab("Modules") 
ggsave(filename="Fig.S5A.tiff", fig5a)

sizeGrWindow(10,7)
pdf(file = "graphs/ModuleTraitRelationships-consensus.pdf", wi = 10, he = 7);
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = consensusCor,
               xLabels = names(Traits[[set]]),
               yLabels = rownames(consensusCor),
               ySymbols = rownames(consensusCor),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste0("Consensus module-trait relationships across\n", setLabels[1], "parental and", setLabels[2], "parental"))
dev.off()
rm(textMatrix)








# PLOT THE CONNECTIVITY kME
##########################################
for (p in setLabels) {names(geneMM[[p]]) <- labels2colors(as.numeric(substring(names(geneMM[[p]]), 3)))} #change names in list
# in plentyofcolors there's what I need for violin plots
module_list <- final.modules$colorME
module_list_plus <- unique(c(module_list, names(table(strega$colorME))))
module_list_plus <- module_list_plus[!module_list_plus %in% "grey"] #remove grey color 

# by population on final.modules colors 
##########################################
modules_plotted_pop <- list()

for (m in 1:length(module_list_plus)) {
  module_to_plot = module_list_plus[[m]]
  
  tmp1.geneMM <- as.data.frame(geneMM[["bi"]][, module_to_plot])
  setnames(tmp1.geneMM, colnames(tmp1.geneMM[1]), "kME")
  tmp1.geneMM$"gene_id" <- rownames(geneMM[["bi"]])
  tmp1.geneMM$"pop" <- "bi"
  
  tmp2.geneMM <- as.data.frame(geneMM[["multi"]][, module_to_plot])
  setnames(tmp2.geneMM, colnames(tmp2.geneMM[1]), "kME")
  tmp2.geneMM$"gene_id" <- rownames(geneMM[["multi"]])
  tmp2.geneMM$"pop" <- "multi"
  
  dataLong <- rbind(tmp1.geneMM, tmp2.geneMM)
  dataLong$"pop" <- factor(dataLong$"pop")
  dataLong <- dataLong[!is.na(dataLong$"pop"), ]
  
  dataLong <- left_join(dataLong, plentyofcolors, by="gene_id")
  dataLong <- dataLong[which(dataLong$"colorME" == module_to_plot) , ]
  
  pl = ggplot(dataLong, aes(x=pop, y=kME)) +
    geom_violin(trim=FALSE, fill='#A4A4A4', color="black") +
    geom_boxplot(width=0.1, fill=module_to_plot) +
    theme_minimal()
  
  #stat_summary(fun=median, geom="point", size=2, color="black")
  
  modules_plotted_pop[[module_to_plot]] <- pl
} #for m

m_pl <- modules_plotted_pop

p.tot <- (m_pl[[1]] | m_pl[[2]] | m_pl[[3]] | m_pl[[4]] | m_pl[[5]] | m_pl[[6]])/
  (m_pl[[7]] | m_pl[[8]] | m_pl[[9]] | m_pl[[10]] | m_pl[[11]] | m_pl[[12]]) /
  (m_pl[[13]] | m_pl[[14]] | m_pl[[15]] | m_pl[[16]] | m_pl[[17]] | m_pl[[18]])
ggsave(paste0(workdir, "/violinplot_kME_pop.tiff"), plot = p.tot)
rm(m, module_to_plot, tmp1.geneMM, tmp2.geneMM, dataLong, pl, m_pl)


# on CONSENSUS geneMM on final.modules colors 
##########################################
#module_list <- final.modules$colorME
modules_plotted_consensus <- list()

# select genes int.eqtl and add the points on violin plot
g <- strega[ , c("gene_id", "colorME")]
g <- g %>% drop_na()
#setnames(c.geneModuleMembership, "gene", "gene_id")
g <- left_join(g, c.geneModuleMembership, by="gene_id")

for (m in 1:length(module_list_plus)) {
  module_to_plot = module_list_plus[[m]]
  
  tmp1.geneMM <- as.data.frame(c.geneModuleMembership[, paste0("MM",module_to_plot)])
  setnames(tmp1.geneMM, colnames(tmp1.geneMM[1]), "kME")
  tmp1.geneMM$"gene_id" <- rownames(c.geneModuleMembership)
  addthispoint <- g[which(g$colorME == module_to_plot), c("gene_id", "colorME", paste0("MM",module_to_plot))]
  addthispoint$"coloME" <- factor(addthispoint$"colorME")
  setnames(addthispoint, paste0("MM",module_to_plot), "kME")

  dataLong <- tmp1.geneMM[!is.na(tmp1.geneMM$"kME"), ]
  dataLong <- left_join(dataLong, plentyofcolors, by="gene_id")
  dataLong[which(dataLong$"colorME" != module_to_plot) , ]$"colorME" <- "outside"
  dataLong$"colorME" <- factor(dataLong$"colorME")
#  dataLong$"colorME" <- factor(dataLong$"colorME", levels=dataLong$"colorME"[order(c(module_to_plot, "outside"))])
#  dataLong <- left_join(dataLong, gTS, by="gene_id")

  pl = ggplot(dataLong, aes(x=reorder(colorME, kME, FUN = median), y=kME)) +
    geom_violin(trim=FALSE, fill='#A4A4A4', color="black") +
    geom_boxplot(width=0.1, fill=module_to_plot) +
    coord_cartesian(ylim = c(-0.8, 1)) + xlab(NULL) +
    geom_point(data = addthispoint, color = "blue") + 
    theme(axis.text.x=element_text(angle=45, hjust=1))
  modules_plotted_consensus[[module_to_plot]] <- pl
} #for m

m_pl <- modules_plotted_consensus

p.tot2 <- (m_pl[[1]] | m_pl[[2]] | m_pl[[3]] | m_pl[[4]] | m_pl[[5]] | m_pl[[6]])/
  (m_pl[[7]] | m_pl[[8]] | m_pl[[9]] | m_pl[[10]] | m_pl[[11]] | m_pl[[12]]) /
  (m_pl[[13]] | m_pl[[14]] | m_pl[[15]] | m_pl[[16]] | m_pl[[17]] | m_pl[[18]]) +
plot_annotation(title="Gene module membership")

ggsave(paste0(workdir, "/violinplot_kME_module_outside.tiff"), plot = p.tot2)
rm(m, module_to_plot, tmp1.geneMM, tmp2.geneMM, dataLong, pl, m_pl, int.g)




# SCATTERPLOT GENE TRAIT SIGNIFICANCE
##########################################
modules_plotted_geneTS <- list()

setnames(c.geneTraitSignificance, "gene", "gene_id")
tmp <- list()
pheno_plotted <- list()

# prepare geneTraitSignificance (consensus) in a format suitable for ggplot
for (ph in phenotypes){
  tmp.tmp <- c.geneTraitSignificance[, c(ph,"gene_id")]
  tmp.tmp$"pheno" <- ph
  setnames(tmp.tmp, ph, "geneTS")
  tmp[[ph]] <- tmp.tmp
}
geneTS <- dplyr::bind_rows(tmp)

for (m in 1:length(module_list)) {
  module_to_plot = module_list[[m]]

  #extract geneMM for alla genes, only in the module_to_plot of interest
  tmp1.geneMM <- as.data.frame(c.geneModuleMembership[, paste0("MM",module_to_plot)])
  setnames(tmp1.geneMM, colnames(tmp1.geneMM[1]), "kME")
  tmp1.geneMM$"gene_id" <- rownames(c.geneModuleMembership)
  
  # merge all datasets, ready for ggplot
  dataLong <- left_join(geneTS, plentyofcolors, by="gene_id")
  #dataLong[which(dataLong$"colorME" != module_to_plot) , ]$"colorME" <- "outside"
  dataLong <- dataLong[which(dataLong$"colorME" == module_to_plot) , ]
  dataLong <- left_join(dataLong, tmp1.geneMM, by="gene_id")
  dataLong$"pheno" <- factor(dataLong$"pheno")
  dataLong$"colorME" <- factor(dataLong$"colorME")
  
  pl = ggplot(dataLong, aes(x=geneTS, y=kME, color=pheno), na.rm=T) + geom_point(size=2, shape=23)

  pheno_plotted[[ph]] <- pl
} #for m




# SCATTERPLOT GENE TRAIT SIGNIFICANCE and kME
##############################################
modules_plotted_geneTS <- list()

setnames(c.geneTraitSignificance, "gene", "gene_id")
tmp <- list()
module_plotted <- list()

# prepare geneTraitSignificance (consensus) in a format suitable for ggplot
for (ph in phenotypes){
  tmp.tmp <- c.geneTraitSignificance[, c(paste0("gTS.", ph),"gene_id")]
  tmp.tmp$"pheno" <- ph
  setnames(tmp.tmp, paste0("gTS.", ph), "geneTS")
  tmp[[ph]] <- tmp.tmp
}
geneTS <- dplyr::bind_rows(tmp)

for (m in 1:length(module_list)) {
  module_to_plot = module_list[[m]]
  
  #extract geneMM for alla genes, only in the module_to_plot of interest
  tmp1.geneMM <- as.data.frame(c.geneMM.mean[, module_to_plot])
  setnames(tmp1.geneMM, colnames(tmp1.geneMM[1]), "kME")
  tmp1.geneMM$"gene_id" <- rownames(c.geneMM.mean)
  
  # merge all datasets, ready for ggplot
  dataLong <- left_join(geneTS, plentyofcolors, by="gene_id")
  #dataLong[which(dataLong$"colorME" != module_to_plot) , ]$"colorME" <- "outside"
  dataLong <- dataLong[which(dataLong$"colorME" == module_to_plot) , ]
  dataLong <- left_join(dataLong, tmp1.geneMM, by="gene_id")
  dataLong$"pheno" <- factor(dataLong$"pheno")
  dataLong$"colorME" <- factor(dataLong$"colorME")
  
  pl = ggplot(dataLong, aes(x=geneTS, y=kME, color=pheno), na.rm=T) + geom_point(size=2, shape=20)
  
  module_plotted[[module_to_plot]] <- pl
} #for m


save(module_plotted, file="module_plotted.rdata")




# plot GENE TRAIT SIGNIFICANCE --> NOT USED!!!
##########################################
#module_list <- final.modules$colorME
modules_plotted_geneTS <- list()

setnames(c.geneTraitSignificance, "gene", "gene_id")
tmp <- list()

for (m in 1:length(module_list)) {
  module_to_plot = module_list[[m]]
  for (ph in phenotypes){
    tmp.tmp <- c.geneTraitSignificance[, c(ph,"gene_id")]
    tmp.tmp$"pheno" <- ph
    setnames(tmp.tmp, ph, "geneTS")
    tmp[[ph]] <- tmp.tmp
  }
  dataLong <- dplyr::bind_rows(tmp)
  dataLong <- left_join(dataLong, plentyofcolors, by="gene_id")
  dataLong[which(dataLong$"colorME" != module_to_plot) , ]$"colorME" <- "outside"
  dataLong$"colorME" <- factor(dataLong$"colorME")
  dataLong$"pheno" <- factor(dataLong$"pheno")
  
  ggplot2.boxplot(data=df, xName='dose',yName='len', groupName='supp', 
                  position=position_dodge(0.8),
                  backgroundColor="white", groupColors=c('#999999','#E69F00'),
                  legendPosition="top")  
  
  pl = ggplot(dataLong, aes(pheno, geneTS)) +
    geom_boxplot(outlier.size=0, alpha=0.6, fill="grey", na.rm = TRUE) +
    facet_wrap(~colorME) + theme_bw()
  pheno_plotted[[ph]] <- pl
  
  
  pl = ggplot(dataLong, aes(x=colorME, y=kME)) +
    
    geom_violin(trim=FALSE, fill='#A4A4A4', color="black") +
    geom_boxplot(width=0.1, fill=module_to_plot) +
    coord_cartesian(ylim = c(-0.8, 1)) + 
    theme(axis.text.x=element_text(angle=45, hjust=1))
  modules_plotted_geneTS[[module_to_plot]] <- pl
} #for m

m_pl <- modules_plotted_geneTS

p.tot2 <- (m_pl[[1]] | m_pl[[2]] | m_pl[[3]] | m_pl[[4]] | m_pl[[5]] | m_pl[[6]]) /
  (m_pl[[7]] | m_pl[[8]] | m_pl[[9]] | m_pl[[10]] | m_pl[[11]] | m_pl[[12]])

ggsave(paste0(workdir, "/violinplot_kME_module_outside.tiff"), plot = p.tot2)
rm(m, module_to_plot, tmp1.geneMM, tmp2.geneMM, dataLong, pl, m_pl)



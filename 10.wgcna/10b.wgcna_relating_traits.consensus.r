########################################################################################################################
# NETWORK ANALYSIS - WGCNA - Consensus network analysis
# 4 - Relating CONSENSUS modules to phenotypic data
########################################################################################################################
R

library(WGCNA)
library(data.table)
library(dplyr)
library(tidyr)
library(reshape2)

workdir = "/home/maize.vib/network_splitted_pop/consensus"
#workdir = "/Users/mara/work_analysis/maize.vib/network_splitted_pop/"
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

exprSize = checkSets(multiExpr);
nSets = exprSize$nSets

Traits <-  PhenoValues
Pheno <- PhenoValues

# multiExpr and Traits containing the expression and trait data
# The consensus network analysis results are represented by the variables consMEs, moduleLabels, moduleColors, and consTree




# RELATE CONSENSUS MODULE TO EXTERNAL PHENOTYIPIC DATA 
########################################################################################################################
# module eigengenes used to relate consensus modules to phenotypic data
# link phenotypic data to consensus module eigengenes in each of the two populations
# for each gene there is a single module assignment, while there is a consensus module eigengenes for each population
# Similarly, there are phenotypic data separately for biparental and multiparental

# Set up variables to contain the module-trait correlations
moduleTraitCor = list();
moduleTraitPvalue = list();

# Calculate the correlations
for (set in 1:length(setLabels)) {
  moduleTraitCor[[set]] = cor(consMEs[[set]]$data, Traits[[set]], use = "p");
  moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], exprSize$nSamples[set]);
}
names(moduleTraitCor) <-  setLabels
names(moduleTraitPvalue) <-  setLabels
rm(set)

# display the module-trait relationships using a color-coded table. 
# print the correlations and the corresponding p-values (entris color-coded by p-value significance).

# Convert numerical lables to colors for labeling of modules in the plot
MEColors = labels2colors(as.numeric(substring(names(consMEs[[1]]$data), 3)));
MEColorNames = paste("ME", MEColors, sep="");



# PLOT the module-trait relationship table for each population 
########################################################################################################################
for (set in setLabels) {
  sizeGrWindow(10,7) # Open a suitably sized window (the user should change the window size if necessary)
  pdf(file = paste0("graphs/ModuleTraitRelationships-", set, "parental.pdf"), wi = 10, he = 7);
  textMatrix = paste(signif(moduleTraitCor[[set]], 2), "\t(", signif(moduleTraitPvalue[[set]], 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor[[set]])
  
  par(mar = c(6, 8.8, 3, 2.2));
  labeledHeatmap(Matrix = moduleTraitCor[[set]], # 39x4
                 xLabels = names(Traits[[set]]),
                 yLabels = MEColorNames, #39
                 ySymbols = MEColorNames, #39
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste0("Module-trait relationships in ", set, "parental"))
  dev.off();
} # for set



# SUMMARIZE MODULE-TRAIT RELATIONSHIP
########################################################################################################################
# conservative way of measure module-trait relationship: for each module-trait pair take the correlation that has the lower absolute value in the two sets if the two correlations have the same sign, and zero relationship if the two correlations have opposite signs.

# Initialize matrices to hold the consensus correlation and p-value
consensusCor = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]));
consensusPvalue = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]));

# --> take modules that are significative in both pops and that have the same sign/direction

# Find consensus negative correlations
negative = moduleTraitCor[[1]] < 0 & moduleTraitCor[[2]] < 0;
consensusCor[negative] = pmax(moduleTraitCor[[1]][negative], moduleTraitCor[[2]][negative]); # very conservative
consensusPvalue[negative] = pmax(moduleTraitPvalue[[1]][negative], moduleTraitPvalue[[2]][negative]);

# Find consensus positive correlations
positive = moduleTraitCor[[1]] > 0 & moduleTraitCor[[2]] > 0;
consensusCor[positive] = pmin(moduleTraitCor[[1]][positive], moduleTraitCor[[2]][positive]); # very conservative
consensusPvalue[positive] = pmax(moduleTraitPvalue[[1]][positive], moduleTraitPvalue[[2]][positive]);
rm(negative,positive)

# display the consensus module–trait relationships again using a color-coded table:
textMatrix = paste0(signif(consensusCor, 2), "\t(", signif(consensusPvalue, 1), ")");
#textMatrix <-  gsub("NA\t(NA)", "NA", textMatrix)
dim(textMatrix) = dim(moduleTraitCor[[set]])

rownames(consensusCor) <- paste0("ME", labels2colors(as.numeric(substring(rownames(moduleTraitCor[[1]]), 3))))
colnames(consensusCor) <- colnames(moduleTraitCor[[1]])
rownames(consensusPvalue) <- paste0("ME", labels2colors(as.numeric(substring(rownames(moduleTraitCor[[1]]), 3))))
colnames(consensusPvalue) <- colnames(moduleTraitCor[[1]])



sizeGrWindow(10,7)
pdf(file = "graphs/ModuleTraitRelationships-consensus.pdf", wi = 10, he = 7);
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = consensusCor,
               xLabels = names(Traits[[set]]),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste0("Consensus module-trait relationships across\n", setLabels[1], "parental and", setLabels[2], "parental"))
dev.off()
rm(textMatrix)

# The advantage of the consensus relationship table is that it isolates the module-trait relationships that are present in both sets, and hence may be in a sense considered validated



################################################################################################
# GENE SIGNIFICANCE and MODULE MEMBERSHIP - Gene relationship to phenotype and important modules
################################################################################################
# quantify associations of individual genes with the phenotype through Gene Significance (GS):
# GS --> (the absolute value of) is defined as the correlation between the gene and the trait. 
# For each module --> define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile. 
# This permits to quantify the similarity of all genes on the array to every module.

# Define variable for a Trait of interest containing the corresponding column of pheData (phenotypes matrix, pheno measure for each sample)



# GENE-TRAIT SIGNIFICANCE
#################################
geneTraitSignificance <- list() 

for (set in setLabels) {
  geneTS <- data.frame(row.names=names(multiExpr[[set]]$data), stringsAsFactors = F)
  
  for(t in phenotypes){
    datExpr <- multiExpr[[set]]$data
    
    pheno.tmp = as.data.frame(Pheno[[set]][ , t]); # measurement of the phenotype for each individual
    names(pheno.tmp) = t
    rownames(pheno.tmp) <- rownames(Pheno[[set]])
    
    geneTS[ , t] <-  cor(datExpr, pheno.tmp, use = "p")
  } # for t in traits
  setnames(geneTS, paste0("GS.", names(geneTS))) 
  geneTraitSignificance[[set]] <- geneTS
} # for set 
rm(set, geneTS, datExpr,pheno.tmp)


# gene-trait significance pvalue
GSPvalue <- list()

for (set in 1:length(setLabels)) {
  GSPvalue[[set]] <- data.frame(row.names=names(multiExpr[[set]]$data), stringsAsFactors = F)
  GSPvalue[[set]] <-  corPvalueStudent(as.matrix(geneTraitSignificance[[set]]), nSamples[set]);
  # setnames(GSPvalue[[set]], paste0("p.GS.", names(GSPvalue[[set]])))
}

# Initialize matrices to hold the consensus correlation and p-value
c.geneTraitSignificance = matrix(NA, nrow = ncol(multiExpr[[set]]$data), ncol = ncol(moduleTraitCor[[1]]))
c.geneTS.pvalue = matrix(NA, nrow = ncol(multiExpr[[set]]$data), ncol = ncol(moduleTraitCor[[1]]))

# Find consensus negative correlations
negative = geneTraitSignificance[[1]] < 0 & geneTraitSignificance[[2]] < 0;
c.geneTraitSignificance[negative] = pmax(geneTraitSignificance[[1]][negative], geneTraitSignificance[[2]][negative]); # very conservative
c.geneTS.pvalue[negative] = pmax(GSPvalue[[1]][negative], GSPvalue[[2]][negative]);

# Find consensus positive correlations
positive = geneTraitSignificance[[1]] > 0 & geneTraitSignificance[[2]] > 0;
c.geneTraitSignificance[positive] = pmax(geneTraitSignificance[[1]][positive], geneTraitSignificance[[2]][positive]); # very conservative
c.geneTS.pvalue[positive] = pmax(GSPvalue[[1]][positive], GSPvalue[[2]][positive]);

c.geneTraitSignificance <- as.data.frame(c.geneTraitSignificance)
rownames(c.geneTraitSignificance) <- colnames(multiExpr[[set]]$data)
names(c.geneTraitSignificance) <- colnames(moduleTraitCor[[1]])

c.geneTS.pvalue <- as.data.frame(c.geneTS.pvalue)
rownames(c.geneTS.pvalue) <- colnames(multiExpr[[set]]$data)
names(c.geneTS.pvalue) <- colnames(moduleTraitCor[[1]])

rm(positive, negative)




################################################################################################
# SELECTION of HUB GENES and INTERESTING GENES
# !!! HUB GENES NOT INCLUDED IN THE PAPER 
################################################################################################

# #chooseOneHubInEachModule 
# ################################################################################################
# net_type = "signed"
# softPower = 8

# colorh <- labels2colors(net$colors) # it is the object moduleColors

# hubs <- list()
# one.hubs <- list()
# for (set in setLabels) {
#   hubs[[set]]  = chooseTopHubInEachModule(multiExpr[[set]]$data, colorh,  power = softPower, type= net_type)
#   one.hubs[[set]]  = chooseOneHubInEachModule(multiExpr[[set]]$data, colorh,  power = softPower, type= net_type)
# }


# choose genes with high connectivity
# (intersection between bi and multi kTotal)
dieci.bi <- kk[[1]] %>% group_by(colorME) %>% filter(kTotal > quantile(kTotal, 0.9))
dieci.multi <- kk[[2]] %>% group_by(colorME) %>% filter(quantile(kTotal, 0.9)< kTotal)
dieci <- dieci.bi[which(dieci.bi$gene_id %in% dieci.multi$gene_id), "gene_id"]
dieci <- left_join(dieci, kk[[1]][,c("gene_id", "colorME")], by="gene_id")

dieci.dieci <- list()
c.geneTraitSignificance$"gene_id" <- rownames(c.geneTraitSignificance)
c.geneMM.mean$"gene_id" <- rownames(c.geneMM.mean)
c.geneTS.pvalue$"gene_id" <- rownames(c.geneTS.pvalue)
module_list_tot <- names(geneMM[[1]])
module_list_tot <- module_list_tot[!module_list_tot %in% c("gene_id", "gene", "grey")]
dieci.ph <- list()
for (ph in phenotypes){
  field <- paste0("gTS.", ph)
  tmp <- left_join(dieci, c.geneTraitSignificance[, c("gene_id", field)], by="gene_id")
  tmp <- tmp[which(abs(tmp[, field]) > 0.2), ] 
  tmp <- left_join(tmp, c.geneMM.mean, by="gene_id")
  # select genes with correlation ≥ 0.2 for the trait
  for (m in module_list_tot) {
    filt.mm <- tmp %>% drop_na %>% filter(colorME == m)
    filt.mm <- left_join(filt.mm, c.geneMM.mean[, c(m, "gene_id")], by="gene_id")
    
    # selection kME ≥ 0.8
    if (nrow(filt.mm) > 0) {
      filt.g <- filt.mm[which(abs(filt.mm[, m]) >= 0.8), ]$"gene_id"
      dieci.dieci[[m]] <- filt.g
    }
    dieci.ph[[ph]] <- dieci.dieci
  } # m
}
  

library(reshape2)
dieci <- reshape2::melt(dieci.ph) # these are the hub genes 
unique(dieci$value) %in% unique(highMM.mean$gene_id)



# choose INTERESTING genes
################################################################################################
# GS > 0.2 AND kME > 0.8 and Gene Sign < 0.05
phenotypes = c("leaf.length", "LER", "LED", "DZ.size")

filtered.genes <- list()
highMM.genes.ph <- list()
c.geneTraitSignificance$"gene_id" <- rownames(c.geneTraitSignificance)
c.geneModuleMembership$"gene_id" <- rownames(c.geneModuleMembership)
module_list_tot <- names(geneMM[[1]])
module_list_tot <- module_list_tot[!module_list_tot %in% c("gene_id", "gene", "grey")]

#c.geneTS.pvalue$"gene_id" <- rownames(c.geneTS.pvalue)

for (ph in phenotypes) {
  tmpg <- which(c.geneTS.pvalue[, ph] < 0.05) #gene significative for the trait in the consensus
  filt <- c.geneTraitSignificance[tmpg, c("gene_id", ph)]
  # select genes with correlation ≥ 0.2 for the trait
  filt <- filt[which(abs(filt[, ph]) > 0.2), ] 
  
  # select genes for each module with geneMM > 0.8
  for (m in module_list_tot) {
    filt.m <- left_join(filt, c.geneModuleMembership[, c(paste0("MM", m), "gene_id")], by="gene_id")
    # selection kME ≥ 0.8
    filt.m <- filt.m[which(abs(filt.m[, paste0("MM", m)]) >= 0.8), ]$"gene_id"
    filtered.genes[[m]] <- filt.m
  } # m
  highMM.genes.ph[[ph]] <- filtered.genes
} # set

library(reshape2)
highMM <- reshape2::melt(highMM.genes.ph) # these are the hub genes 
high.genes <- unique(highMM$value)
write.table(high.genes, file="high.genes.txt", row.names = F, col.names = F, quote = F)
setnames(highMM, "value", "gene_id")
highMM <- left_join(highMM, c.geneTraitSignificance[, c("gene_id",phenotypes)], by="gene_id")
for (ph in phenotypes) {setnames(highMM, ph, paste0("GS.",ph))}
setnames(highMM, "L2", "mod_color")
setnames(highMM, "L1", "phenotype")



# find interesting genes from the geneMM MEAN (to see what changes!)
phenotypes = c("leaf.length", "LER", "LED", "DZ.size")

filtered.genes <- list()
highMM.genes.ph <- list()
c.geneTraitSignificance$"gene_id" <- rownames(c.geneTraitSignificance)
c.geneMM.mean$"gene_id" <- rownames(c.geneMM.mean)
module_list_tot <- names(geneMM[[1]])
module_list_tot <- module_list_tot[!module_list_tot %in% c("gene_id", "gene", "grey")]
highMM.gTS <- list()
for (ph in phenotypes) {
  tmpg <- which(c.geneTS.pvalue[, ph] < 0.05) #gene significative for the trait in the consensus
  filt <- c.geneTraitSignificance[tmpg, c("gene_id", paste0("gTS.", ph))]
  # select genes with correlation ≥ 0.2 for the trait
  filt <- filt[which(abs(filt[, paste0("gTS.", ph)]) > 0.2), ] 
  
  # select genes for each module with geneMM > 0.8
  for (m in module_list_tot) {
    filt.mm <- left_join(filt, c.geneMM.mean[, c(m, "gene_id")], by="gene_id")
    # selection kME ≥ 0.8
    filt.g <- filt.mm[which(abs(filt.mm[, m]) >= 0.8), ]$"gene_id"
    filtered.genes[[m]] <- filt.g
  highMM.genes.ph[[ph]] <- filtered.genes
} # set

library(reshape2)
highMM.mean <- reshape2::melt(highMM.genes.ph) # these are the hub genes 
setnames(highMM.mean, "value", "gene_id")
highMM.mean <- left_join(highMM.mean, c.geneTraitSignificance[, c("gene_id",phenotypes)], by="gene_id")
for (ph in phenotypes) {setnames(highMM.mean, ph, paste0("GS.",ph))}
setnames(highMM.mean, "L2", "mod_color")
setnames(highMM.mean, "L1", "phenotype")
g <- unique(highMM.mean$gene_id)
write.table(g, file="high.genes.mean.txt", row.names = F, col.names = F, quote = F)



highMM.mean$gene_id %in% highMM$gene_id
highMM$gene_id %in% highMM.mean$gene_id

super.hub <- names(which(table(highMM.mean$gene_id) >1))
write.table(super.hub, file="superhub.txt", row.names = F, col.names = F, quote = F)




# look the enrichemnt 
# look if eqtl genes are here 
eqtl_file = "/home/maize.vib/eqtl_results_2/me.all.pca10.all.reduced.3.Rdata"
load(eqtl_file)
setnames(eqtl, "gene", "gene_id")
setnames(highMM, "value", "gene_id")



# PLOT
sizeGrWindow(8,6)
par(mfrow=c(2,2))
which.color="turquoise";
restrictGenes=colorh==which.color
verboseScatterplot(kIM[["bi"]]$kWithin[ restrictGenes],
                   (geneMM[["bi"]][restrictGenes, which.color])^softPower,
                   col=which.color,
                   xlab="Intramodular Connectivity",
                   ylab="(Module Membership)^8")



# TAKE only the SIGN module for each trait  and create a data.frame with the correlation values
pvalue = 0.05 #threshold
final.modules <- matrix(nrow=nrow(consensusCor), ncol = ncol(consensusCor))
for (t in 1:ncol(consensusCor)) {
  tmp <- consensusCor[, t]
  tmp[which(consensusPvalue[,t]>pvalue)] <- NA
  final.modules[ , t] <- tmp
}
final.modules <- as_tibble(final.modules)
colnames(final.modules) <- colnames(consensusCor)
final.modules$"colorME" <- substring(rownames(consensusCor),3)
#rownames(final.modules) <- rownames(consensusCor)
final.modules <- final.modules[which(rowSums(is.na(final.modules[, phenotypes])) < ncol(final.modules[,phenotypes])), ] # remove NA rows

rm(t, set, tmp)

tmp <- dplyr::bind_rows(hubs)
names(tmp) <- paste0("hub.", names(tmp))
#tmp$"colorME" <- paste0("ME",names(hubs[[1]]))
tmp$"colorME" <- names(hubs[[1]])
final.modules <- left_join(final.modules, tmp, by="colorME")

plentyofcolors <- labels2colors(net$color)
#names(plentyofcolors) <- names(net$colors)
plentyofcolors <- as_tibble(plentyofcolors)
plentyofcolors$"gene" <- names(net$colors)
names(plentyofcolors) <- c("colorME", "gene_id")

#tableCol <- as_tibble(table(plentyofcolors$"colorME"))
tableCol <- as.data.frame(table(plentyofcolors$"colorME"))
names(tableCol) <- c("colorME", "ngenes")

final.modules <- left_join(final.modules, tableCol, by="colorME")
colnames(final.modules)[1:4] <- paste0("corr.",colnames(final.modules[1:4]))

write.table(final.modules, file = paste0(workdir, "/final.modules.txt"), append = FALSE, quote = F, sep = "\t", dec = ".", row.names = F, col.names = TRUE)




######################################
# SAVE RESULTS
######################################
save(moduleTraitCor, moduleTraitPvalue, MEColors, consensusCor, consensusPvalue,
     c.geneModuleMembership, MMPvalue, geneTraitSignificance, c.geneTraitSignificance, GSPvalue, c.geneTS.pvalue,
     hubs, one.hubs, plentyofcolors, final.modules, 
     file=paste0(workdir, "/consensus_relating_traits.RData"))

#save(c.geneModuleMembership, geneMM, highMM.genes.ph, kIM, kIM.tot, kIM.scaled,
#     file=paste0(workdir, "/connectivities.RData"))



################################################################################################
# MODULE EIGENGENE 
################################################################################################
setLabels = c("bi", "multi")
softpower = 8

modE <- list()
for (pop in setLabels) {
  modE[[pop]] <- moduleEigengenes(multiExpr[[pop]]$data, 
                                  moduleColors,
                                  impute = TRUE, 
                                  nPC = 10,
                                  align = "along average",
                                  excludeGrey = TRUE,
                                  grey = if (is.numeric(colors)) 0 else "grey",
                                  subHubs = TRUE,
                                  trapErrors = FALSE,
                                  returnValidOnly = FALSE, 
                                  softPower = softpower,
                                  scale = TRUE,
                                  verbose = 0, indent = 0)
} # for pop

#MEs <- multiSetMEs(multiExpr, universalColors = moduleColors, excludeGrey = TRUE);

save(modE, file="module.eigengene.rdata")

extMEs <- list()
for (set in setLabels) {
  extMEs[[set]] = list(data = cbind(PhenoValues[[set]], MEs[[set]]$data))
}

# Check how many eigengenes there are
#lapply(MEs, checkSets) do not work
lapply(extMEs, checkSets)


# Number of modules in each set. (Subtract 1 for the grey module)
nModules =  length(MEs[[set]]$data);
# Initialize lists to hold association statistics
MEsignif = list();
MEsignif.p = list();
MEsignif.Z = list();
MEsignif.metaZ = list();
MEsignif.metap = list();

# Calculate associations
nTraits <- length(names(PhenoValues[[1]]))
MEsignif = MEsignif.p = MEsignif.Z = array(NA, dim = c(nModules, nAnalysisSets, nTraits))
for (set in 1:length(setLabels)) {
  cp = corAndPvalue(MEs[[ set ]]$data, PhenoValues[[set]],alternative = "greater")
  MEsignif[ , set, ] = cp$cor;
  MEsignif.p[ , set, ] = cp$p;
  MEsignif.Z[ , set, ] = cp$Z;
}

dimnames(MEsignif) = dimnames(MEsignif.p) = dimnames(MEsignif.Z) =
  list(colnames(MEs[[ 1 ]] $ data), setLabels, names(PhenoValues[[set]]))

MEsignif.metaZ = apply(MEsignif.Z, c(1,3), sum, na.rm = TRUE)/ sqrt(apply(!is.na(MEsignif.Z), c(1,3), sum))
MEsignif.metap = 2*pnorm(abs(MEsignif.metaZ), lower.tail = FALSE);

collectGarbage()

#modE, extMEs, MEsignif, MEsignif.p, MEsignif.Z, MEsignif.metaZ, MEsignif.metap


save(highMM.mean, high.genes, highMM, highMM.genes.ph, highMM.gTS, highMM.mean2, file= "high.genes.rdata")


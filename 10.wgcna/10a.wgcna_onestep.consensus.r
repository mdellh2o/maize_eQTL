########################################################################################################################
# NETWORK ANALYSIS - WGCNA - Consensus network analysis
# 2 One-step automatic network construction and module detect
########################################################################################################################
R

library(WGCNA)

workdir = "/home/maize.vib/network_splitted_pop/consensus"
#workdir = "/Users/mara/work_analysis/maize.vib/network_splitted_pop/"
setwd(workdir)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

# Load the data saved in the first part
load(file = "/home/maize.vib/network_splitted_pop/biparental/RIL_network_dataInput.bi.RData");
load(file = "/home/maize.vib/network_splitted_pop/multiparental/RIL_network_dataInput.multi.RData");

pop = c("bi", "multi")
nSets = length(pop)
net_type = "signed"


########################################################################################################################
# 1.b Rudimentary data cleaning
########################################################################################################################
# copy expression data into a list
multiExpr = vector(mode = "list", length = nSets)
#multiExpr[[1]] <- exprOnly.bi
#multiExpr[[2]] <- exprOnly.multi
#names(multiExpr) <- pop

genes <- intersect(names(exprOnly.bi), names(exprOnly.multi))
setLabels = c("bi", "multi")

multiExpr[[1]] = list(data = exprOnly.bi[ , genes]);
names(multiExpr[[1]]$data) = names(exprOnly.bi[ , genes]);
rownames(multiExpr[[1]]$data) = rownames(exprOnly.bi[ , genes]);
multiExpr[[2]] = list(data = exprOnly.multi[ , genes]);
names(multiExpr[[2]]$data) = names(exprOnly.multi[ , genes]);
rownames(multiExpr[[2]]$data) = rownames(exprOnly.multi[ , genes]);
names(multiExpr) <- setLabels

exprSize = checkSets(multiExpr) # Check that the data has the correct format for many functions operating on multiple sets:

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK


if (!gsg$allOK) {
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], collapse = ", ")))
  for (set in 1:exprSize$nSets) {
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples", paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes]
  } # for
  # Update exprSize
  exprSize = checkSets(multiExpr)
} #if



# We now cluster the samples on their Euclidean distance, separately in each set
sampleTrees = list()
for (set in 1:nSets) {sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")}


pdf(file = "graphs/SampleClustering.pdf", width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]), xlab="", sub="", cex = 0.7);
dev.off();

# WE DO NOT NEED TO CUT SOME SAMPLES (already did!)



#1.c Loading trait data
###########################################
# Form a multi-set structure that will hold phenotypes
PhenoValues = vector(mode="list", length = nSets);

PhenoValues <- list()
# bi
setSamples = rownames(multiExpr[["bi"]]$data);
traitRows = match(setSamples, pheno.bi$"samples");
PhenoValues[["bi"]] = pheno.bi[traitRows, -1];
#rownames(PhenoValues[[set]]$data) = pheno.bi[traitRows, ];

# multi
setSamples = rownames(multiExpr[["multi"]]$data);
traitRows = match(setSamples, pheno.multi$"samples");
PhenoValues[["multi"]] = pheno.multi[traitRows, -1];
#rownames(PhenoValues[[set]]$data) = pheno.multi[traitRows, ];


# Define data set dimensions
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples


# save results
save(nGenes, nSamples, nSets, PhenoValues, multiExpr, sampleTrees, file = paste0(workdir, "/Consensus-dataInput.RData"))



#####################################################################
#2 Network construction and module detection - ONE-STEP CONSTRUCTION
########################################################################################################################

# RECALCULATE SOFT POWER
# # Choosing the soft-thresholding power: analysis of network topology 
# to choose the soft thresholding power based on the criterion of approximate scale-free topology.
#################################################################################################
powers = c(seq(4,10,by=1), seq(12,20, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets) {
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers, verbose = 2, networkType = net_type)[[2]]);
}
names(powerTables) <- pop
collectGarbage();


# Plot the results:
colors = c("black", "red")
plotCols = c(2,5,6,7) # Will plot these columns of the returned scale free analysis tables "SFT.R.sq", "mean.k", "median.k", "max.k"

pdf(paste0(workdir, "/graphs/choose.softPower.pdf"))
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity", "Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets) {
  for (col in 1:length(plotCols)) {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }#col
} #set

# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)
#pdf(file = "Plots/scaleFreeAnalysis.pdf", wi = 8, he = 6)
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = pop, col = colors, pch = 20) ;
  } else
    legend("topright", legend = pop, col = colors, pch = 20) ;
}
dev.off();



# 2.a.2 NETWORK construction and CONSENSUS module detection
#################################################################################################
softpower = 8
nThreads = 22

net = blockwiseConsensusModules(
  multiExpr, power = softpower, minModuleSize = 30, deepSplit = 2,
  networkType = net_type,
  pamRespectsDendro = FALSE,
  mergeCutHeight = 0.2, numericLabels = TRUE,
  minKMEtoStay = 0,
  saveTOMs = TRUE, verbose = 5, nThreads = nThreads, maxBlockSize = 30000)


consMEs = net$multiMEs;
moduleLabels = net$colors;
# Convert the numeric labels to color labels
moduleColors = labels2colors(moduleLabels)
consTree = net$dendrograms[[1]]

# plot the gene dendrogram and the corresponding module colors
pdf(file = "graphs/ConsensusDendrogram-auto.pdf", wi = 8, he = 6)
plotDendroAndColors(consTree, moduleColors, "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
dev.off()


save(net, net_type, powerTables, consMEs, moduleLabels, moduleColors, consTree, file = "Consensus-NetworkConstruction-auto.RData")




################################################################################################
# 3.b  GENE MODULE MEMBERSHIP - and PVALUE
################################################################################################
# For each module --> define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile. 
# This allows to quantify the similarity of all genes on the array to every module.

# CONNECTIVITY! kME = ModuleMembership
########################################
c.geneModuleMembership = matrix(NA, nrow = ncol(multiExpr[[1]]$data), ncol= nrow(moduleTraitCor[[1]]))
geneMM <- list()
for (set in setLabels) {geneMM[[set]] <- as.data.frame(cor(multiExpr[[set]]$data, consMEs[[set]]$data, use = "p"))}

# very conservative to find consensus geneMM
negative = geneMM[[1]] < 0 & geneMM[[2]] < 0;
c.geneModuleMembership[negative] = pmax(geneMM[[1]][negative], geneMM[[2]][negative]); # very conservative
positive = geneMM[[1]] > 0 & geneMM[[2]] > 0;
c.geneModuleMembership[positive] = pmin(geneMM[[1]][positive], geneMM[[2]][positive]); # very conservative
c.geneModuleMembership <- as.data.frame(c.geneModuleMembership)
setnames(c.geneModuleMembership, gsub("ME", "MM", MEColorNames))
rownames(c.geneModuleMembership) <- names(multiExpr[[1]]$data)

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(c.geneModuleMembership), nSamples)); # Pvalue, each gene with the module
names(MMPvalue) =  gsub("ME", "p.MM", MEColorNames)
rm(negative,positive)
c.geneMM.conservative <-  c.geneModuleMembership


# authors in the pubblication PlosONE 2013 say that consensus geneMM could be calculated as average between sets
geneMM[[1]]$gene_id <- rownames(geneMM[[1]])
geneMM[[2]]$gene_id <- rownames(geneMM[[2]])
dt <- data.frame(row.names=rownames(geneMM[[1]]))
names(dt) <- names(geneMM[[1]])

for (i in names(table(moduleColors))){
  dt[,i] <- apply(cbind(geneMM[[1]][,i],geneMM[[2]][,i]),1,mean)
}
c.geneMM.mean <- dt
rm(dt)




############## Module membership analysis, kME
colorh <- labels2colors(net$colors) # it is the object moduleColors
# colorh --> colors of the genes 
# moduleEigengenes(multiExpr[["bi"]]$data, colorh)$eigengenes  already done data in --> consMEs object
# geneMM is the kME


######################################### 
# Intramodular connectivity (kIM)
######################################### 
colorh <- labels2colors(net$colors) # it is the object moduleColors
softpower = 8
kIM <- list()
#kIM.scaled <- list()
adj <- list()
netk <- list() <-# to store whole network connectivity
  kME <- list()
k <- list()

for (set in setLabels) {
  
  #adj[[set]] <- adjacency(multiExpr[[set]]$data, type = net_type, power=softPower)
  #netk.tmp <- rowSums(adj[[set]])-1 # total network connectivity --> is the same of kIM[[set]]$kTotal, the sum of the adjancency
  #netk[[set]] <- as.data.frame(netk.tmp)
  
  # kTotal, the within module connectivity kWithin, kOut=kTotal-kWithin, and kDiff=kIn-kOut=2*kIN-kTotal
  kIM[[set]]=intramodularConnectivity.fromExpr(multiExpr[[set]]$data, colorh, networkType = net_type,power=softpower)
  #kIM.scaled[[set]]=intramodularConnectivity.fromExpr(multiExpr[[set]]$data, colorh, scaleByMax=TRUE, networkType = net_type,power=softpower)
  
  # kME[[set]] <-  signedKME(multiExpr[[set]]$data, net$multiMEs[[set]]$data) # module membership --> is the same of geneMM[[set]]
  k[[set]] = cbind(kIM[[set]],geneMM[[set]])
  
}



# calculate mean kTOT and kDiff
###############################
plentyofcolors <- labels2colors(net$color)
plentyofcolors <- as_tibble(plentyofcolors)
plentyofcolors$"gene" <- names(net$colors)
names(plentyofcolors) <- c("colorME", "gene_id")
#tableCol <- as_tibble(table(plentyofcolors$"colorME"))
tableCol <- as.data.frame(table(plentyofcolors$"colorME"))
names(tableCol) <- c("colorME", "ngenes")


library(plyr)
# calculate mean kTOT and kDiff
kk <- list()
for (set in setLabels) {kk[[set]] <- cbind(plentyofcolors, kIM[[set]])}
k.mean <- list()
for (set in setLabels) {
  kTotal <- dplyr::summarise(group_by(kk[[set]], colorME), mean=mean(kTotal), sd=sd(kTotal))
  setnames(kTotal, "mean", "kTotal.mean")
  kWithin <- dplyr::summarise(group_by(kk[[set]], colorME), mean=mean(kWithin), sd=sd(kWithin))
  setnames(kWithin, "mean", "kWithin.mean")
  kOut <- dplyr::summarise(group_by(kk[[set]], colorME), mean=mean(kOut), sd=sd(kOut))
  setnames(kOut, "mean", "kOut.mean")
  kDiff <- dplyr::summarise(group_by(kk[[set]], colorME), mean=mean(kDiff), sd=sd(kDiff))
  setnames(kDiff, "mean", "kDiff.mean")
  k.mean[[set]] <- plyr::join_all(list(kTotal,kWithin,kOut,kDiff, tableCol), by='colorME', type='left')
  k.mean[[set]] <- k.mean[[set]] %>% drop_na
}# set
dplyr::summarise(kk[[1]][which(plentyofK[[1]]$colorME != "grey"), ], mean=mean(kTotal), sd=sd(kTotal))
dplyr::summarise(kk[[1]][which(plentyofK[[1]]$colorME != "grey"), ], mean=mean(kWithin), sd=sd(kWithin))
dplyr::summarise(kk[[1]][which(plentyofK[[1]]$colorME != "grey"), ], mean=mean(kOut), sd=sd(kOut))

plot(k.mean[["bi"]]$kTotal.mean, k.mean[["bi"]]$ngenes)
with(text(k.mean[["bi"]]$kTotal.mean, k.mean[["bi"]]$ngenes, labels = k.mean[["bi"]]$colorME))
plot(k.mean[["multi"]]$kTotal.mean, k.mean[["multi"]]$ngenes)
with(text(k.mean[["multi"]]$kTotal.mean, k.mean[["multi"]]$ngenes, labels = k.mean[["bi"]]$colorME))

plot(k.mean[["bi"]]$kDiff.mean, k.mean[["bi"]]$ngenes)
with(text(k.mean[["bi"]]$kDiff.mean, k.mean[["bi"]]$ngenes, labels = k.mean[["bi"]]$colorME))
plot(k.mean[["multi"]]$kDiff.mean, k.mean[["multi"]]$ngenes)
with(text(k.mean[["multi"]]$kDiff.mean, k.mean[["multi"]]$ngenes, labels = k.mean[["bi"]]$colorME))


plentyofK <- list()
plentyofK[["bi"]] <- cbind(plentyofcolors, kIM[[1]])
plentyofK[["multi"]] <- cbind(plentyofcolors, kIM[[2]])
summary(plentyofK[[2]][which(plentyofK[[1]]$colorME != "grey"), "kTotal"])
summary(plentyofK[[1]][which(plentyofK[[1]]$colorME != "grey"), "kTotal"])
summary(plentyofK[[2]][which(plentyofK[[1]]$colorME != "grey"), "kWithin"])
summary(plentyofK[[1]][which(plentyofK[[1]]$colorME != "grey"), "kWithin"])
summary(plentyofK[[2]][which(plentyofK[[1]]$colorME != "grey"), "kOut"])
summary(plentyofK[[1]][which(plentyofK[[1]]$colorME != "grey"), "kOut"])

hist(plentyofK[[2]][which(plentyofK[[1]]$colorME != "grey"), "kDiff"])
hist(plentyofK[[1]][which(plentyofK[[1]]$colorME != "grey"), "kDiff"])

hist(k.mean[[1]][which(k.mean[[1]]$colorME %in% final.modules$colorME), c(1,2,4,6,8,10)])
#'%!in%' <- function(x,y)!('%in%'(x,y))
hist(k.mean[[1]][-c(which(k.mean[[1]]$colorME %in% final.modules$colorME)), c(1,2,4,6,8,10)])


#inModule = which((colorh=="magenta"))
#kIM_magenta=kIM.tot[inModule, ]
#write.table(kIM_yellow,"kIM_yellow.txt")

collectGarbage()

# SCATTER PLOT INTRAMODULAR CONNECTIVITY

save(c.geneModuleMembership, geneMM, highMM.genes.ph, kIM, kIM.tot, kIM.scaled,
     file=paste0(workdir, "/connectivities.RData"))


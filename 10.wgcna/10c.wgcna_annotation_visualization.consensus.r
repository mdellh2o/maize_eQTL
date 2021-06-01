########################################################################################################################
# NETWORK ANALYSIS - WGCNA - ANNOTATION
# add annotation and export network for visualization
# add infos about eQTL
########################################################################################################################
R

library(WGCNA)
library(data.table)
library(dplyr)
library(stringr)

workdir = "/home/maize.vib/network_splitted_pop/consensus"
setwd(workdir)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

# Load the data saved in the first part
load(file = "/home/maize.vib/network_splitted_pop/consensus/Consensus-dataInput.RData") 
load(file = "/home/maize.vib/network_splitted_pop/consensus/Consensus-NetworkConstruction-auto.RData")
load(file = "/home/maize.vib/network_splitted_pop/consensus/consensus_relating_traits.RData")

load(file = "/home/maize.vib/network_splitted_pop/consensus/consensus_relating_traits.RData") # file NETWORK from previous step, #4 consensus
load(file = "/home/maize.vib/eqtl_results_2/me.all.pca10.all.reduced.3.Rdata") # eQTL data 


######################################
# SAVE GENE ID OF EACH MODULE
######################################
# we need the list of the genes and the color assignment --> object plentyofcolors
all.colors <- unique(labels2colors(net$colors))
module_final <- final.modules$"colorME" # modules that are significative correlated with traits
genes_each_module <- list()

for (m in all.colors) {
  tmp.m <- as.data.frame(plentyofcolors[which(plentyofcolors$"colorME" == m), "gene_id"])
  setnames(tmp.m, "gene_id", paste0("gene_id_", m))
  genes_each_module[[m]] <- tmp.m
  
  write.table(tmp.m, file=paste0(workdir, "/gene_id_consensus_", m, ".txt"), quote=F, row.names = F)
} 


# LOAD ANNOTATION
########################################################
# use a gene annotation file to connect gene IDs to gene names and universally recognized identification numbers
# gene ID, module color, gene significance for trait, and module membership and p-values in all modules.
# The modules will be ordered by their significance for trait, with the most significant ones to the left

annotation_file = "/home/maize.vib/annotation/Zea_mays.AGPv4.37.chr.gtf"
#library(rtracklayer)
#annot = readGFF(annotation_file, version = 1)

gtf = fread(annotation_file, data.table=F, stringsAsFactors =F, skip="#!gen")
names(gtf) <- c("seqid", "gene_source", "type", "start", "end", "unkn", "strand", "phase", "group")
gtf = gtf[which(gtf$"seqid" != "Mt" & gtf$"seqid" != "Pt"), c("seqid", "type", "start", "end", "strand", "group")]
#annot$group <- droplevels(annot$group)
gtf$"gene_id" = str_split_fixed(gtf$"group", ";", 2)[ , 1]
gtf$"gene_id" = gsub("gene_id ", "", gtf$"gene_id")
gtf$"gene_id" <-gsub("\"", "", gtf$"gene_id")

group <- str_split(gtf$"group", ";", simplify = TRUE)
gtf$"gene_biotype" <- group[grep("gene_biotype", group)]
gtf$"gene_biotype" = gsub("gene_biotype ", "", gtf$"gene_biotype")
gtf$"gene_biotype" <- gsub("\"", "", gtf$"gene_biotype")

gtf <- gtf[which(gtf$"type" == "gene"), c("gene_id", "type", "seqid", "start", "end", "strand", "gene_biotype")]
names(gtf) = gsub("seqid", "gene_chr",names(gtf))
names(gtf) = gsub("start", "gene_start",names(gtf))
names(gtf) = gsub("end", "gene_end",names(gtf))
gtf <- unique(gtf)

rm(group)

# LOAD EXCEL ANNOTATION and merge
#library(gdata)
#library(dplyr)
#genelistpisa = read.xls("/home/maize.vib/genelistpisa180505.xlsx")
#names(genelistpisa) = gsub("v4_gene_model", "gene_id",names(genelistpisa))
#g <- genelistpisa[ , c("gene_id", "v4_chr", "v4_start", "v4_end", "locus", "entrez_gene", "Transcription.factor.name", "other.gene..", "arabido..orthologs")]
#genelist <- left_join(gtf, g, by = "gene_id") #merge both annotations
#genelist <- genelist[ , c("gene_id", "type", "seqid", "start", "end", "strand", "gene_biotype", "Transcription.factor.name")]
#names(genelist) = gsub("Transcription.factor.name", "transcription_factor",names(genelist))
#genelist <-  unique(genelist)

probes = names(moduleLabels) # names of genes considered in the network construction
probes2annot = match(probes, gtf[which(gtf$"type" == "gene"),]$"gene_id")
p2a = which(is.na(probes2annot)) # number or probes without annotation

gTS <- c.geneTraitSignificance
#names(gTS) <- paste0("geneTS.", names(gTS))
#gTS$"gene_id" <- rownames(c.geneTraitSignificance)

gTSp <- c.geneTS.pvalue
names(gTSp) <- paste0("GSPvalue.", names(c.geneTS.pvalue))
gTSp$"gene_id" <- rownames(c.geneTS.pvalue)


gtf <- left_join(as.data.frame(plentyofcolors), gtf, by="gene_id", all.y=F)
gtf <- left_join(gtf, gTS, by="gene_id", all.y=F)
gtf <- left_join(gtf, gTSp, by="gene_id", all.y=F)

rm(gTS, gTSp)


# add INFO on HUB genes (not used in the paper)
########################################################
colnames(gtf)[2] <- "gene"
strega <- left_join(int.eqtl, gtf[, c(1:2,7:16)], by="gene")  



# EXPORTING results of the network analysis
########################################################################################################################
# Put together a data frame that summarizes the results of network analysis, namely the gene significances (GS) and module memberships (also known as kME) of all probes
c.geneModuleMembership$"gene" <- rownames(c.geneModuleMembership)
c.geneTraitSignificance$"gene" <- rownames(c.geneTraitSignificance)
export <- left_join(c.geneModuleMembership, c.geneTraitSignificance, by="gene") 



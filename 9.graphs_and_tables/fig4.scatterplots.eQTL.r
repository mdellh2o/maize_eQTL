################################################################################################################
# MAKE SCATTER PLOTS FROM eQTL RESULTS - FIG 4
################################################################################################################
R
library(data.table)
library(ggplot2)

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

traits<-c("leaf.length", "DZ.size", "LER", "LED")
setwd(workdir)

load(eqtl_file_name) # objects: me, eqtl_maf, pcc.eqtl, gwa.eqtl (Bonferroni threshold)
#eqtl <- fread(eqtl_df, data.table=F) # or eqtl <- me.all.pca10.2$all$eqtls


# add gene positions at eqtl file
#genepos <- fread(gene_position_file)
#names(genepos) <- c("gene", "gene.chr", "gene.start", "gene.end")
#eqtl <- merge(eqtl, genepos, by = "gene")
#eqtl <- as.data.frame(eqtl)

#load colors
#load("/home/maize.vib/PUBBLICATION/colors.by.phenotype.Rdata")
load("/Users/mara/work_analysis/colors.by.phenotype.Rdata")
phenocol<-phenocol[c(2,1,4,3),]
stopifnot(all(phenocol[,1]==traits))
phenocol <- as.data.frame(phenocol)

# prepare data
gwa.leaf.length <- gwa.eqtl[which(gwa.eqtl$gwa.leaf.length > 0), ]
gwa.LER <- gwa.eqtl[which(gwa.eqtl$gwa.LER > 0), ]
gwa.LED <- gwa.eqtl[which(gwa.eqtl$gwa.LED > 0), ]
gwa.DZ.size <- gwa.eqtl[which(gwa.eqtl$gwa.DZ.size > 0), ]
# prepare data
pcc.leaf.length <- pcc.eqtl[which(pcc.eqtl$pcc.leaf.length > 0), ]
pcc.LER <- pcc.eqtl[which(pcc.eqtl$pcc.LER > 0), ]
pcc.LED <- pcc.eqtl[which(pcc.eqtl$pcc.LED > 0), ]
pcc.DZ.size <- pcc.eqtl[which(pcc.eqtl$pcc.DZ.size > 0), ]



#FIX POSITIONS
####################################
#get in chr length info
#chrlen<-read.table("/home/maize.vib/annotation/chr.dimensions.ref4.txt", header=T)
chrlen<-read.table("C:/Users/Mara/OneDrive/eQTL_maize/data/chr.dimensions.ref4.txt", header=T)

#make all in Mb
if (max(chrlen[,"length"])>1e6){
  chrlen[,"length"]<-chrlen[,"length"]*1e-6
}

#get cumulative bp position by chr
cum<-cumsum(chrlen[,2])
chrlen$startingpos<-c(0,cum[1:(length(cum)-1)])
chrlen$endpos<-cumsum(chrlen[,2])


#produce xy scatterplots for PCC-eQTL
gwa.pos.list <-  list()
gwa.pos.list[["leaf.length"]] <- gwa.leaf.length
gwa.pos.list[["LER"]] <- gwa.LER
gwa.pos.list[["LED"]] <- gwa.LED
gwa.pos.list[["DZ.size"]] <- gwa.DZ.size

pcc.pos.list <- list()
pcc.pos.list[["leaf.length"]] <- pcc.leaf.length
pcc.pos.list[["LER"]] <- pcc.LER
pcc.pos.list[["LED"]] <- pcc.LED
pcc.pos.list[["DZ.size"]] <- pcc.DZ.size



scat.plot <- list()
#pdf("Fig.4.pdf")
tiff("Fig.4.tiff", width = 9.5, height = 10, units = 'in', res = 300)
#grid::grid.newpage()
par(mfrow=c(2,2))
for (f in traits){ #f=1
  pcc.pos <- pcc.pos.list[[f]]
  gwa.pos <- gwa.pos.list[[f]]
  tmptrait <- f
  print(tmptrait)
  #bring everyting in Mb
  if (max(pcc.pos[,"snp.pos"])>1e6){
    pcc.pos[,"snp.pos"]<-pcc.pos[,"snp.pos"]*1e-6
    pcc.pos[,"gene.start"]<-pcc.pos[,"gene.start"]*1e-6
    pcc.pos[,"gene.end"]<-pcc.pos[,"gene.end"]*1e-6
  }
  
  #convert positions on a genomic scale
  #line by line
  snpgeno<-c()
  genegeno<-c()
  for (i in 1:nrow(pcc.pos)){
    snpgeno[i]<-pcc.pos[i,"snp.pos"]+chrlen[pcc.pos[i,"snp.chr"],"startingpos"]
    genegeno[i]<-pcc.pos[i,"gene.start"]+chrlen[pcc.pos[i,"gene.chr"],"startingpos"]
  }#for i
  toplot.pcc<-data.frame(pcc.pos[,"snps"], snpgeno, pcc.pos[,"gene"], genegeno,
                     pcc.pos[,"pvalue"], pcc.pos[,"beta"], pcc.pos[,"cistrans"])
  toplot.pcc<-toplot.pcc[!duplicated(toplot.pcc[,c(1,3)]),]
  
  
  if (max(gwa.pos[,"snp.pos"])>1e6){
    gwa.pos[,"snp.pos"]<-gwa.pos[,"snp.pos"]*1e-6
    gwa.pos[,"gene.start"]<-gwa.pos[,"gene.start"]*1e-6
    gwa.pos[,"gene.end"]<-gwa.pos[,"gene.end"]*1e-6
  }
  
  snpgeno<-c()
  genegeno<-c()
  for (i in 1:nrow(gwa.pos)){
    snpgeno[i]<-gwa.pos[i,"snp.pos"]+chrlen[gwa.pos[i,"snp.chr"],"startingpos"]
    genegeno[i]<-gwa.pos[i,"gene.start"]+chrlen[gwa.pos[i,"gene.chr"],"startingpos"]
  }#for i
  toplot.gwa<-data.frame(gwa.pos[,"snps"], snpgeno, gwa.pos[,"gene"], genegeno,
                     gwa.pos[,"pvalue"], gwa.pos[,"beta"], gwa.pos[,"cistrans"])
  
  #get info about genes shared with PCC.eQTL so to plot them
  #tmpshare<-shareg[which(shareg[,1]==tmptrait),]
  #toplotshare<-toplot[toplot[,3] %in% tmpshare[,2],]
  
  #remove duplicated hits from both dataframes
  toplot.gwa<-toplot.gwa[!duplicated(toplot.gwa[,c(1,3)]),]
  #toplotshare<-toplotshare[!duplicated(toplotshare[,c(1,3)]),]
  
  #make the plot
  #pdf(paste(tmptrait, ".xyscatter.pcc_eQTL.pdf", sep=""))
  #pdf(NULL)
  #dev.control(displaylist="enable")
  
  plot(x=toplot.pcc[,"snpgeno"],y=toplot.pcc[,"genegeno"],ylab="Genomic position of Transcript (Mb)", 
       xlab="Genomic position of eQTL (Mb)", main=tmptrait, type="n",
       xlim=c(0, max(chrlen[,4])), ylim=c(0, max(chrlen[,4])),  xaxs="i",yaxs="i")
  rect(xleft=chrlen[1,4], xright=chrlen[2,4], ybottom=0, ytop=max(chrlen[,4]), col = "gray96", border = NA)
  rect(xleft=chrlen[3,4], xright=chrlen[4,4], ybottom=0, ytop=max(chrlen[,4]), col = "gray96", border = NA)
  rect(xleft=chrlen[5,4], xright=chrlen[6,4], ybottom=0, ytop=max(chrlen[,4]), col = "gray96", border = NA)
  rect(xleft=chrlen[7,4], xright=chrlen[8,4], ybottom=0, ytop=max(chrlen[,4]), col = "gray96", border = NA)
  rect(xleft=chrlen[9,4], xright=chrlen[10,4], ybottom=0, ytop=max(chrlen[,4]), col = "gray96", border = NA)
  abline(v=chrlen[1:9,4], col="gray96", lty=3)
  #add pccs hits
  lines(x = c(0, max(chrlen[,4])), y = c(0, max(chrlen[,4])), col="gray60")
  points(x=toplot.pcc[,2],y=toplot.pcc[,4], pch=20, cex=0.8, col="gray40")
  color.trait <- as.character(phenocol[which(phenocol$phenotype == f), 2])
  points(x=toplot.gwa[,2],y=toplot.gwa[,4], pch=21, col=color.trait, cex=1.5)
  
  box()
  
  #dev.off()
  #scat.plot[[f]] <- recordPlot()
  #invisible(dev.off())
  
}#for f
dev.off()




#produce xy scatterplots for TOT-eQTL
#####################################
eqtl.pos  <- eqtl[which(eqtl$tot.gwa > 0 | eqtl$tot.pcc > 0 ), ]


tiff("Fig.S5.tiff", width = 5, height = 5, units = 'in', res = 300)
#grid::grid.newpage()
par(mfrow=c(1,1))
  #bring everyting in Mb
  if (max(eqtl.pos[,"snp.pos"])>1e6){
    eqtl.pos[,"snp.pos"]<-eqtl.pos[,"snp.pos"]*1e-6
    eqtl.pos[,"gene.start"]<-eqtl.pos[,"gene.start"]*1e-6
    eqtl.pos[,"gene.end"]<-eqtl.pos[,"gene.end"]*1e-6
  }
  
  #convert positions on a genomic scale
  #line by line
  snpgeno<-c()
  genegeno<-c()
  for (i in 1:nrow(eqtl.pos)){
    snpgeno[i]<-eqtl.pos[i,"snp.pos"]+chrlen[eqtl.pos[i,"snp.chr"],"startingpos"]
    genegeno[i]<-eqtl.pos[i,"gene.start"]+chrlen[eqtl.pos[i,"gene.chr"],"startingpos"]
  }#for i
  toplot.eqtl<-data.frame(eqtl.pos[,"snps"], snpgeno, eqtl.pos[,"gene"], genegeno,
                         eqtl.pos[,"pvalue"], eqtl.pos[,"beta"], eqtl.pos[,"cistrans"])
  toplot.eqtl<-toplot.eqtl[!duplicated(toplot.eqtl[,c(1,3)]),]
  
  
  plot(x=toplot.eqtl[,"snpgeno"],y=toplot.eqtl[,"genegeno"],ylab="Genomic position of Transcript (Mb)", 
       xlab="Genomic position of eQTL (Mb)", main="significative GWA.eQTL and PCC.eQTL", type="n",
       xlim=c(0, max(chrlen[,4])), ylim=c(0, max(chrlen[,4])),  xaxs="i",yaxs="i")
  rect(xleft=chrlen[1,4], xright=chrlen[2,4], ybottom=0, ytop=max(chrlen[,4]), col = "gray96", border = NA)
  rect(xleft=chrlen[3,4], xright=chrlen[4,4], ybottom=0, ytop=max(chrlen[,4]), col = "gray96", border = NA)
  rect(xleft=chrlen[5,4], xright=chrlen[6,4], ybottom=0, ytop=max(chrlen[,4]), col = "gray96", border = NA)
  rect(xleft=chrlen[7,4], xright=chrlen[8,4], ybottom=0, ytop=max(chrlen[,4]), col = "gray96", border = NA)
  rect(xleft=chrlen[9,4], xright=chrlen[10,4], ybottom=0, ytop=max(chrlen[,4]), col = "gray96", border = NA)
  abline(v=chrlen[1:9,4], col="gray96", lty=3)
  #add eqtls hits
  lines(x = c(0, max(chrlen[,4])), y = c(0, max(chrlen[,4])), col="gray80")
  points(x=toplot.eqtl[,2],y=toplot.eqtl[,4], pch=20, cex=0.8)

  box()
dev.off()

# eqtl are already selected by Bonferroni threshold
if (max(eqtl[,"snp.pos"])>1e6){
  eqtl$"snp.pos.cor"<-eqtl[,"snp.pos"]*1e-6
  eqtl$"gene.start.cor"<-eqtl[,"gene.start"]*1e-6
  eqtl$"gene.end.cor"<-eqtl$"gene.end"*1e-6
}


#convert positions on a genomic scale
#snpgeno<-c()
#genegeno<-c()
eqtl$"snp.pos.cor" <- eqtl[, "snp.pos.cor"] + chrlen[eqtl[,"snp.chr"],"startingpos"]
eqtl$"gene.start.cor" <- eqtl[, "gene.start.cor"]+chrlen[eqtl[,"gene.chr"],"startingpos"]
eqtl$"gene.end.cor" <- eqtl[, "gene.end.cor"]+chrlen[eqtl[,"gene.chr"],"startingpos"]


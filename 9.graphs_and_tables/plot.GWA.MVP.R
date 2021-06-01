maindir<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/MAIZE_ghent/analyses/3.GWAS.MVP/best.model.out/from.server"
setwd(maindir)
library(MVP)

#set thrs
llth4r<-6.770000E-07
dzthr<-5.150000E-07
lerthr<-1.070000E-06
ledthr<-4.630000E-07

imMVP1<-read.csv("MVP.leaf.length.FarmCPU.csv")
imMVP2<-read.csv("MVP.DZ.size.FarmCPU.csv")
imMVP3<-read.csv("MVP.LER.FarmCPU.csv")
imMVP4<-read.csv("MVP.LED.FarmCPU.csv")

imMVP<-data.frame(imMVP1[,c(1,2,3,5)], imMVP2[,5], imMVP3[,5],imMVP4[,5])
colnames(imMVP)[4:ncol(imMVP)]<-c("Leaf length", "DZ size", "LER", "LED")

#plot GWAS 
MVP.Report(imMVP, plot.type="m", LOG10=TRUE, ylim=NULL,
		threshold=ledthr,
		threshold.lty=2,
		multracks=TRUE,
        col=c("grey60","grey30"), 
		threshold.lwd=1	, 
		threshold.col="black", 
		amplify=TRUE,
        chr.den.col=c("darkgreen", "yellow", "red"),
		bin.size=1e6,
		signal.col="red",
        signal.cex=1,
		signal.pch=19,
		file="tiff",
		memo="",
		dpi=900)

		
#plot density
MVP.Report(imMVP[,1:3], plot.type="d", col=c("darkgreen", "yellow", "red"), file="tiff", dpi=900)
		
#plot QQ plots
MVP.Report(imMVP,plot.type="q",col=c("#b6ba27","#0095f8","#ff6459","#c0258c"),
		threshold=lerthr,
        signal.pch=19,signal.cex=1.5,signal.col="red",
		conf.int.col="grey",box=FALSE,multracks=TRUE,
		file="tiff",memo="",dpi=900)

#plot a PCA
genotype <- attach.big.matrix("../mvp.hmp.geno.desc")
pcafull <- prcomp(t(as.matrix(genotype)))
pca<-pcafull$x[, 1:3]

phenotype <- read.table("../mvp.hmp.phe",head=TRUE)
#create plotting object
rownames(pca)<-phenotype[,1]
pca<-pca[rev(order(rownames(pca))),]
pca<-pca[c(2:nrow(pca), 1),]
colzpca<-c(rep("#5f0079",103), rep("#b32618", 94),  rep("#46c35b",8))

PoV <- pcafull$sdev^2/sum(pcafull$sdev^2)
#plot it
pdf("PCA.plot.genotypes.pdf")
	plot(pca[,1], pca[,2], col=colzpca, pch=19, xlab=paste("PC1 ", round(PoV[1]*100,1), "%", sep=""), ylab=paste("PC2 ", round(PoV[2]*100,1), "%", sep=""))
	legend("bottomleft", legend=c("2-way", "MAGIC", "Parentals"), col=c("#5f0079", "#b32618", "#46c35b"), pch=19, cex=0.8)
dev.off()


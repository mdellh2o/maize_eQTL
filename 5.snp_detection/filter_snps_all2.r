# extract genotypes
#######################################################################################################################################################################################
#https://knausb.github.io/vcfR_documentation/visualization_1.html
R
#install.packages('vcfR')
library(vcfR)
library(data.table)
library(dplyr)

# VARIABLES
dna_file="/home/m.miculan/reference/zea_mays/Zea_mays.AGPv4.dna.toplevel.fa"
gtf_file="/home/m.miculan/reference/zea_mays/annotation/Zea_mays.AGPv4.37.chr.gtf"	

ril.multi=c("NG-6860_R_1","NG-6860_R_10","NG-6860_R_11","NG-6860_R_12","NG-6860_R_13","NG-6860_R_14","NG-6860_R_15","NG-6860_R_16","NG-6860_R_18","NG-6860_R_19","NG-6860_R_2","NG-6860_R_20","NG-6860_R_22","NG-6860_R_23","NG-6860_R_24","NG-6860_R_26","NG-6860_R_27","NG-6860_R_28","NG-6860_R_29","NG-6860_R_3","NG-6860_R_30","NG-6860_R_31","NG-6860_R_32","NG-6860_R_33","NG-6860_R_34","NG-6860_R_35","NG-6860_R_36","NG-6860_R_37","NG-6860_R_38","NG-6860_R_39","NG-6860_R_4","NG-6860_R_40","NG-6860_R_41","NG-6860_R_42","NG-6860_R_44","NG-6860_R_45","NG-6860_R_46","NG-6860_R_47","NG-6860_R_48","NG-6860_R_49","NG-6860_R_5","NG-6860_R_50","NG-6860_R_51","NG-6860_R_52","NG-6860_R_6","NG-6860_R_7","NG-6860_R_8","NG-6860_R_9","NG-7084_R_53","NG-7084_R_54","NG-7084_R_55","NG-7084_R_57","NG-7084_R_58","NG-7084_R_59","NG-7084_R_60","NG-7084_R_61","NG-7084_R_62","NG-7084_R_63","NG-7084_R_64","NG-7084_R_65","NG-7084_R_66","NG-7084_R_67","NG-7084_R_68","NG-7084_R_69","NG-7084_R_70","NG-7084_R_71","NG-7084_R_72","NG-7084_R_73","NG-7084_R_74","NG-7084_R_75","NG-7084_R_76","NG-7084_R_77","NG-7084_R_78","NG-7084_R_79","NG-7084_R_80","NG-7084_R_81","NG-7084_R_82","NG-7084_R_83","NG-7084_R_84","NG-7084_R_85","NG-7084_R_86","NG-7084_R_87","NG-7084_R_88","NG-7084_R_89","NG-7084_R_90","NG-7084_R_91","NG-7084_R_92","NG-7084_R_93","NG-7084_R_94","NG-7084_R_95","NG-7084_R_96","NG-7084_R_97","NG-7084_R_98","NG-7084_R_99")
ril.bi=c("NG-6280_RIL138_P1_lib15140","NG-6280_RIL83_P1_lib15144","NG-6280_RIL94_P1_lib15138","NG-6280_RIL96_P1_lib15142", "NG-6564_RIL100_lib19904","NG-6564_RIL103_lib19905","NG-6564_RIL104_lib19906","NG-6564_RIL105_lib19907","NG-6564_RIL106_lib19908","NG-6564_RIL107_lib19909","NG-6564_RIL109_lib19910","NG-6564_RIL111_lib19911","NG-6564_RIL113_lib19912","NG-6564_RIL114_lib19913","NG-6564_RIL115_lib19914","NG-6564_RIL116_lib19915","NG-6564_RIL117_lib19916","NG-6564_RIL118_lib19917","NG-6564_RIL119_lib19918","NG-6564_RIL11_lib19839","NG-6564_RIL120_lib19919","NG-6564_RIL121_lib19920","NG-6564_RIL122","NG-6564_RIL123_lib19922","NG-6564_RIL124_lib19923","NG-6564_RIL125_lib19924","NG-6564_RIL126_lib19925","NG-6564_RIL127_lib19926","NG-6564_RIL128_lib19927","NG-6564_RIL130_lib19928","NG-6564_RIL131_lib19929","NG-6564_RIL132_lib19930","NG-6564_RIL133_lib19931","NG-6564_RIL134_lib19932","NG-6564_RIL135_lib19933","NG-6564_RIL136_lib19934","NG-6564_RIL137_lib19935","NG-6564_RIL139_lib19936","NG-6564_RIL13_lib19840","NG-6564_RIL140_lib19937","NG-6564_RIL141_lib19938","NG-6564_RIL16_lib19841","NG-6564_RIL17_lib19842","NG-6564_RIL18_lib19843","NG-6564_RIL19_lib19844","NG-6564_RIL23_lib19845","NG-6564_RIL25_lib19846","NG-6564_RIL26_lib19847","NG-6564_RIL27_lib19848","NG-6564_RIL29_lib19849","NG-6564_RIL30_lib19850","NG-6564_RIL32_lib19851","NG-6564_RIL33_lib19852","NG-6564_RIL36_lib19853","NG-6564_RIL37_lib19854","NG-6564_RIL38_lib19855","NG-6564_RIL39_lib19856","NG-6564_RIL40_lib19857","NG-6564_RIL42_lib19858","NG-6564_RIL43_lib19859","NG-6564_RIL44_lib19860","NG-6564_RIL45_lib19861","NG-6564_RIL46_lib19862","NG-6564_RIL47_lib19863","NG-6564_RIL48_lib19864","NG-6564_RIL49_lib19865","NG-6564_RIL50_lib19866","NG-6564_RIL51_lib19867","NG-6564_RIL52_lib19868","NG-6564_RIL53_lib19869","NG-6564_RIL54_lib19870","NG-6564_RIL55_lib19871","NG-6564_RIL57_lib19872","NG-6564_RIL58_lib19873","NG-6564_RIL59_lib19874","NG-6564_RIL60_lib19875","NG-6564_RIL61_lib19876","NG-6564_RIL62_lib19877","NG-6564_RIL63_lib19878","NG-6564_RIL64_lib19879","NG-6564_RIL65_lib19880","NG-6564_RIL66_lib19881","NG-6564_RIL67_lib19882","NG-6564_RIL69_lib19883","NG-6564_RIL70_lib19884","NG-6564_RIL71_lib19885","NG-6564_RIL72_lib19886","NG-6564_RIL73_lib19887","NG-6564_RIL74_lib19888","NG-6564_RIL75","NG-6564_RIL76_lib19890","NG-6564_RIL80_lib19891","NG-6564_RIL81_lib19892","NG-6564_RIL82_lib19893","NG-6564_RIL84_lib19894","NG-6564_RIL85_lib19895","NG-6564_RIL88_lib19896","NG-6564_RIL89_lib19897","NG-6564_RIL90_lib19898","NG-6564_RIL91_lib19899","NG-6564_RIL92_lib19900","NG-6564_RIL93_lib19901","NG-6564_RIL95_lib19902","NG-6564_RIL99_lib19903")
parents=c("A632", "B73", "CML91", "F7", "H99", "HP301", "Mo17", "W153R")

vcf_file_name="/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/vcf_tmp/snps.all.UG.goodQUAL.08.filtered.tmp2.vcf"

setwd("/home/m.miculan/zea_mays_projects/eQTL/SNPs")


vcf=read.vcfR(vcf_file_name) 

gt=extract.gt(vcf, element = "GT", mask = FALSE, as.numeric = FALSE, return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE,convertNA = TRUE)

# substitute genotypes
gt=gsub("0/0", "0", gt)
gt=gsub("0/1", "1", gt)
gt=gsub("1/1", "2", gt)
gt=as.data.frame(gt)


# align names to be coherent with expression file 
names(gt)=gsub("_1178_5", "", colnames(gt))
names(gt)=gsub("_1418_3", "", colnames(gt))
names(gt)=gsub("_1418_4", "", colnames(gt))
names(gt)=gsub("_1418_5", "", colnames(gt))
names(gt)=gsub("_1418_6", "", colnames(gt))
names(gt)=gsub("_1418_7", "", colnames(gt))
names(gt)=gsub("_1479_1", "", colnames(gt))
names(gt)=gsub("_1479_2", "", colnames(gt))
names(gt)=gsub("_1479_3", "", colnames(gt))
names(gt)=gsub("_1479_4", "", colnames(gt))
names(gt)=gsub("_1479_6", "", colnames(gt))


gt<-cbind(row.names(gt),gt)
names(gt)[1]<-"SNP"
gt=gt[c("SNP", parents, ril.bi, ril.multi)]


# VALIDATION MULTI- marker (through BIPARENTAL POPULATION) - SAVE GENOTYPES MULTI and FILTERED, BIALLELIC
#######################################################################################################################################################################################
output_file_name="/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/genotypes.chr.all.filtered.txt"

# select BIPARENTAL
gt.bi=gt[ , c("SNP", ril.bi) ]
gt.bi[,1]<-NULL
gt.bi=as.data.frame(t(gt.bi))


# count genotype class each marker, consider even NA because we need the percentage 
count0<-apply(apply(gt.bi ,2,"==",0),2,sum,na.rm=T)
count1<-apply(apply(gt.bi,2,"==",1),2,sum,na.rm=T)
count2<-apply(apply(gt.bi,2,"==",2),2,sum,na.rm=T)
countNA<-apply(apply(gt.bi,2,is.na),2,sum,na.rm=T)
conte<-cbind(count0,count1,count2,countNA)

# proportion each genotype class
prop=as.data.frame(prop.table(as.matrix(conte), margin = 1))
prop<-cbind(row.names(prop),prop)
names(prop)[1]<-"SNP"
names(prop)=gsub("count", "prop", colnames(prop))

gt.all=as.data.frame(cbind(gt, prop))

# delete the ones that should be homoz reference
# non funziona

#no.good=gt.all[(gt.all$"B73"==0 & gt.all$"H99"==0 & (gt.all$"prop1" + gt.all$"prop2" > 0.02)), ]
#no.good=no.good$"SNP"

# null!

gt.all=gt.all[ , c("SNP", parents, ril.bi, ril.multi) ]

write.table(gt.all, file = output_file_name , sep="\t", na="NA", row.names=F, col.names=T, quote = F)



# BIPARENTAL POPULATION
#######################################################################################################################################################################################
# input_file_name="/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/genotypes.chr.all.filtered.txt"


# take advantage of calculation made in the previous step about genotypes frequency calculation
# use objects: 
# prop --> is the proportion of class of genotypes
# gt.bi --> is a traspose of gt.all only with ril.bi population
# gt.all--> all genotypes


riprova0<-prop$"prop0"<=0.1
riprova1<-prop$"prop1"<=0.1
riprova2<-prop$"prop2"<=0.1

gt.all$"riprova0"<-riprova0
gt.all$"riprova1"<-riprova1
gt.all$"riprova2"<-riprova2


for (aaa in (length(parents)+2):ncol(gt.all))
{
    cat(aaa,"\n")
    gt.all[,aaa]<-replace(gt.all[,aaa], riprova0&gt.all[,aaa]==0,NA)
    gt.all[,aaa]<-replace(gt.all[,aaa], riprova1&gt.all[,aaa]==1,NA)
    gt.all[,aaa]<-replace(gt.all[,aaa], riprova2&gt.all[,aaa]==2,NA)
}


# calculate again proportions of genotypes after masking rare genotypes
gt.bi=gt.all[c("SNP", ril.bi)]
row.names(gt.bi)<-gt.bi[,1]
gt.bi[,1]<-NULL
gt.bi=as.data.frame(t(gt.bi))

newcount0<-apply(apply(gt.bi ,2,"==",0),2,sum,na.rm=T)
newcount1<-apply(apply(gt.bi,2,"==",1),2,sum,na.rm=T)
newcount2<-apply(apply(gt.bi,2,"==",2),2,sum,na.rm=T)
newcountNA<-apply(apply(gt.bi,2,is.na),2,sum,na.rm=T)
newconte<-cbind(newcount0,newcount1,newcount2,newcountNA)
newprop=as.data.frame(prop.table(as.matrix(newconte), margin = 1))
names(newprop)=gsub("count", "prop", colnames(newprop))


gt.bi2=cbind(gt.all[c("SNP", "B73", "H99", ril.bi)], newprop)

geno5=gt.bi2[gt.bi2$"newpropNA"<0.2, ]


# delete only one genotipic class, no informative
geno5=geno5[(geno5$"newprop0" + geno5$"newpropNA")!=1, ]
geno5=geno5[(geno5$"newprop1" + geno5$"newpropNA")!=1, ]
geno5=geno5[(geno5$"newprop2" + geno5$"newpropNA")!=1, ]
geno5=geno5[geno5$"newprop2">0, ] # should be at least some markers 1/1 


# select BIPARENTAL population:
# delete markers that are heterozygous in parents
bip=geno5[ , c("SNP", "B73", "H99")]
bip[,1]<-NULL
bip=as.data.frame(t(bip))
count0<-apply(apply(bip ,2,"==",0),2,sum,na.rm=T)
count1<-apply(apply(bip ,2,"==",1),2,sum,na.rm=T)
count2<-apply(apply(bip ,2,"==",2),2,sum,na.rm=T)
countNA<-apply(apply(bip,2,is.na),2,sum,na.rm=T)


biparental=cbind(geno5 [ ,c("SNP", "B73", "H99" , ril.bi)],count0,count1,count2,countNA)
biparental=biparental[ biparental$"count0"!=2, ]
biparental=biparental[ biparental$"count1"!=2, ]
biparental=biparental[ biparental$"count2"!=2, ]
biparental=biparental[ biparental$"countNA"!=2, ]
#biparental=biparental[ biparental$"countNA"!=2, ]

biparental=biparental[ , c("SNP", "B73", "H99", ril.bi)]




# in teoria quello sotto non mi serve, ma controllare!!!
# e copiare file., controlalre di aver tolto tutto, salvare biparental nel file output e poi salvare estrarre le coordinate! 



# delete het markers in PARENTS 
#?????
p=geno5[c("SNP", parents)]
row.names(p)<-p[,1]
p[,1]<-NULL
p=as.data.frame(t(p))
count_p0<-apply(apply(p ,2,"==",0),2,sum,na.rm=T)
count_p1<-apply(apply(p,2,"==",1),2,sum,na.rm=T)
count_p2<-apply(apply(p,2,"==",2),2,sum,na.rm=T)
count_pNA<-apply(apply(p,2,is.na),2,sum,na.rm=T)
newconte<-cbind(count_p0,count_p1,count_p2,count_pNA)
prop_par=as.data.frame(prop.table(as.matrix(newconte), margin = 1))
names(prop_par)=gsub("count", "prop", colnames(prop_par))

geno5=cbind(geno5, prop_par)

geno5=geno5[(geno5$"prop_p1" + geno5$"prop_pNA")!=1, ]
geno5=geno5[(geno5$"prop_p0" + geno5$"prop_pNA")!=1, ]



# SELECT markers:
write.table(geno5, file = "/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/genotypes.chr.multi.filtered.txt" , sep="\t", na="NA", row.names=F, col.names=T, quote = F)

output_file_name="/home/m.miculan/zea_mays_projects/eQTL/SNPs/unifiedgenotyper/genotypes.chr.bi.filtered2.txt"
write.table(biparental, file = output_file_name , sep="\t", na="NA", row.names=F, col.names=T, quote = F)


quit()



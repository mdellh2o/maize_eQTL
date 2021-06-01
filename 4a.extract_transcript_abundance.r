#######################################################################################################################################################################################
# 6.  Extract FPKM/TPM from StringTie. Write a table with genes name and unstranded raw counts.
#######################################################################################################################################################################################
R

#install.packages("dplyr")
library(data.table)
library(dplyr)

#infiles_exp1=c("NG-6280_RIL138_P1","NG-6280_RIL138_P2","NG-6280_RIL83_P1","NG-6280_RIL83_P2","NG-6280_RIL94_P1","NG-6280_RIL94_P2","NG-6280_RIL96_P1","NG-6280_RIL96_P2","NG-6564_RIL100_lib19904","NG-6564_RIL103_lib19905","NG-6564_RIL104_lib19906","NG-6564_RIL105_lib19907","NG-6564_RIL106_lib19908","NG-6564_RIL107_lib19909","NG-6564_RIL109_lib19910","NG-6564_RIL111_lib19911","NG-6564_RIL113_lib19912","NG-6564_RIL114_lib19913","NG-6564_RIL115_lib19914","NG-6564_RIL116_lib19915","NG-6564_RIL117_lib19916","NG-6564_RIL118_lib19917","NG-6564_RIL119_lib19918","NG-6564_RIL11_lib19839","NG-6564_RIL120_lib19919","NG-6564_RIL121_lib19920","NG-6564_RIL122","NG-6564_RIL123_lib19922","NG-6564_RIL124_lib19923","NG-6564_RIL125_lib19924","NG-6564_RIL126_lib19925","NG-6564_RIL127_lib19926","NG-6564_RIL128_lib19927","NG-6564_RIL130_lib19928","NG-6564_RIL131_lib19929","NG-6564_RIL132_lib19930","NG-6564_RIL133_lib19931","NG-6564_RIL134_lib19932","NG-6564_RIL135_lib19933","NG-6564_RIL136_lib19934","NG-6564_RIL137_lib19935","NG-6564_RIL139_lib19936","NG-6564_RIL13_lib19840","NG-6564_RIL140_lib19937","NG-6564_RIL141_lib19938","NG-6564_RIL16_lib19841","NG-6564_RIL17_lib19842","NG-6564_RIL18_lib19843","NG-6564_RIL19_lib19844","NG-6564_RIL23_lib19845","NG-6564_RIL25_lib19846","NG-6564_RIL26_lib19847","NG-6564_RIL27_lib19848","NG-6564_RIL29_lib19849","NG-6564_RIL30_lib19850","NG-6564_RIL32_lib19851","NG-6564_RIL33_lib19852","NG-6564_RIL36_lib19853","NG-6564_RIL37_lib19854","NG-6564_RIL38_lib19855","NG-6564_RIL39_lib19856","NG-6564_RIL40_lib19857","NG-6564_RIL42_lib19858","NG-6564_RIL43_lib19859","NG-6564_RIL44_lib19860","NG-6564_RIL45_lib19861","NG-6564_RIL46_lib19862","NG-6564_RIL47_lib19863","NG-6564_RIL48_lib19864","NG-6564_RIL49_lib19865","NG-6564_RIL50_lib19866","NG-6564_RIL51_lib19867","NG-6564_RIL52_lib19868","NG-6564_RIL53_lib19869","NG-6564_RIL54_lib19870","NG-6564_RIL55_lib19871","NG-6564_RIL57_lib19872","NG-6564_RIL58_lib19873","NG-6564_RIL59_lib19874","NG-6564_RIL60_lib19875","NG-6564_RIL61_lib19876","NG-6564_RIL62_lib19877","NG-6564_RIL63_lib19878","NG-6564_RIL64_lib19879","NG-6564_RIL65_lib19880","NG-6564_RIL66_lib19881","NG-6564_RIL67_lib19882","NG-6564_RIL69_lib19883","NG-6564_RIL70_lib19884","NG-6564_RIL71_lib19885","NG-6564_RIL72_lib19886","NG-6564_RIL73_lib19887","NG-6564_RIL74_lib19888","NG-6564_RIL75","NG-6564_RIL76_lib19890","NG-6564_RIL80_lib19891","NG-6564_RIL81_lib19892","NG-6564_RIL82_lib19893","NG-6564_RIL84_lib19894","NG-6564_RIL85_lib19895","NG-6564_RIL88_lib19896","NG-6564_RIL89_lib19897","NG-6564_RIL90_lib19898","NG-6564_RIL91_lib19899","NG-6564_RIL92_lib19900","NG-6564_RIL93_lib19901","NG-6564_RIL95_lib19902","NG-6564_RIL99_lib19903")
infiles_exp1=c("NG-6280_RIL138_P1","NG-6280_RIL83_P1","NG-6280_RIL94_P1","NG-6280_RIL96_P1","NG-6564_RIL100_lib19904","NG-6564_RIL103_lib19905","NG-6564_RIL104_lib19906","NG-6564_RIL105_lib19907","NG-6564_RIL106_lib19908","NG-6564_RIL107_lib19909","NG-6564_RIL109_lib19910","NG-6564_RIL111_lib19911","NG-6564_RIL113_lib19912","NG-6564_RIL114_lib19913","NG-6564_RIL115_lib19914","NG-6564_RIL116_lib19915","NG-6564_RIL117_lib19916","NG-6564_RIL118_lib19917","NG-6564_RIL119_lib19918","NG-6564_RIL11_lib19839","NG-6564_RIL120_lib19919","NG-6564_RIL121_lib19920","NG-6564_RIL122","NG-6564_RIL123_lib19922","NG-6564_RIL124_lib19923","NG-6564_RIL125_lib19924","NG-6564_RIL126_lib19925","NG-6564_RIL127_lib19926","NG-6564_RIL128_lib19927","NG-6564_RIL130_lib19928","NG-6564_RIL131_lib19929","NG-6564_RIL132_lib19930","NG-6564_RIL133_lib19931","NG-6564_RIL134_lib19932","NG-6564_RIL135_lib19933","NG-6564_RIL136_lib19934","NG-6564_RIL137_lib19935","NG-6564_RIL139_lib19936","NG-6564_RIL13_lib19840","NG-6564_RIL140_lib19937","NG-6564_RIL141_lib19938","NG-6564_RIL16_lib19841","NG-6564_RIL17_lib19842","NG-6564_RIL18_lib19843","NG-6564_RIL19_lib19844","NG-6564_RIL23_lib19845","NG-6564_RIL25_lib19846","NG-6564_RIL26_lib19847","NG-6564_RIL27_lib19848","NG-6564_RIL29_lib19849","NG-6564_RIL30_lib19850","NG-6564_RIL32_lib19851","NG-6564_RIL33_lib19852","NG-6564_RIL36_lib19853","NG-6564_RIL37_lib19854","NG-6564_RIL38_lib19855","NG-6564_RIL39_lib19856","NG-6564_RIL40_lib19857","NG-6564_RIL42_lib19858","NG-6564_RIL43_lib19859","NG-6564_RIL44_lib19860","NG-6564_RIL45_lib19861","NG-6564_RIL46_lib19862","NG-6564_RIL47_lib19863","NG-6564_RIL48_lib19864","NG-6564_RIL49_lib19865","NG-6564_RIL50_lib19866","NG-6564_RIL51_lib19867","NG-6564_RIL52_lib19868","NG-6564_RIL53_lib19869","NG-6564_RIL54_lib19870","NG-6564_RIL55_lib19871","NG-6564_RIL57_lib19872","NG-6564_RIL58_lib19873","NG-6564_RIL59_lib19874","NG-6564_RIL60_lib19875","NG-6564_RIL61_lib19876","NG-6564_RIL62_lib19877","NG-6564_RIL63_lib19878","NG-6564_RIL64_lib19879","NG-6564_RIL65_lib19880","NG-6564_RIL66_lib19881","NG-6564_RIL67_lib19882","NG-6564_RIL69_lib19883","NG-6564_RIL70_lib19884","NG-6564_RIL71_lib19885","NG-6564_RIL72_lib19886","NG-6564_RIL73_lib19887","NG-6564_RIL74_lib19888","NG-6564_RIL75","NG-6564_RIL76_lib19890","NG-6564_RIL80_lib19891","NG-6564_RIL81_lib19892","NG-6564_RIL82_lib19893","NG-6564_RIL84_lib19894","NG-6564_RIL85_lib19895","NG-6564_RIL88_lib19896","NG-6564_RIL89_lib19897","NG-6564_RIL90_lib19898","NG-6564_RIL91_lib19899","NG-6564_RIL92_lib19900","NG-6564_RIL93_lib19901","NG-6564_RIL95_lib19902","NG-6564_RIL99_lib19903")
infiles_exp2=c("NG-6860_R_1","NG-6860_R_10","NG-6860_R_11","NG-6860_R_12","NG-6860_R_13","NG-6860_R_14","NG-6860_R_15","NG-6860_R_16","NG-6860_R_18","NG-6860_R_19","NG-6860_R_2","NG-6860_R_20","NG-6860_R_22","NG-6860_R_23","NG-6860_R_24","NG-6860_R_26","NG-6860_R_27","NG-6860_R_28","NG-6860_R_29","NG-6860_R_3","NG-6860_R_30","NG-6860_R_31","NG-6860_R_32","NG-6860_R_33","NG-6860_R_34","NG-6860_R_35","NG-6860_R_36","NG-6860_R_37","NG-6860_R_38","NG-6860_R_39","NG-6860_R_4","NG-6860_R_40","NG-6860_R_41","NG-6860_R_42","NG-6860_R_44","NG-6860_R_45","NG-6860_R_46","NG-6860_R_47","NG-6860_R_48","NG-6860_R_49","NG-6860_R_5","NG-6860_R_50","NG-6860_R_51","NG-6860_R_52","NG-6860_R_6","NG-6860_R_7","NG-6860_R_8","NG-6860_R_9","NG-7084_R_53","NG-7084_R_54","NG-7084_R_55","NG-7084_R_57","NG-7084_R_58","NG-7084_R_59","NG-7084_R_60","NG-7084_R_61","NG-7084_R_62","NG-7084_R_63","NG-7084_R_64","NG-7084_R_65","NG-7084_R_66","NG-7084_R_67","NG-7084_R_68","NG-7084_R_69","NG-7084_R_70","NG-7084_R_71","NG-7084_R_72","NG-7084_R_73","NG-7084_R_74","NG-7084_R_75","NG-7084_R_76","NG-7084_R_77","NG-7084_R_78","NG-7084_R_79","NG-7084_R_80","NG-7084_R_81","NG-7084_R_82","NG-7084_R_83","NG-7084_R_84","NG-7084_R_85","NG-7084_R_86","NG-7084_R_87","NG-7084_R_88","NG-7084_R_89","NG-7084_R_90","NG-7084_R_91","NG-7084_R_92","NG-7084_R_93","NG-7084_R_94","NG-7084_R_95","NG-7084_R_96","NG-7084_R_97","NG-7084_R_98","NG-7084_R_99")
infiles=c(infiles_exp1, infiles_exp2)
workdir="/home/m.miculan/zea_mays_projects/eQTL/stringtie_output"
out_fpkm="/home/m.miculan/zea_mays_projects/eQTL/stringtie_output/ril.fpkm.chr.txt"
out_tpm="/home/m.miculan/zea_mays_projects/eQTL/stringtie_output/ril.tpm.chr.txt"

setwd(workdir)
genes=read.table("/home/m.miculan/zea_mays_projects/eQTL/genes.txt")
setnames(genes, "gene_id")

# FPKM
for(i in 1:length(infiles))
{
     cat(i,"\n")
     ril_file_name=paste0(workdir, "/", infiles[i], "/", infiles[i], ".StringTie.gene_abundance.tab")
     ril=fread(ril_file_name, sep="\t", data.table=F)
     ril=ril[ril$"Reference"!="Mt", ]
     ril=ril[ril$"Reference"!="Pt", ]
     ril=ril[order(ril[,1]), c(1,8)]
     setnames(ril,c("gene_id",infiles[i]))

     #ril=read.table(paste0(infiles[i],"/",infiles[i],".StringTie.gene_abundance.tab"),sep="\t", header=TRUE)
    
        if(i==1) {
            allril=merge(genes, ril, by="gene_id", all.x=T, sort=F)
#            allril=inner_join(genes, ril, by="gene_id")
                }
        if(i!=1) {
           allril=merge(allril, ril, by="gene_id", all.x=T, sort=F)
#            allril=inner_join(allril, ril, by="gene_id")
                   }    
 }
allril$"NG-6280_RIL138_P2"=NULL
allril$"NG-6280_RIL83_P2"=NULL
allril$"NG-6280_RIL94_P2"=NULL
allril$"NG-6280_RIL96_P2"=NULL
write.table(allril,out_fpkm,sep="\t",quote=FALSE,row.names=FALSE)


# TPM
for(i in 1:length(infiles))
{
     cat(i,"\n")
     ril_file_name=paste0(workdir, "/", infiles[i], "/", infiles[i], ".StringTie.gene_abundance.tab")
     ril=fread(ril_file_name, sep="\t", data.table=F)
     ril=ril[ril$"Reference"!="Mt", ]
     ril=ril[ril$"Reference"!="Pt", ]
     ril=ril[order(ril$"Gene ID"), c(1,9)]
     setnames(ril,c("gene_id",infiles[i]))

     #ril=read.table(paste0(infiles[i],"/",infiles[i],".StringTie.gene_abundance.tab"),sep="\t", header=TRUE)
    
        if(i==1) {
#            allril=merge(genes, ril, by="gene_id", all.x=T, sort=F) #maybe we do not need
#            allril=inner_join(genes, ril, by="gene_id")
             allril=ril
                }
        if(i!=1) {
           allril=merge(allril, ril, by="gene_id", all.x=T, sort=F)
#            allril=inner_join(allril, ril, by="gene_id")
                   }    
 }
allril$"NG-6280_RIL138_P2"=NULL
allril$"NG-6280_RIL83_P2"=NULL
allril$"NG-6280_RIL94_P2"=NULL
allril$"NG-6280_RIL96_P2"=NULL
write.table(allril,out_tpm,sep="\t",quote=FALSE,row.names=FALSE)


# PARENTS
#######################################################################################################################################################################################
R
library(data.table)
library(dplyr)

#infiles_parents=c("A632_repl1", "B73_repl1", "CML91_repl1", "F7_repl1", "H99_repl1", "HP301_repl1", "Mo17_repl1", "W153R_repl1")
infiles_parents=c("A632", "B73", "CML91", "F7", "H99", "HP301", "Mo17", "W153R")
workdir="/home/m.miculan/zea_mays_projects/eQTL/stringtie_output"
out_tpm_parents="/home/m.miculan/zea_mays_projects/eQTL/stringtie_output/parents.tpm.chr.txt"

setwd(workdir)
genes=read.table("/home/m.miculan/zea_mays_projects/eQTL/genes.txt")
setnames(genes, "gene_id")


# TPM
for(i in 1:length(infiles_parents))
{
     cat(i,"\n")
     parents_file_name=paste0(workdir, "/", infiles_parents[i], "/", infiles_parents[i], ".StringTie.gene_abundance.tab")
     parents=fread(parents_file_name, sep="\t", data.table=F)
     parents=parents[parents$"Reference"!="Mt", ]
     parents=parents[parents$"Reference"!="Pt", ]
     parents=parents[order(parents$"Gene ID"), c(1,9)]
     setnames(parents,c("gene_id",infiles_parents[i]))
     dim(parents)

     #parents=read.table(paste0(infiles[i],"/",infiles[i],".StringTie.gene_abundance.tab"),sep="\t", header=TRUE)
    
        if(i==1) {
             allparents=merge(genes, parents, by="gene_id", all.x=T, sort=F) #maybe we do not need
#            allparents=inner_join(genes, parents, by="gene_id")
             #allparents=parents
                }
        if(i!=1) {
           allparents=merge(allparents, parents, by="gene_id", all.x=T, sort=F)
#            allparents=inner_join(allparents, parents, by="gene_id")
                   }    
 }
write.table(allparents,out_tpm_parents,sep="\t",quote=FALSE,row.names=FALSE)



# DIVIDE BIPARENTAL POPULATION and MULTIPARENTAL POPULATION
#######################################################################################################################################################################################
infiles_parents=c("A632", "B73", "CML91", "F7", "H99", "HP301", "Mo17", "W153R")
tpm_parents_file="/home/m.miculan/zea_mays_projects/eQTL/stringtie_output/parents.tpm.chr.txt"
tpm_ril_file="/home/m.miculan/zea_mays_projects/eQTL/stringtie_output/ril.tpm.chr.txt"

infiles_bi=c("NG-6280_RIL138_P1","NG-6280_RIL83_P1","NG-6280_RIL94_P1","NG-6280_RIL96_P1","NG-6564_RIL100_lib19904","NG-6564_RIL103_lib19905","NG-6564_RIL104_lib19906","NG-6564_RIL105_lib19907","NG-6564_RIL106_lib19908","NG-6564_RIL107_lib19909","NG-6564_RIL109_lib19910","NG-6564_RIL111_lib19911","NG-6564_RIL113_lib19912","NG-6564_RIL114_lib19913","NG-6564_RIL115_lib19914","NG-6564_RIL116_lib19915","NG-6564_RIL117_lib19916","NG-6564_RIL118_lib19917","NG-6564_RIL119_lib19918","NG-6564_RIL11_lib19839","NG-6564_RIL120_lib19919","NG-6564_RIL121_lib19920","NG-6564_RIL122","NG-6564_RIL123_lib19922","NG-6564_RIL124_lib19923","NG-6564_RIL125_lib19924","NG-6564_RIL126_lib19925","NG-6564_RIL127_lib19926","NG-6564_RIL128_lib19927","NG-6564_RIL130_lib19928","NG-6564_RIL131_lib19929","NG-6564_RIL132_lib19930","NG-6564_RIL133_lib19931","NG-6564_RIL134_lib19932","NG-6564_RIL135_lib19933","NG-6564_RIL136_lib19934","NG-6564_RIL137_lib19935","NG-6564_RIL139_lib19936","NG-6564_RIL13_lib19840","NG-6564_RIL140_lib19937","NG-6564_RIL141_lib19938","NG-6564_RIL16_lib19841","NG-6564_RIL17_lib19842","NG-6564_RIL18_lib19843","NG-6564_RIL19_lib19844","NG-6564_RIL23_lib19845","NG-6564_RIL25_lib19846","NG-6564_RIL26_lib19847","NG-6564_RIL27_lib19848","NG-6564_RIL29_lib19849","NG-6564_RIL30_lib19850","NG-6564_RIL32_lib19851","NG-6564_RIL33_lib19852","NG-6564_RIL36_lib19853","NG-6564_RIL37_lib19854","NG-6564_RIL38_lib19855","NG-6564_RIL39_lib19856","NG-6564_RIL40_lib19857","NG-6564_RIL42_lib19858","NG-6564_RIL43_lib19859","NG-6564_RIL44_lib19860","NG-6564_RIL45_lib19861","NG-6564_RIL46_lib19862","NG-6564_RIL47_lib19863","NG-6564_RIL48_lib19864","NG-6564_RIL49_lib19865","NG-6564_RIL50_lib19866","NG-6564_RIL51_lib19867","NG-6564_RIL52_lib19868","NG-6564_RIL53_lib19869","NG-6564_RIL54_lib19870","NG-6564_RIL55_lib19871","NG-6564_RIL57_lib19872","NG-6564_RIL58_lib19873","NG-6564_RIL59_lib19874","NG-6564_RIL60_lib19875","NG-6564_RIL61_lib19876","NG-6564_RIL62_lib19877","NG-6564_RIL63_lib19878","NG-6564_RIL64_lib19879","NG-6564_RIL65_lib19880","NG-6564_RIL66_lib19881","NG-6564_RIL67_lib19882","NG-6564_RIL69_lib19883","NG-6564_RIL70_lib19884","NG-6564_RIL71_lib19885","NG-6564_RIL72_lib19886","NG-6564_RIL73_lib19887","NG-6564_RIL74_lib19888","NG-6564_RIL75","NG-6564_RIL76_lib19890","NG-6564_RIL80_lib19891","NG-6564_RIL81_lib19892","NG-6564_RIL82_lib19893","NG-6564_RIL84_lib19894","NG-6564_RIL85_lib19895","NG-6564_RIL88_lib19896","NG-6564_RIL89_lib19897","NG-6564_RIL90_lib19898","NG-6564_RIL91_lib19899","NG-6564_RIL92_lib19900","NG-6564_RIL93_lib19901","NG-6564_RIL95_lib19902","NG-6564_RIL99_lib19903")
output_bi_file="/home/m.miculan/zea_mays_projects/eQTL/stringtie_output/ril.tpm.chr.bi.txt"

infiles_multi=c("NG-6860_R_1","NG-6860_R_10","NG-6860_R_11","NG-6860_R_12","NG-6860_R_13","NG-6860_R_14","NG-6860_R_15","NG-6860_R_16","NG-6860_R_18","NG-6860_R_19","NG-6860_R_2","NG-6860_R_20","NG-6860_R_22","NG-6860_R_23","NG-6860_R_24","NG-6860_R_26","NG-6860_R_27","NG-6860_R_28","NG-6860_R_29","NG-6860_R_3","NG-6860_R_30","NG-6860_R_31","NG-6860_R_32","NG-6860_R_33","NG-6860_R_34","NG-6860_R_35","NG-6860_R_36","NG-6860_R_37","NG-6860_R_38","NG-6860_R_39","NG-6860_R_4","NG-6860_R_40","NG-6860_R_41","NG-6860_R_42","NG-6860_R_44","NG-6860_R_45","NG-6860_R_46","NG-6860_R_47","NG-6860_R_48","NG-6860_R_49","NG-6860_R_5","NG-6860_R_50","NG-6860_R_51","NG-6860_R_52","NG-6860_R_6","NG-6860_R_7","NG-6860_R_8","NG-6860_R_9","NG-7084_R_53","NG-7084_R_54","NG-7084_R_55","NG-7084_R_57","NG-7084_R_58","NG-7084_R_59","NG-7084_R_60","NG-7084_R_61","NG-7084_R_62","NG-7084_R_63","NG-7084_R_64","NG-7084_R_65","NG-7084_R_66","NG-7084_R_67","NG-7084_R_68","NG-7084_R_69","NG-7084_R_70","NG-7084_R_71","NG-7084_R_72","NG-7084_R_73","NG-7084_R_74","NG-7084_R_75","NG-7084_R_76","NG-7084_R_77","NG-7084_R_78","NG-7084_R_79","NG-7084_R_80","NG-7084_R_81","NG-7084_R_82","NG-7084_R_83","NG-7084_R_84","NG-7084_R_85","NG-7084_R_86","NG-7084_R_87","NG-7084_R_88","NG-7084_R_89","NG-7084_R_90","NG-7084_R_91","NG-7084_R_92","NG-7084_R_93","NG-7084_R_94","NG-7084_R_95","NG-7084_R_96","NG-7084_R_97","NG-7084_R_98","NG-7084_R_99")
output_multi_file="/home/m.miculan/zea_mays_projects/eQTL/stringtie_output/ril.tpm.chr.multi.txt" 

# read input files
tpm.parents=fread(tpm_parents_file ,sep="\t", data.table=F)
tpm.ril=fread(tpm_ril_file ,sep="\t", data.table=F)



# take BIPARENTAL
##################################################################################
#individuals=c("B73_repl1", "H99_repl1", infiles_bi )
individuals=c("B73", "H99", infiles_bi )
tpm.bi=merge(tpm.parents, tpm.ril, by="gene_id", all.x=T, sort=F)
tpm.bi=tpm.bi[ , c("gene_id", individuals)]
write.table(tpm.bi,output_bi_file,sep="\t",quote=FALSE,row.names=FALSE)



# take MULTIPARENTAL
##################################################################################
individuals=c(infiles_parents, infiles_multi )
tpm.multi=merge(tpm.parents, tpm.ril, by="gene_id", all.x=T, sort=F)
tpm.multi=tpm.multi[ , c("gene_id", individuals) ]
write.table(tpm.multi,output_multi_file,sep="\t",quote=FALSE,row.names=FALSE)


# put parents together!
#######################################################################################################################################################################################
# since we have replicates, we merge together raw counts, then we normalize


R
library(data.table)

setwd("/home/m.miculan/zea_mays_projects/eQTL/mappings")

for (i in c("A632", "B73", "CML91", "F7", "H99", "HP301", "Mo17", "W153R"))
{
    cat(i,"\n")
    
    output_name=paste0(i, "/", i, ".ReadsPerGene.out.tab")
    gene_abundance_file1=paste0("/home/m.miculan/zea_mays_projects/eQTL/mappings/", i, "_repl1/", i, "_repl1.ReadsPerGene.out.tab" )
    gene_abundance_file2=paste0("/home/m.miculan/zea_mays_projects/eQTL/mappings/", i, "_repl2/", i, "_repl2.ReadsPerGene.out.tab" )
    
    ab1=fread(gene_abundance_file1, sep="\t", data.table=F)
    ab1=ab1[ 5:nrow(ab1), 1:2]
    names(ab1)=c("gene_id","count_r1")

    ab2=fread(gene_abundance_file2, sep="\t", data.table=F)
    ab2=ab2[ 5:nrow(ab2), 1:2]
    names(ab2)=c("gene_id","count_r2")
    
    abM=merge(ab1, ab2)
    abM$"count"= abM$"count_r1" + abM$"count_r2"
    abM=abM[ , c("gene_id", "count")]
    write.table(abM,output_name,quote=F,sep="\t",row.names=F)
} 

quit()


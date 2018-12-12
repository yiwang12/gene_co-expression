#convert ensemble id to gene symbols using biomaRt package getBM function
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#genes <- c("ENSG00000012048","ENSG00000214179","ENSG00000200632")
#df=read.table("/Users/yiwang/Downloads/GD660.GeneQuantRPKM.txt",header=T)
df=read.table("/users/ywang/12_11/data/GD660.GeneQuantRPKM.txt",header=T)
dim(df)
# 53934   664

genes_original <- df$Gene_Symbol
list_gid=strsplit(as.character(genes_original),split=".", fixed=TRUE)
genes=unlist(list_gid)[ c(TRUE,FALSE) ]
df$ensembl_gene_id=genes
#df<-df[,-4]
#G_list <- getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id","hgnc_symbol"),values=genes,mart= mart)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)

merged_list = merge(df,G_list,by.x="ensembl_gene_id")
#,by.y="ensembl_gene_id")
save(merged_list,file="/users/ywang/12_11/out/merged_list_GD660_GeneQuantRPKM.RData")

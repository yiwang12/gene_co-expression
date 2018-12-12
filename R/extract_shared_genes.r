#extract shared genes in UMI-scRNA and bulk RNA-seq data
df_UMI=read.table("/users/ywang/12_11/data/GSM3044892_GeneExp.UMIs.10X2.txt",header=T)#load UMI data
dim(df_UMI)
#9474 6484
#df_UMI=read.table("/Users/yiwang/Dropbox/project_Kasper/12_11/data/GSM3044892_GeneExp.UMIs.10X2_small.txt",header=T)
#df_UMI$genes
#ATRX  TCOF1 NSRP1 SPPL3 OPA3  OPA1  ITGA4
df_UMI$hgnc_symbol=df_UMI$genes # add one column hgnc_symbol to UMI count data

load("/users/ywang/12_11/data/merged_list_GD660_GeneQuantRPKM.RData")#load merged_list of bulk data
merged_list_UMI_bulk = merge(df_UMI,merged_list,by.x="hgnc_symbol") #merge bulk and UMI data together
save(merged_list_UMI_bulk,file="/users/ywang/12_11/out/merged_list_UMI_bulk_10X2_GD660_GeneQuantRPKM.RData")




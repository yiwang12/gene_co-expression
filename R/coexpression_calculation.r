#### calculate co-expression in UMI and bulk data separately

#### filt high expression genes
load("/users/ywang/12_11/out/merged_list_UMI_bulk_10X2_GD660_GeneQuantRPKM.RData") ## load merged bulk and UMI data : merged_list_UMI_bulk
#> dim(merged_list_UMI_bulk)
#[1] 8781 7150
df_UMI_subsample = read.table("/users/ywang/12_11/data/GSM3044892_GeneExp.UMIs.10X2_small.txt",header=T)#load UMI data
#df_UMI_subsample = read.table("/Users/yiwang/Dropbox/project_Kasper/12_11/data/GSM3044892_GeneExp.UMIs.10X2_small.txt",header=T)#load UMI data

col_numbers_UMI = dim(df_UMI_subsample)[2]+1
UMI_count = merged_list_UMI_bulk[, 3:col_numbers_UMI]
col_numbers_total = dim(merged_list_UMI_bulk)[2]
bulk_count = merged_list_UMI_bulk[, (col_numbers_UMI+1+5):col_numbers_total]
save(bulk_count,file="/users/ywang/12_11/out/bulk_count.RData")
save(UMI_count,file="/users/ywang/12_11/out/UMI_count.RData")

## normalization of UMI data bu cpm
#lib.size_all=colSums(UMI_count)
#UMI_count_norm = (UMI_count+1)*1000000/(lib.size_all+1)
#save(UMI_count_norm,file="/users/ywang/12_11/out/UMI_count_norm.RData")

#####################################
mean_bulk = rowMeans(bulk_count)
mean_UMI = rowMeans(UMI_count_norm)

order_mean_bulk = order(mean_bulk)
order_mean_UMI = order(mean_UMI)
order_overall = order(-(order_mean_bulk + order_mean_UMI)) # order of overall expression levels, from high to low

save(order_overall,file="/users/ywang/12_11/out/order_overall.RData")
#select overall high expression 1000 genes
bulk_count_filt = bulk_count[order_overall[1:1000],] 
#UMI_count_norm_filt = UMI_count_norm[order_overall[1:1000],]
merged_list_UMI_bulk_filt = merged_list_UMI_bulk[order_overall[1:1000],]
UMI_count_raw_filt = UMI_count[order_overall[1:1000],]

save(UMI_count_raw_filt,file="/users/ywang/12_11/out/UMI_count_raw_filt.RData")
save(bulk_count_filt,file="/users/ywang/12_11/out/bulk_count_filt.RData")
save(merged_list_UMI_bulk_filt,file="/users/ywang/12_11/out/merged_list_UMI_bulk_filt.RData")

## normalization of UMI data bu cpm
lib.size_all=colSums(UMI_count_raw_filt)
UMI_count_norm_filt = (UMI_count_raw_filt+1)*1000000/(lib.size_all+1)
save(UMI_count_norm_filt,file="/users/ywang/12_11/out/UMI_count_norm_filt.RData")


##### correlation calculation!
load("/users/ywang/12_11/out/bulk_count_filt.RData")
bulk_count_filt_log2=log2(bulk_count_filt+1)#log2 (RPKM+1)
save(bulk_count_filt_log2,file="/users/ywang/12_11/out/bulk_count_filt_log2.RData")
cor_bulk = cor(t(bulk_count_filt_log2),method = c("pearson"))
save(cor_bulk,file="/users/ywang/12_11/out/cor_bulk.RData")

#> dim(cor_bulk)
#[1] 1000 1000

load("/users/ywang/12_11/out/UMI_count_norm_filt.RData")
UMI_count_norm_filt_log2 = log2(UMI_count_norm_filt)
save(UMI_count_norm_filt_log2,file="/users/ywang/12_11/out/UMI_count_norm_filt_log2.RData")

cor_UMI = cor(t(UMI_count_norm_filt_log2),method = c("pearson"))
#> dim(cor_UMI)
#[1] 1000 1000
save(cor_UMI,file="/users/ywang/12_11/out/cor_UMI.RData")

##### correlation comparison
gene_number=1000
cor_UMI_vector=array(dim=c(gene_number*(gene_number-1)/2))
cor_bulk_vector=array(dim=c(gene_number*(gene_number-1)/2))

count_tmp=0
for(i in 1:(gene_number-1)){
	for(j in (i+1):gene_number){
		count_tmp=count_tmp+1
		cor_UMI_vector[count_tmp]=cor_UMI[i,j]
		cor_bulk_vector[count_tmp]=cor_bulk[i,j]
	}
}
save(cor_UMI_vector,file="/users/ywang/12_11/out/cor_UMI_vector.RData")
save(cor_bulk_vector,file="/users/ywang/12_11/out/cor_bulk_vector.RData")

#### run in PC
load("/Users/yiwang/Dropbox/project_Kasper/12_11/out/cor_bulk_vector.RData")
load("/Users/yiwang/Dropbox/project_Kasper/12_11/out/cor_UMI_vector.RData")

cor_UMI_bulk=round(cor(cor_bulk_vector,cor_UMI_vector),2)
png("/Users/yiwang/Dropbox/project_Kasper/12_11/figures/cor_UMI_bulk.png")
plot(cor_bulk_vector,cor_UMI_vector,xlab="bulk RNA-seq log2(RPKM+1)",ylab="scRNA-seq log2(CPM)",cex.lab=1.3, cex.axis=1.5)
lines(lowess(cor_bulk_vector,cor_UMI_vector),col="red")
mtext(paste0("R = ",cor_UMI_bulk),side=3,cex = 1.5)
dev.off()
#cor.test(cor_bulk_vector,cor_UMI_vector)
#p-value < 2.2e-16
#> cor(cor_bulk_vector,cor_UMI_vector)
#[1] 0.1791646




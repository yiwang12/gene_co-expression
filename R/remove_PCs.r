load("/users/ywang/12_11/out/UMI_count_raw_filt_saver_norm_log2.RData")#load saver recovered UMI data(after normalize):

library(WGCNA)
UMI_saver_calibrated=removePrincipalComponents(UMI_count_raw_filt_saver_norm_log2,4)#remove top 4 PCs
save(UMI_saver_calibrated,file="/users/ywang/12_11/out/UMI_saver_calibrated.RData")

load("/users/ywang/12_11/out/bulk_count_filt_log2.RData")#remove top 4 PCs
bulk_count_filt_log2=as.matrix(bulk_count_filt_log2)
bulk_calibrated=removePrincipalComponents(bulk_count_filt_log2,4)
save(bulk_calibrated,file="/users/ywang/12_11/out/bulk_calibrated.RData")

load("/users/ywang/12_11/out/UMI_count_norm_filt_log2.RData")#remove top 4 PCs
UMI_count_norm_filt_log2=as.matrix(UMI_count_norm_filt_log2)
UMI_original_calibrated=removePrincipalComponents(UMI_count_norm_filt_log2,4)
save(UMI_original_calibrated,file="/users/ywang/12_11/out/UMI_original_calibrated.RData")

###############
cor_bulk_calibrated = cor(t(bulk_calibrated),method = c("pearson"))
cor_UMI_original_calibrated = cor(t(UMI_original_calibrated),method = c("pearson"))
cor_UMI_saver_calibrated = cor(t(UMI_saver_calibrated),method = c("pearson"))

save(cor_UMI_original_calibrated,file="/users/ywang/12_11/out/cor_UMI_original_calibrated.RData")
save(cor_UMI_saver_calibrated,file="/users/ywang/12_11/out/cor_UMI_saver_calibrated.RData")
save(cor_bulk_calibrated,file="/users/ywang/12_11/out/cor_bulk_calibrated.RData")

###############
gene_number=1000
cor_UMI_vector_calibrated=array(dim=c(gene_number*(gene_number-1)/2))
cor_bulk_vector_calibrated=array(dim=c(gene_number*(gene_number-1)/2))
cor_UMI_saver_vector_calibrated=array(dim=c(gene_number*(gene_number-1)/2))

count_tmp=0
for(i in 1:(gene_number-1)){
	for(j in (i+1):gene_number){
		count_tmp=count_tmp+1
		cor_UMI_vector_calibrated[count_tmp]=cor_UMI_original_calibrated[i,j]
		cor_bulk_vector_calibrated[count_tmp]=cor_bulk_calibrated[i,j]
		cor_UMI_saver_vector_calibrated[count_tmp]=cor_UMI_saver_calibrated[i,j]
	}
}
save(cor_UMI_vector_calibrated,file="/users/ywang/12_11/out/cor_UMI_vector_calibrated.RData")
save(cor_bulk_vector_calibrated,file="/users/ywang/12_11/out/cor_bulk_vector_calibrated.RData")
save(cor_UMI_saver_vector_calibrated,file="/users/ywang/12_11/out/cor_UMI_saver_vector_calibrated.RData")

###############run on PC
load("/Users/yiwang/Dropbox/project_Kasper/12_11/out/cor_UMI_vector_calibrated.RData")
load("/Users/yiwang/Dropbox/project_Kasper/12_11/out/cor_bulk_vector_calibrated.RData")
cor_UMI_bulk=round(cor(cor_bulk_vector_calibrated,cor_UMI_vector_calibrated),2)
png("/Users/yiwang/Dropbox/project_Kasper/12_11/figures/cor_UMI_bulk_calibrated.png")
plot(cor_bulk_vector_calibrated,cor_UMI_vector_calibrated,xlab="bulk RNA-seq log2(RPKM+1)",ylab="scRNA-seq log2(CPM)",cex.lab=1.3, cex.axis=1.5)
lines(lowess(cor_bulk_vector_calibrated,cor_UMI_vector_calibrated),col="red")
mtext(paste0("R = ",cor_UMI_bulk),side=3,cex = 1.5)
dev.off()

load("/Users/yiwang/Dropbox/project_Kasper/12_11/out/cor_UMI_saver_vector_calibrated.RData")
cor_UMI_saver_bulk=round(cor(cor_bulk_vector_calibrated,cor_UMI_saver_vector_calibrated),2)
png("/Users/yiwang/Dropbox/project_Kasper/12_11/figures/cor_UMI_saver_bulk_calibrated.png")
plot(cor_bulk_vector_calibrated,cor_UMI_saver_vector_calibrated,xlab="bulk RNA-seq log2(RPKM+1)",ylab="scRNA-seq log2(CPM) - SAVER - cor",cex.lab=1.3, cex.axis=1.5)
lines(lowess(cor_bulk_vector_calibrated,cor_UMI_saver_vector_calibrated),col="red")
mtext(paste0("R = ",cor_UMI_saver_bulk),side=3,cex = 1.5)
dev.off()

#load("/Users/yiwang/Dropbox/project_Kasper/12_11/out/cor_UMI_vector_calibrated.RData")#load correlation matrix of bulk data
cor_UMI_saver_original=round(cor(cor_UMI_vector_calibrated,cor_UMI_saver_vector_calibrated),2)
png("/Users/yiwang/Dropbox/project_Kasper/12_11/figures/cor_UMI_saver_original_calibrated.png")
plot(cor_UMI_vector_calibrated,cor_UMI_saver_vector_calibrated,xlab="scRNA-seq log2(CPM)",ylab="scRNA-seq log2(CPM) - SAVER - cor",cex.lab=1.3, cex.axis=1.5)
lines(lowess(cor_UMI_vector_calibrated,cor_UMI_saver_vector_calibrated),col="red")
mtext(paste0("R = ",cor_UMI_saver_original),side=3,cex = 1.5)
dev.off()

######
png("/Users/yiwang/Dropbox/project_Kasper/12_11/figures/correlation_distribution_calibrated.png")
d <- density(cor_UMI_vector_calibrated) # returns the density data 
plot(d,cex.lab=1.5, cex.axis=1.5,xlim=c(-1,1),col="blue",xlab="correlation",ylab="density",main="")#,main="Distribution of gene correlation")
d <- density(cor_bulk_vector_calibrated) # returns the density data 
lines(d,col="black")
d <- density(cor_UMI_saver_vector_calibrated) # returns the density data 
lines(d,col="red")
#mtext(paste0("R = ",cor_UMI_saver_original),side=3,cex = 1.5)
#legend("topleft", legend=c("bulk RNA", "UMI scRNA"),
#       col=c("black", "blue"), lty=1, cex=1.2)
legend("topleft", legend=c("bulk RNA", "scRNA","scRNA-SAVER"),
       col=c("black", "blue","red"), lty=1, cex=1.2)
dev.off()



> round(mean(abs(cor_UMI_saver_vector_calibrated)),2)
[1] 0.37
> round(mean(abs(cor_UMI_vector_calibrated)),2)
[1] 0.15
> round(mean(abs(cor_bulk_vector_calibrated)),2)
[1] 0.22








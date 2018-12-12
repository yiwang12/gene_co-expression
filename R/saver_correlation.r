###  calculate saver-recovered correlation and compare it with bulk data/original UMI data correlation
### nomalized saver-recovered UMI data
load("/users/ywang/12_11/out/UMI_count_raw_filt_saver_filt.RData")#load saver recovered raw UMI data:UMI_count_raw_filt_saver
#load("/users/ywang/12_11/out/UMI_count_raw_saver.RData")

lib.size_all=colSums(UMI_count_raw_filt_saver$estimate)
dim(UMI_count_raw_filt_saver$estimate)
#1000 6483
UMI_count_raw_filt_saver_norm = (UMI_count_raw_filt_saver$estimate+1)*1000000/(lib.size_all+1)  #cpm normalize 
save(UMI_count_raw_filt_saver_norm,file="/users/ywang/12_11/out/UMI_count_raw_filt_saver_norm.RData")
##################### correlation comparison between bulk and UMI-saver-recovered data
#load("/users/ywang/12_11/out/UMI_count_raw_filt_saver_norm")
UMI_count_raw_filt_saver_norm_log2 = log2(UMI_count_raw_filt_saver_norm)
save(UMI_count_raw_filt_saver_norm_log2,file="/users/ywang/12_11/out/UMI_count_raw_filt_saver_norm_log2.RData")

cor_UMI_saver = cor(t(UMI_count_raw_filt_saver_norm_log2),method = c("pearson"))
save(cor_UMI_saver,file="/users/ywang/12_11/out/cor_UMI_saver.RData")

### correlation calculation using SAVER package
library(SAVER)
cor_UMI_saver_package = cor.genes(UMI_count_raw_filt_saver)
save(cor_UMI_saver_package,file="/users/ywang/12_11/out/cor_UMI_saver_package.RData")


####### convert correlation matrix to vector
gene_number=1000 #1000 genes
cor_UMI_vector_saver=array(dim=c(gene_number*(gene_number-1)/2)) 
cor_UMI_vector_saver_package=array(dim=c(gene_number*(gene_number-1)/2)) 

count_tmp=0 
for(i in 1:(gene_number-1)){
	for(j in (i+1):gene_number){
		count_tmp=count_tmp+1
		cor_UMI_vector_saver[count_tmp]=cor_UMI_saver[i,j]
		cor_UMI_vector_saver_package[count_tmp]=cor_UMI_saver_package[i,j]
	}
}
save(cor_UMI_vector_saver,file="/users/ywang/12_11/out/cor_UMI_vector_saver.RData")
save(cor_UMI_vector_saver_package,file="/users/ywang/12_11/out/cor_UMI_vector_saver_package.RData")

###run on PC
load("/Users/yiwang/Dropbox/project_Kasper/12_11/out/cor_bulk_vector.RData")#load correlation matrix of bulk data
load("/Users/yiwang/Dropbox/project_Kasper/12_11/out/cor_UMI_vector_saver.RData")#load correlation matrix of UMI-saver-recovered data
load("/Users/yiwang/Dropbox/project_Kasper/12_11/out/cor_UMI_vector_saver_package.RData")#load correlation matrix of UMI-saver-recovered data

cor_UMI_saver_bulk=round(cor(cor_bulk_vector,cor_UMI_vector_saver_package),2)
png("/Users/yiwang/Dropbox/project_Kasper/12_11/figures/cor_UMI_saver_bulk.png")
plot(cor_bulk_vector,cor_UMI_vector_saver_package,xlab="bulk RNA-seq log2(RPKM+1)",ylab="scRNA-seq log2(CPM) - SAVER ",cex.lab=1.3, cex.axis=1.5)
lines(lowess(cor_bulk_vector,cor_UMI_vector_saver_package),col="red")
mtext(paste0("R = ",cor_UMI_saver_bulk),side=3,cex = 1.5)
dev.off()

##################### correlation comparison between UMI originate and UMI-saver-recovered data
###run on PC
load("/Users/yiwang/Dropbox/project_Kasper/12_11/out/cor_UMI_vector.RData")#load correlation matrix of bulk data
load("/Users/yiwang/Dropbox/project_Kasper/12_11/out/cor_UMI_vector_saver.RData")#load correlation matrix of bulk data
cor_UMI_saver_original=round(cor(cor_UMI_vector,cor_UMI_vector_saver_package),2)
png("/Users/yiwang/Dropbox/project_Kasper/12_11/figures/cor_UMI_saver_original.png")
plot(cor_UMI_vector,cor_UMI_vector_saver_package,xlab="UMI scRNA-seq log2(CPM+1)",ylab="scRNA-seq log2(CPM) - SAVER",cex.lab=1.3, cex.axis=1.5)
lines(lowess(cor_UMI_vector,cor_UMI_vector_saver_package),col="red")
mtext(paste0("R = ",cor_UMI_saver_original),side=3,cex = 1.5)
dev.off()

cor_UMI_saver_original=round(cor(cor_UMI_vector_saver_package,cor_UMI_vector_saver),2)
png("/Users/yiwang/Dropbox/project_Kasper/12_11/figures/cor_UMI_saver_package_saver_cor.png")
plot(cor_UMI_vector_saver_package,cor_UMI_vector_saver,xlab="scRNA-seq log2(CPM) - SAVER - cor.genes",ylab="scRNA-seq log2(CPM) - SAVER - cor",cex.lab=1.3, cex.axis=1.5)
#lines(lowess(cor_UMI_vector_saver_package,cor_UMI_vector_saver),col="red")
#mtext(paste0("R = ",cor_UMI_saver_original),side=3,cex = 1.5)
abline(0,1,col="skyblue")
dev.off()



############################ plot the distribution of correlation
load("/Users/yiwang/Dropbox/project_Kasper/12_11/out/cor_bulk_vector.RData")#load correlation matrix of bulk data
png("/Users/yiwang/Dropbox/project_Kasper/12_11/figures/correlation_distribution.png")
d <- density(cor_UMI_vector) # returns the density data 
plot(d,cex.lab=1.5, cex.axis=1.5,xlim=c(-1,1),col="blue",xlab="correlation",ylab="density",main="")#,main="Distribution of gene correlation")
d <- density(cor_bulk_vector) # returns the density data 
lines(d,col="black")
d <- density(cor_UMI_vector_saver_package) # returns the density data 
lines(d,col="red")
#mtext(paste0("R = ",cor_UMI_saver_original),side=3,cex = 1.5)
legend("topleft", legend=c("bulk RNA", "scRNA","scRNA-SAVER"),
       col=c("black", "blue","red"), lty=1, cex=1.2)
dev.off()



> round(mean(abs(cor_UMI_vector_saver)),2)
[1] 0.28
> round(mean(abs(cor_UMI_vector)),2)
[1] 0.12
> round(mean(abs(cor_bulk_vector)),2)
[1] 0.15


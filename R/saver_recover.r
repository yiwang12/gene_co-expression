library(devtools)
library(SAVER)
#install.packages("SAVER")

#load("/users/ywang/12_11/out/bulk_count_filt.RData")#load bulk RNA RPKM counts data

### run saver to recover UMI counts 
#load("/users/ywang/12_11/out/UMI_count_raw_filt.RData") #load UMI raw counts data
load("/users/ywang/12_11/out/UMI_count.RData")
UMI_count_raw_saver <- saver(as.matrix(UMI_count), ncores = 40) #run saver
save(UMI_count_raw_saver,file="/users/ywang/12_11/out/UMI_count_raw_saver.RData")


load("/users/ywang/12_11/out/order_overall.RData")#order_overall
#select overall high expression 1000 genes
UMI_count_raw_filt_saver = UMI_count_raw_saver[order_overall[1:1000],] 
save(UMI_count_raw_filt_saver,file="/users/ywang/12_11/out/UMI_count_raw_filt_saver.RData")




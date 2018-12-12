library(devtools)
library(SAVER)
#install.packages("SAVER")

#load("/users/ywang/12_11/out/bulk_count_filt.RData")#load bulk RNA RPKM counts data

### run saver to recover UMI counts 
load("/users/ywang/12_11/out/UMI_count_raw_filt.RData") #load UMI raw counts data
#load("/users/ywang/12_11/out/UMI_count.RData")
UMI_count_raw_filt_saver <- saver(as.matrix(UMI_count_raw_filt), ncores = 20) #run saver
save(UMI_count_raw_filt_saver,file="/users/ywang/12_11/out/UMI_count_raw_filt_saver_filt.RData")





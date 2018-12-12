load("/users/ywang/12_11/out/UMI_count_raw_filt_saver_norm_log2.RData")#load saver-recovered UMI data :UMI_count_norm_filt_saver


cor_UMI_saver = cor(t(UMI_count_raw_filt_saver_norm_log2),method = c("pearson"))# calculate correlation matrix of saver-recovered UMI data
save(cor_UMI_saver,file="/users/ywang/12_11/out/cor_UMI_saver.RData")


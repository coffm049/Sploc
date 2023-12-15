setwd("/panfs/jay/groups/4/miran045/shared/projects/BWAS_toolbox/experiments/Burden_Sploc/Results")


pvec_type1_sploc_perm=c()
pvec_type1_burden_perm=c()
for(i in 1:10){
  data=read.csv(paste0("Sploc_Burden_null_pvals_i10_j10_r10_es1_p0_",i,".csv"))
  data=as.data.frame(data)
  pvec_type1_sploc_perm=c(pvec_type1_sploc_perm,data$pvec_type1_sploc_perm)
  pvec_type1_burden_perm=c(pvec_type1_burden_perm,data$pvec_type1_burden_perm)
}

type1err_sploc_perm=sum(ifelse(pvec_type1_sploc_perm<0.05,1,0))/length(pvec_type1_sploc_perm)
type1err_burden_perm=sum(ifelse(pvec_type1_burden_perm<0.05,1,0))/length(pvec_type1_burden_perm)

#########################################################################################################
pvec_power_sploc_perm=c()
pvec_power_burden_perm=c()
for(i in 1:10){
  data=read.csv(paste0("Sploc_Burden_power_binones_i50_j10_r10_es1_p_0.5_",i,".csv"))
  data=as.data.frame(data)
  pvec_power_sploc_perm=c(pvec_power_sploc_perm,data$pvec_power_sploc_perm)
  pvec_power_burden_perm=c(pvec_power_burden_perm,data$pvec_power_burden_perm)
}

powererr_sploc_perm=sum(ifelse(pvec_power_sploc_perm<0.05,1,0))/length(pvec_power_sploc_perm)
powererr_burden_perm=sum(ifelse(pvec_power_burden_perm<0.05,1,0))/length(pvec_power_burden_perm)


#####################################################################################

pvec_power_sploc_perm_binones=c()
pvec_power_burden_perm_binones=c()
for(i in c(1:3,6:10)){
  data=read.csv(paste0('Sploc_Burden_alt_pvals_binones_1_',i,'.csv'))
  data=as.data.frame(data)
  pvec_power_sploc_perm_binones=c(pvec_power_sploc_perm_binones,data$pvec_power_sploc_perm)
  pvec_power_burden_perm_binones=c(pvec_power_burden_perm_binones,data$pvec_power_burden_perm)
}

powererr_sploc_perm_binones=sum(ifelse(pvec_power_sploc_perm_binones<0.05,1,0))/length(pvec_power_sploc_perm_binones)
powererr_burden_perm_binones=sum(ifelse(pvec_power_burden_perm_binones<0.05,1,0))/length(pvec_power_burden_perm_binones)

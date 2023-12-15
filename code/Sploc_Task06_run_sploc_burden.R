job_num = commandArgs(trailingOnly=TRUE)[1]
setwd("/panfs/jay/groups/4/miran045/shared/projects/BWAS_toolbox/experiments/Burden_Sploc")
source("Sploc_Task03_permutation_test.R")
source("Sploc_Task02_burden_test.R")
source("Sploc_Task05_Sploc_test.R")
#data=read.csv("data_Sploc_SIM.csv")

n_iter= 10
flag=c()


pvec_type1_sploc_perm=c()
pvec_type1_burden_perm=c()

pvec_power_sploc_perm=c()
pvec_power_burden_perm=c()
set.seed(job_num)

for (iter in 1:n_iter){
  source("Sploc_Task01_create_data.R")
  source("functional_dist.R")
  # estimate \sigma and V from null
  model <- lmer(Y_sim~(1|Site/Family),data=data,REML=TRUE)
  #summary(model)
  
  sigma_gamma_hat=as.data.frame(VarCorr(model))$sdcor[2]
  sigma_tau_hat=as.data.frame(VarCorr(model))$sdcor[1]
  sigma_err_hat=as.data.frame(VarCorr(model))$sdcor[3]
  fixed_eff= coef(summary(model))[1]
  # age_eff_hat= coef(summary(model))[2]
  # sex_eff_hat= coef(summary(model))[3]
  
  ## since there is no sex and age effect now, we set them to zero
  age_eff=0
  sex_eff=0
  age_eff_hat= NA
  sex_eff_hat= NA
  
  hat_vec=c(sigma_gamma_hat,sigma_tau_hat,sigma_err_hat,age_eff_hat,sex_eff_hat)#(size effect est,family efft est, err est, intercept est)
  true_vec=c(sigma_gamma,sigma_tau,sigma_err,age_eff,sex_eff)
  
  p=1000  # Number of Permutations for the test
  # pval=Sploc_test(data,fixed_eff,sigma_err_hat,sigma_gamma_hat,sigma_tau_hat,fam_size,p)
  # pvec_power=c(pvec_power,pval)
  # if (pval>0.05){
  #   flag_power=c(flag_power,0) # we fail to reject null
  #   
  # }else {
  #   flag_power=c(flag_power,1) # we reject the null
  # }
  pval_sploc=Sploc_test(data,fixed_eff,sigma_err_hat,sigma_gamma_hat,sigma_tau_hat,fam_size,p)
  
  if (prob_beta==0){
    pvec_type1_sploc_perm=c(pvec_type1_sploc_perm,pval_sploc$sploc_pval)
    pvec_type1_burden_perm=c(pvec_type1_burden_perm,pval_sploc$burden_pval)
    cat("Iteration ", iter, " complete.\n")
    type1_sploc_perm=sum(ifelse(pvec_type1_sploc_perm<0.05,1,0))/(iter)
    type1_burden_perm=sum(ifelse(pvec_type1_burden_perm<0.05,1,0))/(iter)
    
    pvec_type1=cbind(pvec_type1_sploc_perm,pvec_type1_burden_perm) %>% data.frame()
    type1_burd_sploc= data.frame(test_name=c("Sploc_perm","Burden_perm"),type_1_error=c(type1_sploc_perm,type1_burden_perm))
  #   write.csv(type1_burd_sploc,file=paste0('Sploc_Burden_type1_error_i',I,'_j',J,'_r',r,'_es',es,'_p0_',job_num,'.csv'))
  #   write.csv(pvec_type1,file=paste0('Sploc_Burden_null_pvals_i',I,'_j',J,'_r',r,'_es',es,'_p0_',job_num,'.csv')) 
  }else{
    pvec_power_sploc_perm=c(pvec_power_sploc_perm,pval_sploc$sploc_pval)
    pvec_power_burden_perm=c(pvec_power_burden_perm,pval_sploc$burden_pval)
    cat("Iteration ", iter, " complete.\n")
    power_sploc_perm=sum(ifelse(pvec_power_sploc_perm<0.05,1,0))/(iter)
    power_burden_perm=sum(ifelse(pvec_power_burden_perm<0.05,1,0))/(iter)
    
    pvec_power=cbind(pvec_power_sploc_perm,pvec_power_burden_perm) %>% data.frame()
    power_burd_sploc= data.frame(test_name=c("Sploc_perm","Burden_perm"),power=c(power_sploc_perm,power_burden_perm))
    # write.csv(power_burd_sploc,file=paste0('Sploc_Burden_power_allones_',job_num,'.csv'))
    # write.csv(pvec_power,file=paste0('Sploc_Burden_alt_pvals_allones_',job_num,'.csv'))   
  }
}

if (prob_beta==0){
  
  # pvec_type1_sploc_perm=c(pvec_type1_sploc_perm,pval_sploc[[1]][2])
  # pvec_type1_burden_perm=c(pvec_type1_burden_perm,pval_sploc[[1]][1])
  # cat("Iteration ", iter, " complete.\n")
  # type1_sploc_perm=sum(ifelse(pvec_type1_sploc_perm<0.05,1,0))/(iter)
  # type1_burden_perm=sum(ifelse(pvec_type1_burden_perm<0.05,1,0))/(iter)
  # 
  # pvec_type1=cbind(pvec_type1_sploc_perm,pvec_type1_burden_perm) %>% data.frame()
  # type1_burd_sploc= data.frame(test_name=c("Sploc_perm","Burden_perm"),type_1_error=c(type1_sploc_perm,type1_burden_perm))
  setwd("/panfs/jay/groups/4/miran045/shared/projects/BWAS_toolbox/experiments/Burden_Sploc/Results")
  write.csv(type1_burd_sploc,file=paste0('Sploc_Burden_type1_error_i',I,'_j',J,'_r',r,'_es',es,'_p0_',job_num,'.csv'))
  write.csv(pvec_type1,file=paste0('Sploc_Burden_null_pvals_i',I,'_j',J,'_r',r,'_es',es,'_p0_',job_num,'.csv'))
  
}else {
  
  setwd("/panfs/jay/groups/4/miran045/shared/projects/BWAS_toolbox/experiments/Burden_Sploc/Results")
  write.csv(power_burd_sploc,file=paste0('Sploc_Burden_power_i',I,'_j',J,'_r',r,'_es',es,'_p',prob_beta,'_',job_num,'.csv'))
  write.csv(pvec_power,file=paste0('Sploc_Burden_alt_pvals_i',I,'_j',J,'_r',r,'_es',es,'_p',prob_beta,'_',job_num,'.csv'))  
}
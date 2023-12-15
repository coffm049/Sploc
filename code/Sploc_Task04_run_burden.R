source("Sploc_Task03_permutation_test.R")
source("Sploc_Task02_burden_test.R")
#data=read.csv("data_Sploc_SIM.csv")

n_iter= 10
flag=c()

pvec_type2_burd=c()
pvec_type2_sploc=c()

set.seed(125)

for (iter in 1:n_iter){
  source("Sploc_Task01_create_data.R")
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
 
  p=10 # Number of Permutations for the test
  # pval=Burden_test(data,fixed_eff,sigma_err_hat,sigma_gamma_hat,sigma_tau_hat,fam_size,p)
  # pvec_type2=c(pvec_type2,pval)
  # if (pval>0.05){
  #   flag_type2=c(flag_type2,0) # we fail to reject null
  #   
  # }else {
  #   flag_type2=c(flag_type2,1) # we reject the null
  # }
  pval_burden=Burden_test(data,fixed_eff,sigma_err_hat,sigma_gamma_hat,sigma_tau_hat,fam_size,p)
  pvec_power_burd=c(pvec_type2_burd,pval_burden$pval_perm)
  #pvec_type2_sploc=c(pvec_type2_sploc,pvals$sploc_pval)
  cat("Iteration ", iter, " complete.\n")
}

power_burd=sum(ifelse(pvec_type2_burd<0.05,1,0))/(iter-1)
#power_sploc=sum(ifelse(pvec_type2_sploc<0.05,1,0))/(iter-1)
write.csv(power_burd,file='Burden_test_nocov_power.csv')
  
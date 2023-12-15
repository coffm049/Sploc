permutation_test<-function(data,fixed_eff,age_eff_hat,sex_eff_hat,sigma_err_hat,sigma_gamma_hat,sigma_tau_hat,fam_size,p,V_inv,n){
  
  
  Ustat_mat=matrix(0,n,p)
  
  
  for( j in 1:p){
    
    # to shuffle all the rows
    #data$Y_sim= sample(data$Y_sim, length(data$Y_sim),replace=FALSE)
    
    # #to shuffle within sites
    # uniq_y=unique(data$Site)
    # for(temp in 1:length(uniq_y)){
    #   subset <- data[data$Site == uniq_y[temp], ]
    #   permuted_indices <- sample(nrow(subset))
    #   subset$Y_sim= subset[permuted_indices, ]$Y_sim
    #   # Use the permuted indices to reorder the rows of the data frame
    #   data[data$Site == uniq_y[temp], ] <- subset
    #   
    # }
    
    
    #to shuffle within family
    data$Fam_idtfr=as.factor(data$Fam_idtfr)
    uniq_y=unique(data$Fam_idtfr)
    for(temp in 1:length(uniq_y)){
      subset <- data[data$Fam_idtfr == uniq_y[temp], ]
      permuted_indices <- sample(nrow(subset))
      
      ## ONLY NEED TO PERMUTE Y_SIM, AGE AND SEX SINCE THE SITE AND FAMILLY EFFECT AND LABELS ARE 
      ## SAME WHILE SHUFFLING WITHIN FAMILY
      subset$Y_sim= subset[permuted_indices, ]$Y_sim
      if(!is.na(age_eff_hat)){
      subset$age= subset[permuted_indices, ]$age
      }
      if(!is.na(sex_eff_hat)){
      subset$sex= subset[permuted_indices, ]$sex
      }
      # Use the permuted indices to reorder the rows of the data frame
      data[data$Fam_idtfr == uniq_y[temp], ] <- subset
      
    }
    if(is.na(age_eff_hat) & is.na(sex_eff_hat)){
    for (i in 1:n){
      
      brain_feature=data[,(i+1)]  ########## NEEDS TO CHANGE WHEN DIFFERENT # COVARIATES
      Y=data$Y_sim
      Ustat_mat[i,j]=as.numeric(brain_feature%*%V_inv%*%(Y-fixed_eff))
    }
    }
    # else if (is.na(age_eff_hat)){
    #   for (i in 1:n){
    #     
    #     brain_feature=data[,(i+3)]  ########## NEEDS TO CHANGE WHEN DIFFERENT # COVARIATES
    #     Y=data$Y_sim
    #     Ustat_mat[i,j]=as.numeric(brain_feature%*%V_inv%*%(Y-fixed_eff-sex*sex_eff_hat))
    #   }
    # }
    # else if (is.na(sex_eff_hat)){
    #   for (i in 1:n){
    #     
    #     brain_feature=data[,(i+3)]  ########## NEEDS TO CHANGE WHEN DIFFERENT # COVARIATES
    #     Y=data$Y_sim
    #     Ustat_mat[i,j]=as.numeric(brain_feature%*%V_inv%*%(Y-fixed_eff-age*age_eff_hat))
    #   }
    # }
    # else {
    #   for (i in 1:n){
    #     
    #     brain_feature=data[,(i+3)]  ########## NEEDS TO CHANGE WHEN DIFFERENT # COVARIATES
    #     Y=data$Y_sim
    #     Ustat_mat[i,j]=as.numeric(brain_feature%*%V_inv%*%(Y-fixed_eff-age*age_eff_hat-sex*sex_eff_hat))
    #   }
    # }
  }
  return(list(U=Ustat_mat,V=cov(Ustat_mat)))
}



#library("devtools"); install_github("lme4",user="lme4")
library(MASS)
library(lme4)


Sploc_test<-function(data,fixed_eff,sigma_err_hat,sigma_gamma_hat,sigma_tau_hat,fam_size,p){
  #sigma_tau_hat=0.5
  n=dim(data)[2]-6
  Ustat=matrix(0,n,1)
  h=dim(data)[1]
  counts <- fam_size
  
  # create the blocks
  blocks <- lapply(counts, function(x) matrix(1, nrow = x, ncol = x))
  
  # create the block diagonal matrix 
  # Very Important: dsCmatrix can be used to store sparse matrices 
  # helps a lot in compuattional time 
  # easy to invert and multiply
  diag_J <- do.call(bdiag, blocks)  # this creates a sparse matrix
  
  c_hat=(sigma_gamma_hat)^2+(sigma_tau_hat)^2 ########################
  
  V=(sigma_err_hat)^2*diag(h)+c_hat*(diag_J)
  V_inv=solve(V)
  
  Usum=0
  for (i in 1:n){
    
    brain_feature=data[,(i+1)]
    Y=data$Y_sim
    Ustat[i]=brain_feature%*%V_inv%*%(Y-fixed_eff)
    Usum=Usum+Ustat[i]
    
  }
  
  perm_test=permutation_test(data,fixed_eff,NA,NA,sigma_err_hat,sigma_gamma_hat,sigma_tau_hat,fam_size,p,V_inv,n)
  U_perm=perm_test$U # matrix of U vectors after permutation
  V_u=cov(t(U_perm))
  
  
  burden_stat=(Usum)^2/(sum(V_u))
  burden_perm=(colSums(U_perm))^2/(sum(V_u))
  burden_pval=sum(ifelse(burden_perm>burden_stat,1,0))/p 
  
  # prob=0.4
  # ver_nbd=matrix(nrow=r,ncol=r)
  # ver_nbd=apply(ver_nbd, c(1,2), function(x) sample(c(0,1),1,prob = c(1-prob,prob))) 
  # ver_nbd[upper.tri(ver_nbd)]<-0
  # ver_nbd<-ver_nbd+t(ver_nbd)
  # diag(ver_nbd)<-1
  # 
  # ## finished making vertex neighbourhood matrix
  # 
  # # Create a label matrix to keep track of the vertex indices
  # loc_label=matrix(nrow = r,ncol=r)
  # 
  # lab=1;
  # 
  # for (i in 1:r){
  #   for (j in i:c){
  #     loc_label[i,j]=lab;
  #     lab=lab+1;
  #   }
  # }
  # 
  # # symmetrizing the label matrix
  # loc_label=t(loc_label)
  # loc_label[upper.tri(loc_label)]=t(loc_label)[upper.tri(loc_label)] 
  # 
  # 
  # # Create the functional overlap matrix and metrics 
  # n_pairs=r*(r+1)/2
  # pair_nbd=matrix(nrow=n_pairs,ncol=n_pairs)
  # 
  # for (i1 in 1:r){
  #   for(j1 in 1:c){
  #     for(i2 in 1:r){
  #       for(j2 in 1:c){
  #         if (i1==i2 & j1==j2){
  #           next;
  #         } else {pair_nbd[loc_label[i1,j1],loc_label[i2,j2]]=(ver_nbd[i1,i2]+ver_nbd[i1,j2]+ver_nbd[i2,j1]+ver_nbd[j1,j2])/4;
  #         }
  #       }
  #     }
  #   }
  # }
  # 
  # func_overlp=pair_nbd #this is the functional overlap matrix
  # func_dist=1/func_overlp #functional distance is inverse of functional overlap
  # func_dist=round(func_dist,2)
  # 
  # diag(func_dist)=0 #setting diagonal =0, since distance of each pair from itself is 0
  # 
  # func_dist[!is.finite(func_dist)] <- 8
  # #View(func_dist)
  # 
  
  func_dist=read.csv("func_dist.csv")
  Gamma_rad=c(0.5,1.2,1.5,2.5,5,10)
  test_stat=matrix(nrow=sz,ncol=length(Gamma_rad))
  
  
  idx=c()
  temp=1
  dummy=r
  while(temp<=dim(func_dist)[1]){
    idx=c(idx,temp)
    temp=temp+dummy
    dummy=dummy-1
  }
  
  for (rad in 1:length(Gamma_rad)){
    for( k in 1:sz){
      w_1=ifelse( func_dist[k,]<Gamma_rad[rad],1,0) 
      w_1=w_1[-idx]
      
      
      if(w_1%*%V_u%*%w_1>0){
        test_stat[k,rad]=(w_1%*%Ustat)^2/(w_1%*%V_u%*%w_1)
      }else{
        test_stat[k,rad]=-1
      }
    }
  }
  
  sp_stat=max(test_stat)
  
  T_perm= matrix(nrow=p,1)
  
  for (num in 1:p){
    
    
    perm_test_stat=matrix(nrow=sz,ncol=length(Gamma_rad))
    
    for (rad in 1:length(Gamma_rad)){
      for( k in 1:sz){
        w_1=ifelse( func_dist[k,]<Gamma_rad[rad],1,0) 
        w_1=w_1[-idx]
        
        
        if(w_1%*%V_u%*%w_1>0){
          perm_test_stat[k,rad]=(w_1%*%U_perm[,num])^2/(w_1%*%V_u%*%w_1)
        }else{
          perm_test_stat[k,rad]=-1
        }
      }
    }
    
    T_perm[num]=max(perm_test_stat)
    
  }
  
  
  sp_pval=sum(ifelse(T_perm>sp_stat,1,0))/p
  
  #pval=length(which(abs(Usum)<abs(colSums(U_perm))))/ncol(U_perm)
  return(list(burden_pval=burden_pval,sploc_pval=sp_pval))
  
  
  
}
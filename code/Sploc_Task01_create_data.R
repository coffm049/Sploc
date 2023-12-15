##### ---- Packages

# Install packages if not installed, load. 

# p_load checks if a package is installed prior to loading into your work environment. If the package is not installed, it installs the package and then loads it.


# if (("pacman" %in% installed.packages()[,"Package"]) == FALSE) { install.packages("pacman") }
# 
# pacman::p_load(cifti, lme4, stats, tidyverse, ggpubr, R.matlab) 
library(dplyr)

I=10 #sites
J=10 # 10 families in each site
J_vec=rep(J,I) # list of number of families at the I sites (adding 1 to make it positive)


sigma_gamma=1.5
sigma_tau=.5
sigma_err=1

Gamma=rnorm(I,0,sigma_gamma) # generating the site effect gamma
gamma= rep(Gamma,J_vec) # creating repeats of the site effect for the final model
print(Gamma)
print(gamma)




tau <- unlist(lapply(seq_along(rep(1,I)), 
                     function(i) rnorm(J_vec[i], 0, sigma_tau)))
# the above gives a vector of all family effects 


#assuming family size does not depend on the site or the family 

### For different family sizes, use the following code:
###fam_size<- unlist(lapply(seq_along(tau), function(i) rpois(1,J=10)))+1
fam_size<-rep(2,length(tau))
# list of family sizes at all the sites

h=sum(fam_size) # total number of participants in the study

# dimension of fconn is 15x15
r=10
c=r

# df>r to make the simulated matrices invertible wp 1

df=r+5
Sigma=1/df*diag(r)
W= rWishart(1, df, Sigma)


fconn=c(1:(r*c*h))
dim(fconn)=c(r,c,h)

fcov=rWishart(h, df, Sigma);

for (i in 1:h){
  fconn[,,i]=cov2cor(fcov[,,i]) 
}


sz=r*(r-1)/2


########################
# If we assume the betas have a distribution, we could
# simulate from mvrnorm() as below
# But as we are doing a frequentist method
# let us assume beta to be fixed
#beta_0=mvrnorm(n = 1, mu, Sigma_b)
########################


prob_beta=0
es=1 # effect size es, ie, beta_0=(es,es,es,...(r-1 times),0,0,0..0)
#beta_0=c(rep(es,(r-1)),rep(0,sz-r+1)) #
#beta_0=c(rep(es,sz))
beta_0=es*rbinom(sz,1,prob_beta)

#es=0.1 # effect size es, ie, beta_0=(es,es,es,...(r-1 times),0,0,0..0)
#beta_0=c(rep(es,(r-1)),rep(0,sz-r+1)) #
#beta_0=c(rep(es,sz))
#beta_0=c(rep(0,sz))
#beta_0=c(es,rep(0,sz-1))


brain_ft=matrix(0,h,sz)


ind1=1

for (k in 1:h){
  ind1=1
  for (i in 1:(r-1)){
    for (j in (i+1):c){
      brain_ft[k,ind1]=fconn[i,j,k]
      ind1=ind1+1
    }
  }
}



#T=cbind(rep(rep(rep(1,length(J_vec)),J_vec),fam_size),rep(rep(rep(1,length(J_vec)),J_vec),fam_size))
T=c(1,1)
b=cbind( rep(gamma,fam_size),rep(tau,fam_size))

#Y_sim= b%*%T + rnorm(h,0,sigma_err) # Z*alpha+ error
Y_sim=brain_ft%*%beta_0  + b%*%T + rnorm(h,0,sigma_err)#  #X*beta +Z*alpha+ error # alternate true

data=data.frame(Y_sim,brain_ft,b)
#View(data)
colnames(data)[dim(data)[2]-1] <- "site_eff"
colnames(data)[dim(data)[2]] <- "fam_eff"
vec=c(1:I)
vec1=rep(vec,J_vec)
site_label=rep(vec1,fam_size)


sequences <- lapply(J_vec, seq) # create sequences
family_label <- unlist(sequences) # combine sequences
data$Site=site_label

data$Family=rep(family_label,fam_size)

sequences <- lapply(fam_size, seq) # create sequences
# individual_label <- unlist(sequences) # combine sequences
# 
# data$ind=individual_label
# Create a vector with the elements to repeat for family identifier
vec_temp <- rep(1:(sum(J_vec)))

# Define the number of times to repeat the vector
rep_time=fam_size[1] 

# Use lapply to repeat each element of the vector
result <- unlist(lapply(vec_temp, rep, times = rep_time))

#adding a column of unique family identifier over all sites 
data$Fam_idtfr=result

#View(data)
#data.frame(Y_sim,brain_ft,site_eff,fam_eff,Site,Family,Fam_idtfr)
setwd("/panfs/jay/groups/4/miran045/shared/projects/BWAS_toolbox/experiments/Burden_Sploc")
write.csv(data,"data_Sploc_SIM.csv")
getwd()

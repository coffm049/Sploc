library(tidyverse)

I <- 5 # sites
J <- 5 # 10 families in each site
J_vec <- rep(J, I) # list of number of families at the I sites (adding 1 to make it positive)
sigma_gamma <- 1.5
sigma_tau <- .5
sigma_err <- 1
Gamma <- rnorm(I, 0, sigma_gamma) # generating the site effect gamma
gamma <- rep(Gamma, J_vec) # creating repeats of the site effect for the final model
tau <- rnorm(J_vec, 0, sigma_tau) # family effects

# list of family sizes at all the sites
fam_size <- rep(2, I)
h <- sum(fam_size) # total number of participants in the study

# dimension of fconn
r <- 10

# df>r to make the simulated matrices invertible wp 1
df <- r + 5
Sigma <- 1 / df * diag(r)
W <- rWishart(1, df, Sigma)


fconn <- c(1:(r * r * h))
dim(fconn) <- c(r, r, h)

fcov <- rWishart(h, df, Sigma)
for (i in 1:h) {
  fconn[, , i] <- cov2cor(fcov[, , i])
}

sz <- r * (r - 1) / 2

prob_beta <- 0
es <- 1 # effect size es, ie, beta_0=(es,es,es,...(r-1 times),0,0,0..0)
# beta_0=c(rep(es,(r-1)),rep(0,sz-r+1)) #
# beta_0=c(rep(es,sz))
beta_0 <- es * rbinom(sz, 1, prob_beta)

brain_ft <- matrix(0, h, sz)

ind1 <- 1

for (k in 1:h) {
  ind1 <- 1
  for (i in 1:(r - 1)) {
    for (j in (i + 1):r) {
      brain_ft[k, ind1] <- fconn[i, j, k]
      ind1 <- ind1 + 1
    }
  }
}

##############################################
#### QUESTIONS ABOUT THESE LINES##############
##############################################
Tt <- c(1, 1)
# Don't know what's going on here
b <- cbind(rep(gamma, each = fam_size), rep(tau, fam_size))

Y_sim <- brain_ft %*% beta_0 + b %*% Tt + rnorm(h, 0, sigma_err) #  #X*beta +Z*alpha+ error # alternate true
#############################################
#### END QUESTION SECTION ###################
#############################################


data <- data.frame(Y_sim, brain_ft, b)
# View(data)
colnames(data)[dim(data)[2] - 1] <- "site_eff"
colnames(data)[dim(data)[2]] <- "fam_eff"
vec <- c(1:I)
vec1 <- rep(vec, J_vec)
site_label <- rep(vec1, fam_size)


sequences <- lapply(J_vec, seq) # create sequences
family_label <- unlist(sequences) # combine sequences
data$Site <- site_label

data$Family <- rep(family_label, fam_size)

sequences <- lapply(fam_size, seq) # create sequences

# Create a vector with the elements to repeat for family identifier
vec_temp <- rep(1:(sum(J_vec)))

# Define the number of times to repeat the vector
rep_time <- fam_size[1]

# Use lapply to repeat each element of the vector
result <- unlist(lapply(vec_temp, rep, times = rep_time))

# adding a column of unique family identifier over all sites
data$Fam_idtfr <- result

setwd("/panfs/jay/groups/4/miran045/shared/projects/BWAS_toolbox/experiments/Burden_Sploc")
write.csv(data, "data_Sploc_SIM.csv")
getwd()

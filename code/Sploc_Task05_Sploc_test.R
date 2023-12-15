library(MASS)
library(lme4)
ilibrary(Matrix)
source("Sploc_Task03_permutation_test.R")

Sploc_test <- function(data, fixed_eff, sigma_err_hat, sigma_gamma_hat, sigma_tau_hat, fam_size, p) {
  n <- dim(data)[2] - 6
  Ustat <- matrix(0, n, 1)
  h <- dim(data)[1]
  counts <- fam_size

  # dsCmatrix can be used to store sparse matrices
  diag_J <- do.call(bdiag, lapply(counts, function(x) {
    matrix(1, nrow = x, ncol = x)}
  ))

  c_hat <- (sigma_gamma_hat)^2 + (sigma_tau_hat)^2
  V_inv <- solve((sigma_err_hat)^2 * diag(h) + c_hat * (diag_J))

  Usum <- 0
  for (i in 1:n) {
    brain_feature <- data[, (i + 1)]
    Y <- data$Y_sim
    Ustat[i] <- brain_feature %*% V_inv %*% (Y - fixed_eff)
    Usum <- Usum + Ustat[i]
  }

  perm_test <- permutation_test(data, fixed_eff, NA, NA, sigma_err_hat, sigma_gamma_hat, sigma_tau_hat, fam_size, p, V_inv, n)
  U_perm <- perm_test$U # matrix of U vectors after permutation
  V_u <- cov(t(U_perm))

  burden_stat <- (Usum)^2 / (sum(V_u))
  burden_perm <- (colSums(U_perm))^2 / (sum(V_u))
  burden_pval <- sum(ifelse(burden_perm > burden_stat, 1, 0)) / p

  #
  # ## finished making vertex neighbourhood matrix
  #
  # # Create a label matrix to keep track of the vertex indices

  func_dist <- read.csv("func_dist.csv")
  Gamma_rad <- c(0.5, 1.2, 1.5, 2.5, 5, 10)
  test_stat <- matrix(nrow = sz, ncol = length(Gamma_rad))

  idx <- c()
  temp <- 1
  dummy <- r
  while (temp <= dim(func_dist)[1]) {
    idx <- c(idx, temp)
    temp <- temp + dummy
    dummy <- dummy - 1
  }

  for (rad in 1:length(Gamma_rad)) {
    for (k in 1:sz) {
      w_1 <- ifelse(func_dist[k, ] < Gamma_rad[rad], 1, 0)
      w_1 <- w_1[-idx]


      if (w_1 %*% V_u %*% w_1 > 0) {
        test_stat[k, rad] <- (w_1 %*% Ustat)^2 / (w_1 %*% V_u %*% w_1)
      } else {
        test_stat[k, rad] <- -1
      }
    }
  }

  sp_stat <- max(test_stat)

  T_perm <- matrix(nrow = p, 1)

  for (num in 1:p) {
    perm_test_stat <- matrix(nrow = sz, ncol = length(Gamma_rad))

    for (rad in 1:length(Gamma_rad)) {
      for (k in 1:sz) {
        w_1 <- ifelse(func_dist[k, ] < Gamma_rad[rad], 1, 0)
        w_1 <- w_1[-idx]


        if (w_1 %*% V_u %*% w_1 > 0) {
          perm_test_stat[k, rad] <- (w_1 %*% U_perm[, num])^2 / (w_1 %*% V_u %*% w_1)
        } else {
          perm_test_stat[k, rad] <- -1
        }
      }
    }

    T_perm[num] <- max(perm_test_stat)
  }


  sp_pval <- sum(ifelse(T_perm > sp_stat, 1, 0)) / p

  # pval=length(which(abs(Usum)<abs(colSums(U_perm))))/ncol(U_perm)
  return(list(burden_pval = burden_pval, sploc_pval = sp_pval))
}


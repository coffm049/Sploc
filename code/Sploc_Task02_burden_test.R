library(MASS)
library(lme4)
source("Sploc_Task03_permutation_test.R")


Burden_test <- function(data, fixed_eff, sigma_err_hat, sigma_gamma_hat, sigma_tau_hat, fam_size, p) {
  n <- dim(data)[2] - 6
  Ustat <- matrix(0, n, 1)
  h <- dim(data)[1]
  counts <- fam_size

  # create the blocks
  blocks <- lapply(counts, function(x) matrix(1, nrow = x, ncol = x))

  # create the block diagonal matrix
  # Very Important: dsCmatrix can be used to store sparse matrices
  # helps a lot in compuattional time
  # easy to invert and multiply
  diag_J <- do.call(bdiag, blocks) # this creates a sparse matrix

  c_hat <- (sigma_gamma_hat)^2 + (sigma_tau_hat)^2 ########################

  V <- (sigma_err_hat)^2 * diag(h) + c_hat * (diag_J)
  V_inv <- solve(V)

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

  pval_perm <- length(which(abs(Usum) < abs(colSums(U_perm)))) / ncol(U_perm)
  pval_exact <- pchisq(burden_stat, 1, lower.tail = FALSE)

  return(list(pval_perm = pval_perm, pval_exact = pval_exact))
}

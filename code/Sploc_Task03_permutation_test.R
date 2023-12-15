library(tidyverse)

permutation_test <- function(data, fixed_eff, age_eff_hat, sex_eff_hat, sigma_err_hat, sigma_gamma_hat, sigma_tau_hat, fam_size, p, V_inv, n) {
  Ustat_mat <- matrix(0, n, p)

  # permutation j
  for (j in 1:p) {
    # shuffle Y_sim within family
    data %>%
      group_by(Fam_idtfr) %>%
      mutate(Y_sim = sample(Y_sim))
####################################################################################  
##### WHY Does permutation get ignored if effect estimates are known?????###########
####################################################################################

    if (is.na(age_eff_hat) & is.na(sex_eff_hat)) {
      for (i in 1:n) {
        brain_feature <- data[, (i + 1)] ########## NEEDS TO CHANGE WHEN DIFFERENT # COVARIATES
        Y <- data$Y_sim
        Ustat_mat[i, j] <- as.numeric(brain_feature %*% V_inv %*% (Y - fixed_eff))
      }
    }
  }
  return(list(U = Ustat_mat, V = cov(Ustat_mat)))
}

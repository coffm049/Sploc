library(tidyverse)

permutation_test <- function(data, fixed_eff, age_eff_hat, sex_eff_hat, sigma_err_hat, sigma_gamma_hat, sigma_tau_hat, fam_size, p, V_inv, n) {
  Ustat_mat <- matrix(0, n, p)

  # permutation j
  for (j in 1:p) {
    # shuffle Y_sim within family
    Ustat_mat[,j] <- data %>%
      # permutate wihin each family
      group_by(Fam_idtfr) %>%
      mutate(Y_sim = sample(Y_sim)) %>%
      ungroup() %>%
      # summarize the Ustat (matrix mutliplcation between Y_sim and brain features [columns starting with "X_"], separately)
      summarize(across(starts_with("X_"), ~ .x %*% V_inv %*% Y_sim)) %>%
      unlist(use.names= FALSE)

    # Could also generate all permutations at once, then loop over them ussing purrr if data isn't too large
  }
  return(list(U = Ustat_mat, V = cov(Ustat_mat)))
}

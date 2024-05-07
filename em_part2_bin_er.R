em_2_er_bin_fun <- function(sub_dat, xvar_names, yvar_name,
                         p_c_dot, p_an_dot, p_d_dot,
                         epsi = 1e-6, maxit = 100, seed = 1,
                         bar = TRUE, verbose = TRUE) {

# EM algorithm for fitting the experts Q_c11 or Q_c00 via a mixture of experts
# when exclusion restriction holds.
# Here, in the case of a binary outcome Y

#### Arguments:
# sub_dat: (data.frame) a subset of data with T=0 or T=1
#         with p dimensional predictors and the column named t and
#         z for treatment and randomization respectively

# xvar_names: (character) names of the covariates

# yvar_name: (character) name of the binary outcome variable

# p_c_dot: (num) length nrow of sub_dat conditional probability of being
#           a complier given X and T=1 (Pc•1) if sub_dat is T=1
#           a complier given X and T=0 (Pc•0) if sub_dat is T=0

# p_an_dot: (num) length nrow of sub_dat conditional probability of being
#           an always taker given X and T=1 (Pa•1) if sub_dat is T=1
#           a never taker given X and T=0 (Pn•0) if sub_dat is T=0

# p_d_dot: (num) length nrow of sub_dat conditional probability of being
#           a defier given X and T=1 (Pd•1) if sub_dat is T=1
#           a defier given X and T=0 (Pd•0) if sub_dat is T=0

# epsi: (num) stopping for convergence on the gating network parameters
# maxit: (num) max no of EM iterations

# bar: if TRUE, display a progress bar
# verbose: if TRUE, print the number of iterations

#### Output:
# p_1: (num) the predictions from the fitted expert network Q_c11 or Q_c00


p <- length(xvar_names)
y <- as.numeric(sub_dat[yvar_name][, 1])

# Initialize the prior probabilities associated with the nodes of the tree
  g_1 <- p_c_dot
  g_2 <- p_an_dot
  g_3 <- p_d_dot

  set.seed(seed)
  # Initialize the parameters psi of the experts
  psi_1_param <- rnorm(p, 0, 100000)
  psi_2_param <- rnorm(p, 0, 100000)
  psi_3_param <- rnorm(p, 0, 100000)

  # Compute individual predictions from the initiated expert networks
  p_1 <- expit(as.matrix(sub_dat[xvar_names]) %*% psi_1_param)
  p_2 <- expit(as.matrix(sub_dat[xvar_names]) %*% psi_2_param)
  p_3 <- expit(as.matrix(sub_dat[xvar_names]) %*% psi_3_param)

  fitted_params <- rep(0, p)
  # Create a progress bar
    if (bar) {
        pb <- progress_bar$new(total = maxit)
    }
  for (i in 1:maxit) {
    old_params <- fitted_params
    # Compute individual contributions to each expert likelihood
    l_1 <- p_1^y * (1 - p_1)^(1 - y)
    l_2 <- p_2^y * (1 - p_2)^(1 - y)
    l_3 <- p_3^y * (1 - p_3)^(1 - y)

    # Compute posterior probabilities associated with the nodes of the tree
    sum_l <- g_1 * l_1 + g_2 * l_2 + g_3 * l_3
    h_1 <- (g_1 * l_1) / sum_l
    h_2 <- (g_2 * l_2) / sum_l
    h_3 <- (g_3 * l_3) / sum_l

    # Initialize model's formula
    f_psi <- formula(paste0(yvar_name, " ~ -1 +",
                            paste(xvar_names, collapse = " + ")))

    # For each expert network estimate parameters
    # by solving the weigthed least square problem :
    mod_q_1 <- glm(f_psi, weights = h_1, data = sub_dat, family = "binomial", method = "glm.fit") # nolint
    mod_q_2 <- glm(f_psi, weights = h_2, data = sub_dat, family = "binomial", method = "glm.fit") # nolint
    mod_q_3 <- glm(f_psi, weights = h_3, data = sub_dat, family = "binomial", method = "glm.fit") # nolint

    # Update the predictions from each expert network
    p_1 <- predict(mod_q_1, newdata = sub_dat, type = "response")
    p_2 <- predict(mod_q_2, newdata = sub_dat, type = "response")
    p_3 <- predict(mod_q_3, newdata = sub_dat, type = "response")

    # Update the parameters from the expert network mod_q_c
    fitted_params <- coef(mod_q_1)

    if (bar) {
        pb$tick()
    }

    if (sum((fitted_params - old_params)^2) / p < epsi) {
      if (verbose) {
        cat("\nEM no. 2 ER converged at iteration", i, "\n")
      }
      break
    } else if (i == maxit && verbose) {
      cat("\nEM no. 2 ER reached max no of iterations\n")
    }
  }
return(coef(mod_q_1))
}

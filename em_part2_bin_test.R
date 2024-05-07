library(progress)

expit <- function(x) 1/(1+exp(-x))

em_2_bin_fun <- function(sub_dat, xvar_names, yvar_name, pc,
                         epsi = 1e-6, maxit = 100, seed = 1,
                         bar = TRUE, verbose = TRUE) {

# EM algorithm for fitting the experts q_c1 or q_c0 via a mixture of experts
# Here, in the case of a continuous outcome Y

#### Arguments:
# sub_dat: (data.frame) a subset of data with Z=0, T=0 or Z=1, T=1
#         with p dimensional predictors and the column named t and
#         z for treatment and randomization respectively

# xvar_names: (character) names of the covariates

# yvar_name: (character) name of the binary outcome variable

# pc: (num) length nrow of sub_dat Probability of being
#           a complier vs always taker (Pc11) if sub_dat is Z=1, T=1
#           or a compiler vs never taker (Pc00) if sub_dat is Z=0, T=0

# epsi: (num) stopping for convergence on the gating network parameters
# maxit: (num) max no of EM iterations

# bar: if TRUE, display a progress bar
# verbose: if TRUE, print the number of iterations

#### Output:
# fitted_params: (num) the vector of parameters for the compiler expert
#                 q_c1 (if sub_dat is Z=1, T=1)
#                 or q_c0 (if sub_dat is Z=0, T=0)

#
p <- length(xvar_names)
y <- as.numeric(sub_dat[yvar_name][, 1])

# Initialize the prior probabilities associated with the nodes of the tree
  g_1 <- pc
  g_0 <- 1 - g_1

  set.seed(seed)
  # Initialize the parameters psi of the expert q_c(·) at random
  psi_c_param <- rnorm(p, 0, 100000)

  # Initialize the parameters psi of the expert q_nc(·) at random
  psi_nc_param <- rnorm(p, 0, 100000)

  # Compute individual predictions from the initiated expert network q_c1(·)
  p_1 <- expit(as.matrix(sub_dat[xvar_names]) %*% psi_c_param)

  # Compute individual predictions from the initiated expert network πq_c1(·)
  p_0 <- expit(as.matrix(sub_dat[xvar_names]) %*% psi_nc_param)

  fitted_params <- rep(0, 4 * p)
  # Create a progress bar
    if (bar) {
        pb <- progress_bar$new(total = maxit)
    }
  for (i in 1:maxit) {
    old_params <- fitted_params
    # Compute individual contributions to q_c1(·)’s likelihood
    l_1 <- p_1^y * (1 - p_1)^(1 - y)

    # Compute individual contributions to q_c0(·)’s likelihood
    l_0 <-  p_0^y * (1 - p_0)^(1 - y)

    # Compute posterior probabilities associated with the nodes of the tree
    h_1 <- (g_1 * l_1) / (g_0 * l_0 + g_1 * l_1)
    h_0 <- 1 - h_1

    # Initialize model's formula
    f_psi <- formula(paste0(yvar_name, " ~ -1 +",
                            paste(xvar_names, collapse = " + ")))

    # For the expert network q_c1(·) estimate parameters psi_c1
    # by solving the weigthed least square problem :
    mod_q_c <- glm(f_psi, weights = h_1, data = sub_dat, family = "binomial", method = "glm.fit") # nolint

    # For the expert network q_c0(·) estimate parameters psi_c0
    # by solving the weigthed least square problem :
    mod_q_nc <- glm(f_psi, weights = h_0, data = sub_dat, family = "binomial", method = "glm.fit") # nolint

    # Update the predictions from the expert network mod_q_c
    p_1 <- predict(mod_q_c, newdata = sub_dat, type = "response")

    # Update the predictions from the expert network mod_q_nc
    p_0 <- predict(mod_q_nc, newdata = sub_dat, type = "response")

    # Update the parameters from the expert network mod_q_c
    fitted_params <- coef(mod_q_c)

    if (bar) {
        pb$tick()
    }

    if (sum((fitted_params - old_params)^2) / p < epsi) {
      if (verbose) {
        cat("\nEM no. 2 converged at iteration", i, "\n")
      }
      break
    } else if (i == maxit && verbose) {
      cat("\nEM no. 2 reached max no of iterations\n")
    }
  }
return(list(c = coef(mod_q_c), nc = coef(mod_q_nc)))
}

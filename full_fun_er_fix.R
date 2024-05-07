er_fun <- function(eta, res, res_mon, dat, t_name, yvar_name,
            xvar_names_11, xvar_names_00,
            epsi_q, maxit_q, bar, verbose) {


#### Assuming no monotonicty, we take res as input
rho_c <- res$g_mat[, "C"]
rho_a <- res$g_mat[, "A"]
rho_n <- res$g_mat[, "N"]
rho_d <- res$g_mat[, "D"]

# Compute the propensity score (PS)
e <- rho_c * eta + rho_a + rho_d * (1 - eta)

# Bound the PS away from 0 and 1 for numerical stability
e <- pmin(1 - 1e-6, pmax(1e-6, e))

# Compute the conditional probabilities given X and T=1
p_c_dot1 <- (eta / e) * rho_c
p_a_dot1 <- (rho_a / e)
p_d_dot1 <- ((1 - eta) / e) * rho_d

# Normalize the conditional probabilities so they always sum to 1
p_dot1_sum <- p_c_dot1 + p_a_dot1 + p_d_dot1
p_c_dot1 <- p_c_dot1 / p_dot1_sum
p_a_dot1 <- p_a_dot1 / p_dot1_sum
p_d_dot1 <- p_d_dot1 / p_dot1_sum

# Compute the conditional probabilities given X and T=0
p_c_dot0 <- ((1 - eta) / (1 - e)) * rho_c
p_n_dot0 <- (rho_n / (1 - e))
p_d_dot0 <- (eta / (1 - e)) * rho_d

# Normalize the conditional probabilities so they always sum to 1
p_dot0_sum <- p_c_dot0 + p_n_dot0 + p_d_dot0
p_c_dot0 <- p_c_dot0 / p_dot0_sum
p_n_dot0 <- p_n_dot0 / p_dot0_sum
p_d_dot0 <- p_d_dot0 / p_dot0_sum

# EM algorithm for fitting the experts Q_c11
q_11_params <- em_2_er_bin_fun(sub_dat = dat[dat[, t_name] == 1, ],
                xvar_names = xvar_names_11, yvar_name = yvar_name,
                p_c_dot = p_c_dot1[dat[, t_name] == 1],
                p_an_dot = p_a_dot1[dat[, t_name] == 1],
                p_d_dot = p_d_dot1[dat[, t_name] == 1],
                epsi = epsi_q, maxit = maxit_q,
                bar = bar, verbose = verbose)

# EM algorithm for fitting the experts Q_c00
q_00_params <- em_2_er_bin_fun(sub_dat = dat[dat[, t_name] == 0, ],
                xvar_names = xvar_names_00, yvar_name = yvar_name,
                p_c_dot = p_c_dot0[dat[, t_name] == 0],
                p_an_dot = p_n_dot0[dat[, t_name] == 0],
                p_d_dot = p_d_dot0[dat[, t_name] == 0],
                epsi = epsi_q, maxit = maxit_q,
                bar = bar, verbose = verbose)

# Compute the predictions from the fitted Q_c11
q_11 <- expit(as.matrix(dat[xvar_names_11]) %*% q_11_params)

# Compute the predictions from the fitted Q_c00
q_00 <- expit(as.matrix(dat[xvar_names_00]) %*% q_00_params)

# Compute the CACE
p_c <- as.numeric(rho_c)
no_mon_er_res <- sum((q_11 - q_00) * p_c) / sum(p_c)


### Assuming monotonicty, we take res_mon as input
rho_c <- res_mon$g_mat[, "C"]
rho_a <- res_mon$g_mat[, "A"]
rho_n <- res_mon$g_mat[, "N"]
rho_d <- 0

# Compute the propensity score (PS)
e <- rho_c * eta + rho_a

# Bound the PS away from 0 and 1 for numerical stability
e <- pmin(1 - 1e-6, pmax(1e-6, e))

# Compute the conditional probabilities given X and T=1
p_c_dot1 <- (eta / e) * rho_c
p_a_dot1 <- (rho_a / e)

# Normalize the conditional probabilities so they always sum to 1
p_dot1_sum <- p_c_dot1 + p_a_dot1
p_c_dot1 <- p_c_dot1 / p_dot1_sum

# Compute the conditional probabilities given X and T=0
p_c_dot0 <- ((1 - eta) / (1 - e)) * rho_c
p_n_dot0 <- (rho_n / (1 - e))

# Normalize the conditional probabilities so they always sum to 1
p_dot0_sum <- p_c_dot0 + p_n_dot0
p_c_dot0 <- p_c_dot0 / p_dot0_sum

# EM algorithm for fitting the experts Q_c11
q_11_params <- em_2_bin_fun(sub_dat = dat[dat[, t_name] == 1, ],
                xvar_names = xvar_names_11, yvar_name = yvar_name,
                pc = p_c_dot1[dat[, t_name] == 1],
                epsi = epsi_q, maxit = maxit_q,
                bar = bar, verbose = verbose)

q_00_params <- em_2_bin_fun(sub_dat = dat[dat[, t_name] == 0, ],
                xvar_names = xvar_names_00, yvar_name = yvar_name,
                pc = p_c_dot0[dat[, t_name] == 0],
                epsi = epsi_q, maxit = maxit_q,
                bar = bar, verbose = verbose)

# Compute the predictions from the fitted Q_c11
q_11 <- expit(as.matrix(dat[xvar_names_11]) %*% q_11_params$c)

# Compute the predictions from the fitted Q_c00
q_00 <- expit(as.matrix(dat[xvar_names_00]) %*% q_00_params$c)

# Compute the CACE
p_c <- as.numeric(rho_c)
mon_er_res <- sum((q_11 - q_00) * p_c) / sum(p_c)

return(list(no_mon_er_res = no_mon_er_res, mon_er_res = mon_er_res))
}

no_er_fun <- function(res, dat, z_name, t_name, yvar_name,
            xvar_names_11, xvar_names_00,
            epsi_q, maxit_q, bar, verbose) {
    p_c11 <- with(as.data.frame(res$g_mat), C / (C + A))[dat[, z_name] == 1 & dat[, t_name] == 1] # nolint
    p_c00 <- with(as.data.frame(res$g_mat), C / (C + N))[dat[, z_name] == 0 & dat[, t_name] == 0] # nolint

    sub_dat_11 <- dat[dat[, z_name] == 1 & dat[, t_name] == 1, ]
    sub_dat_00 <- dat[dat[, z_name] == 0 & dat[, t_name] == 0, ]

    beta <- em_2_bin_fun(sub_dat = sub_dat_11, xvar_names = xvar_names_11,
                yvar_name = yvar_name, pc = p_c11,
                epsi = epsi_q, maxit = maxit_q,
                bar = bar, verbose = verbose)

    gamma <- em_2_bin_fun(sub_dat = sub_dat_00, xvar_names = xvar_names_00,
                yvar_name = yvar_name, pc = p_c00,
                epsi = epsi_q, maxit = maxit_q,
                bar = bar, verbose = verbose)

    q_11 <- expit(as.matrix(dat[xvar_names_11]) %*% beta$c)
    q_00 <- expit(as.matrix(dat[xvar_names_00]) %*% gamma$c)
    p_c <- as.numeric(res$g_mat[, "C"])

    no_er_res <- sum((q_11 - q_00) * p_c) / sum(p_c)
    return(list(no_er_res = no_er_res, beta = beta, gamma = gamma))
}

cace_bin_fun <- function(dat,
        xvar_names = NULL, z_name = "z", t_name = "ttt", yvar_name = "y",
        xvar_names_cand = NULL,
        xvar_names_11 = NULL, xvar_names_00 = NULL,
        xvar_names_y = NULL,
        xvar_names_z = NULL,
        maxit_cand = 100, epsi_cand = 1e-6,
        maxit_q = 500, epsi_q = 1e-6,
        bar = FALSE, verbose = TRUE) {

if (!is.null(xvar_names)) {
xvar_names_cand <- xvar_names
xvar_names_11 <- xvar_names
xvar_names_00 <- xvar_names
xvar_names_y <- xvar_names
xvar_names_z <- z_name
}
res_mon <- em_1_mon_fun(dat, xvar_names = xvar_names_cand,
            z_name = z_name, t_name = t_name,
            maxit = maxit_cand, epsi = epsi_cand,
            bar = bar, verbose = verbose)

# Initialize models eta formulas
f_eta <- formula(paste0(z_name, " ~ -1 +",
                        paste(xvar_names_z, collapse = " + ")))

# Fit the eta model
mod_eta <- glm(f_eta, data = dat, family = "binomial", method = "glm.fit")

# Predictions from the models
eta <- predict(mod_eta, newdata = dat, type = "response")

#### No monotonicity assumption
res <- em_1_fun(dat, xvar_names = xvar_names_cand,
            z_name = z_name, t_name = t_name,
            maxit = maxit_cand, epsi = epsi_cand,
            bar = bar, verbose = verbose)

### No exclusion restriction
# without monotonicity assumption
no_mon_no_er_res <- no_er_fun(res, dat,
            z_name, t_name, yvar_name,
            xvar_names_11, xvar_names_00,
            epsi_q, maxit_q, bar, verbose)

# with monotonicity assumption
mon_no_er_res <- no_er_fun(res_mon, dat,
            z_name, t_name, yvar_name,
            xvar_names_11, xvar_names_00,
            epsi_q, maxit_q, bar, verbose)

### Exclusion restriction
# without monotonicity assumption
er_res <- er_fun(eta, res, res_mon, dat,
            t_name, yvar_name,
            xvar_names_11, xvar_names_00,
            epsi_q, maxit_q, bar, verbose)


### Basic Wald estimator from Angrist et al. JASA, 1996?
itt_hat <- mean(dat[dat[, z_name] == 1, yvar_name]) - mean(dat[dat[, z_name] == 0, yvar_name]) # nolint
rho_c_hat <- mean(dat[dat[, z_name] == 1, t_name]) - mean(dat[dat[, z_name] == 0, t_name]) # nolint
wald_res <- itt_hat / rho_c_hat

### Propensity score matching estimator
### Equation (14) in Frolich, Journal of Econometrics, 2006
# Initialize models q, eta and e formulas
f_m <- formula(paste0(yvar_name, " ~ -1 + eta"))
f_d <- formula(paste0(t_name, " ~ -1 + eta"))

# Fit the models
dat$eta <- eta
mod_m_pi1 <- glm(f_m, data = dat[dat[, z_name] == 1, ], family = "binomial", method = "glm.fit") # nolint
mod_m_pi0 <- glm(f_m, data = dat[dat[, z_name] == 0, ], family = "binomial", method = "glm.fit") # nolint

mod_mu_pi1 <- glm(f_d, data = dat[dat[, z_name] == 1, ], family = "binomial", method = "glm.fit") # nolint
mod_mu_pi0 <- glm(f_d, data = dat[dat[, z_name] == 0, ], family = "binomial", method = "glm.fit") # nolint

# Predictions from the models
m_pi1 <- predict(mod_m_pi1, newdata = dat, type = "response")
m_pi0 <- predict(mod_m_pi0, newdata = dat, type = "response")

mu_pi1 <- predict(mod_mu_pi1, newdata = dat, type = "response")
mu_pi0 <- predict(mod_mu_pi0, newdata = dat, type = "response")

# Compute the PS matching estimator
ps_matching <- sum(m_pi1 - m_pi0) / sum(mu_pi1 - mu_pi0)

# Put estimates in the right order for the output
cace_estimates <- cbind(
             "no_mon_no_er_res" = no_mon_no_er_res$no_er_res,
             "mon_no_er_res" = mon_no_er_res$no_er_res,
             "no_mon_er_res" = er_res$no_mon_er_res,
             "mon_er_res" = er_res$mon_er_res,
             "wald_res" = wald_res,
             "ps_matching" = ps_matching)

return(cace_estimates)
}

#############
#############

# source("gen_dat_full.R")

# p <- 15
# xvar_names <- paste0("X", 0:p)

# sample sizes to explore c(800, 4000, 20000)
#dat <- gen_data(n = 200, p = p, alpha = .5,
#                bin = TRUE, seed = 1, def = FALSE, er = TRUE)

#temp <- cace_bin_fun(dat = dat, xvar_names = xvar_names,
#            z_name = "z", t_name = "ttt", yvar_name = "y",
#            bar = TRUE, verbose = TRUE,
#            epsi_cand = 1e-3, epsi_q = 1e-3, maxit_cand = 100, maxit_q = 100)

# temp
# with(dat[dat$cand==1, ], mean(y_t1 - y_t0))
# 100 * abs(temp - with(dat[dat$cand==1,], mean(y_t1 - y_t0)))

emp_boot_ci <-  function(t0, t, alpha=.05){
  # Empirical Bootstrap CI
  # as described in John Rice, Mathematical Statistics and Data Analysis,
  # 3rd edition, p. 285.
  temp <- rbind(
    2 * t0 - apply(t, 2, function(x) quantile(x, probs = 1 - alpha / 2)),
    2 * t0 - apply(t, 2, function(x) quantile(x, probs = alpha / 2))
  )
  return(temp)
}

library(boot)

cace_bin_fun_boot <- function(data, indices, ...) {
  counter <<- counter + 1
  frac <- 100 * (counter - 1) / r_boot
  if (counter > 0) {
        cat("Bootstrapping at", round(frac), "%\n")
  } else {
        cat("Original estimation\n")
  }
  bootstrap_sample <- data[indices, ]
  return(cace_bin_fun(dat = bootstrap_sample, ...))
}

cace_bin_boot <- function(dat, r_boot,
                          par, ncpus = parallel::detectCores(logical=FALSE),
                          xvar_names, z_name, t_name, yvar_name,
                          bar, verbose,
                          epsi_cand, epsi_q, maxit_cand, maxit_q) {
r_boot <<- r_boot
counter <<- -1

boot_res <- boot(data = dat, statistic = cace_bin_fun_boot, R = r_boot,
                 parallel = par, ncpus = ncpus,
                 xvar_names = xvar_names, z_name = z_name, t_name = t_name, yvar_name = yvar_name, # nolint
                 bar = bar, verbose = verbose,
                 epsi_cand = epsi_cand, epsi_q = epsi_q, maxit_cand = maxit_cand, maxit_q = maxit_q) # nolint

ci_s <- emp_boot_ci(boot_res$t0, boot_res$t, alpha = .05)

res <- c(c(boot_res$t0), c(ci_s))

names(res) <- c(colnames(boot_res$t0),
paste0(rep(colnames(ci_s), each = 2), rep(c("_lb", "_ub"), times = 6)))

return(res)
}

# dat <- gen_data(n = 10000, p = 6, alpha = .5,
#                bin = TRUE, seed = 1960, def = FALSE, er = TRUE)

# mean(dat$cand == 1)

# start_time <- Sys.time()
# res <- cace_bin_boot(dat = dat, r_boot = 0,
#              xvar_names = paste0("X", 0:6), z_name = "z", t_name = "ttt", yvar_name = "y",
#              par = "no", ncpus = parallel::detectCores(logical=FALSE),
#              bar = TRUE, verbose = TRUE,
#              epsi_cand = 1e-6, epsi_q = 1e-6, maxit_cand = 200, maxit_q = 200)
# end_time <- Sys.time()
# end_time - start_time
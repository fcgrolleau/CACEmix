setwd("/Users/francois/Desktop/github repos/CACEmix/")

source("em_part1.R")
source("em_part1_mon.R")
source("em_part2_bin_test.R")
source("em_part2_bin_er.R")
source("full_fun_er_fix2.R")
source("simulation_funs.R")
source("gen_dat_full.R")

# Save the results in a new directory titled with today's date
sim_dir <- paste0("sim_res/", Sys.time())
dir.create(sim_dir, recursive = TRUE)
setwd(sim_dir)

### define the data generating processes

seeds <- c(238, 238, 238, 238)

dgp_ori_fun <- function(n) {
        gen_data(n, p = 14, icpt = TRUE,
                sd_y = 0, np_noise = 0, eigeinval = NULL,
                alpha = .5, bin = TRUE, seed = seeds[1], verbose = FALSE,
                def = TRUE, er = FALSE)
}

dgp_er_fun <- function(n) {
        gen_data(n, p = 14, icpt = TRUE,
                sd_y = 0, np_noise = 0, eigeinval = NULL,
                alpha = .5, bin = TRUE, seed = seeds[2], verbose = FALSE,
                def = TRUE, er = TRUE)
}

dgp_mon_fun <- function(n) {
        gen_data(n, p = 14, icpt = TRUE,
                sd_y = 0, np_noise = 0, eigeinval = NULL,
                alpha = .5, bin = TRUE, seed = seeds[3], verbose = FALSE,
                def = FALSE, er = FALSE)
}

dgp_mon_er_fun <- function(n) {
        gen_data(n, p = 14, icpt = TRUE,
                sd_y = 0, np_noise = 0, eigeinval = NULL,
                alpha = .5, bin = TRUE, seed = seeds[4], verbose = FALSE,
                def = FALSE, er = TRUE)
}

### define the estimation procedures

it_fun <- function(dat) {
R.utils::withTimeout({
        temp <- cace_bin_boot(dat = dat, r_boot = 0,
                      xvar_names = paste0("X", 0:14),
                      z_name = "z", t_name = "ttt", yvar_name = "y",
                      par = "no", bar = FALSE, verbose = TRUE,
                      epsi_cand = 1e-5, epsi_q = 1e-5,
                      maxit_cand = 200, maxit_q = 200)
}, timeout = 150, onTimeout = "error")
        return(temp)
}

it_mis_fun <- function(dat) {
R.utils::withTimeout({
        temp <- cace_bin_boot(dat = dat, r_boot = 0,
                      xvar_names = c(paste0("X", 0:6), paste0("X", 8:13)),
                      z_name = "z", t_name = "ttt", yvar_name = "y",
                      par = "no", bar = FALSE, verbose = TRUE,
                      epsi_cand = 1e-5, epsi_q = 1e-5,
                      maxit_cand = 200, maxit_q = 200)
}, timeout = 150, onTimeout = "error")
        return(temp)
}

### Run the simulations
n_sim <- 300     # number of simulations
n_samples <- c(2000, 5000, 10000)  # sample sizes

# Correctly specified models

res_ori_n1 <- sim_fun(dgp_ori_fun, it_fun, n_sim = n_sim, n_samples = n_samples[1]); save.image(file="res_ori_n1.RData") # nolint
res_er_n1 <-  sim_fun(dgp_er_fun, it_fun, n_sim = n_sim, n_samples = n_samples[1]); save.image(file="res_er_n1.RData") # nolint
res_mon_n1 <- sim_fun(dgp_mon_fun, it_fun, n_sim = n_sim, n_samples = n_samples[1]) ; save.image(file="res_mon_n1.RData") # nolint
res_mon_er_n1 <- sim_fun(dgp_mon_er_fun, it_fun, n_sim = n_sim, n_samples = n_samples[1]); save.image(file="res_mon_er_n1.RData") # nolint

res_ori_n2 <- sim_fun(dgp_ori_fun, it_fun, n_sim = n_sim, n_samples = n_samples[2]); save.image(file="res_ori_n2.RData") # nolint
res_er_n2 <-  sim_fun(dgp_er_fun, it_fun, n_sim = n_sim, n_samples = n_samples[2]); save.image(file="res_er_n2.RData") # nolint
res_mon_n2 <- sim_fun(dgp_mon_fun, it_fun, n_sim = n_sim, n_samples = n_samples[2]); save.image(file="res_mon_n2.RData") # nolint
res_mon_er_n2 <- sim_fun(dgp_mon_er_fun, it_fun, n_sim = n_sim, n_samples = n_samples[2]); save.image(file="res_mon_er_n2.RData") # nolint

res_ori_n3 <- sim_fun(dgp_ori_fun, it_fun, n_sim = n_sim, n_samples = n_samples[3]); save.image(file="res_ori_n3.RData") # nolint
res_er_n3 <-  sim_fun(dgp_er_fun, it_fun, n_sim = n_sim, n_samples = n_samples[3]); save.image(file="res_er_n3.RData") # nolint
res_mon_n3 <- sim_fun(dgp_mon_fun, it_fun, n_sim = n_sim, n_samples = n_samples[3]); save.image(file="res_mon_n3.RData") # nolint
res_mon_er_n3 <- sim_fun(dgp_mon_er_fun, it_fun, n_sim = n_sim, n_samples = n_samples[3]); save.image(file="res_mon_er_n3.RData") # nolint

# Mispecified models

res_ori_n1_mis <- sim_fun(dgp_ori_fun, it_mis_fun, n_sim = n_sim, n_samples = n_samples[1]); save.image(file="res_ori_n1_mis.RData") # nolint
res_er_n1_mis <-  sim_fun(dgp_er_fun, it_mis_fun, n_sim = n_sim, n_samples = n_samples[1]); save.image(file="res_er_n1_mis.RData") # nolint
res_mon_n1_mis <- sim_fun(dgp_mon_fun, it_mis_fun, n_sim = n_sim, n_samples = n_samples[1]); save.image(file="res_mon_n1_mis.RData") # nolint
res_mon_er_n1_mis <- sim_fun(dgp_mon_er_fun, it_mis_fun, n_sim = n_sim, n_samples = n_samples[1]); save.image(file="res_mon_er_n1_mis.RData") # nolint

res_ori_n2_mis <- sim_fun(dgp_ori_fun, it_mis_fun, n_sim = n_sim, n_samples = n_samples[2]); save.image(file="res_ori_n2_mis.RData") # nolint
res_er_n2_mis <-  sim_fun(dgp_er_fun, it_mis_fun, n_sim = n_sim, n_samples = n_samples[2]); save.image(file="res_er_n2_mis.RData") # nolint
res_mon_n2_mis <- sim_fun(dgp_mon_fun, it_mis_fun, n_sim = n_sim, n_samples = n_samples[2]); save.image(file="res_mon_n2_mis.RData") # nolint
res_mon_er_n2_mis <- sim_fun(dgp_mon_er_fun, it_mis_fun, n_sim = n_sim, n_samples = n_samples[2]); save.image(file="res_mon_er_n2_mis.RData") # nolint

res_ori_n3_mis <- sim_fun(dgp_ori_fun, it_mis_fun, n_sim = n_sim, n_samples = n_samples[3]); save.image(file="res_ori_n3_mis.RData") # nolint
res_er_n3_mis <-  sim_fun(dgp_er_fun, it_mis_fun, n_sim = n_sim, n_samples = n_samples[3]); save.image(file="res_er_n3_mis.RData") # nolint
res_mon_n3_mis <- sim_fun(dgp_mon_fun, it_mis_fun, n_sim = n_sim, n_samples = n_samples[3]); save.image(file="res_mon_n3_mis.RData") # nolint
res_mon_er_n3_mis <- sim_fun(dgp_mon_er_fun, it_mis_fun, n_sim = n_sim, n_samples = n_samples[3]); save.image(file="res_mon_er_n3_mis.RData") # nolint

### Save the results

save(res_ori_n1, res_er_n1, res_mon_n1, res_mon_er_n1,
      res_ori_n2, res_er_n2, res_mon_n2, res_mon_er_n2,
      res_ori_n3, res_er_n3, res_mon_n3, res_mon_er_n3,
      res_ori_n1_mis, res_er_n1_mis, res_mon_n1_mis, res_mon_er_n1_mis,
      res_ori_n2_mis, res_er_n2_mis, res_mon_n2_mis, res_mon_er_n2_mis,
      res_ori_n3_mis, res_er_n3_mis, res_mon_n3_mis, res_mon_er_n3_mis,
 file = "simulations_results.RData")

### Check the results

length(unique(res_ori_n1[, "it"]))
100 * var(res_ori_n1[, "true_sample_cace"])
res_ori_n1[rownames(res_ori_n1)  %in% c("wald_res", "ps_matching"), c("it", "est")] # nolint
res_ori_n1[rownames(res_ori_n1) == "no_mon_no_er_res", c("it", "est", "true_sample_cace", "true_cace")] # nolint
100 * (mean(res_ori_n1[rownames(res_ori_n1) == "no_mon_no_er_res", "true_sample_cace"]) - res_ori_n1[1, "true_cace"]) # nolint best possible result
100 * (mean(res_ori_n1[rownames(res_ori_n1) == "no_mon_no_er_res", "est"]) - res_ori_n1[1, "true_cace"]) # nolint result

for (i in 43) {
res <- gen_data(1000, p = 14, icpt = TRUE,
                sd_y = 0, np_noise = 0, eigeinval = NULL,
                alpha = .5, bin = TRUE, seed = i, verbose = FALSE,
                def = TRUE, er = FALSE)

prob_c <- c(mean(res$rho_1), sd(res$rho_1) / sqrt(length(res)))
prob_d <- c(mean(res$rho_4), sd(res$rho_4) / sqrt(length(res)))
er_1 <- mean(with(res, y_a11 != y_a01))
er_0 <- mean(with(res, y_n10 != y_n00))
er_cond <- er_1 + er_0

condition <- prob_c[1] > .4 & prob_c[1]  < .6 & prob_d[1] > .1 & prob_d[1] < .2 & er_cond > 1.1 # nolint

print(
        paste0("i=", i, " ", condition,
        " ", round(prob_c[1], 2), " ", round(prob_d[1], 2), " ", round(er_cond, 2)) # nolint
        )
if (condition) break
}

#temp_dat <- gen_data(2000, p = 14, icpt = TRUE,
#                sd_y = 0, np_noise = 0, eigeinval = NULL,
#                alpha = .5, bin = TRUE, seed = 1960, verbose = FALSE,
#                def = TRUE, er = TRUE)

#temp_res <- cace_bin_boot(dat = temp_dat, r_boot = 0,
#                      xvar_names = paste0("X", 0:14),
#                      z_name = "z", t_name = "ttt", yvar_name = "y",
#                      par = "no", bar = TRUE, verbose = TRUE,
#                      epsi_cand = 1e-5, epsi_q = 1e-5,
#                      maxit_cand = 200, maxit_q = 200)

#with(temp_dat[temp_dat$cand==1,], mean(y_t1 - y_t0))
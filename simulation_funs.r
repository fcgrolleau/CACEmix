library(foreach)
library(doParallel)

# define a function to generate the index of the samples for each simulation
sim_n_fun <- function(n_sim = 100, n_samples = 1000) {
    return(
    sapply(0 : (n_sim - 1),
    function(i) (i * n_samples + 1) : (i * n_samples + n_samples))
    )
}

sim_fun <- function(gen_data_fun, it_fun, n_sim, n_samples) {
# Prepare the simulation
sim_mat <- sim_n_fun(n_sim = n_sim, n_samples = n_samples)
dat <- gen_data_fun(n = n_sim * n_samples)
true_cace <- with(dat[dat$cand == 1, ], mean(y_t1 - y_t0))
com_prop <- mean(dat$cand == 1)
def_prop <- mean(dat$cand == 4)
er1 <- with(dat, mean((y_a01 - y_a11)^2))
er0 <- with(dat, mean((y_n00 - y_n10)^2))
true_ate <- with(dat, mean(y_t1 - y_t0))

# Prepare to parallelize
registerDoParallel(cores = parallel::detectCores(logical = FALSE) - 1)

# Run the simulation
sim_res <- foreach(i = 1:n_sim, .combine = rbind,
                .errorhandling = "remove") %dopar% {
    if (i == 1 || (i %% 10) == 0) {
        cat("sample no.", i, "/", n_sim, "begins\n")}
    dat_i <- dat[sim_mat[, i], ]
    true_sample_cace <- with(dat_i[dat_i$cand == 1, ], mean(y_t1 - y_t0))
    com_prop_sample <- mean(dat_i$cand == 1)
    def_prop_sample <- mean(dat_i$cand == 4)
    er1_sample <- with(dat_i, mean((y_a01 - y_a11)^2))
    er0_sample <- with(dat_i, mean((y_n00 - y_n10)^2))
    true_sample_ate <- with(dat_i, mean(y_t1 - y_t0))

    cbind("it" = i,
          "est" = it_fun(dat_i),
          "true_sample_cace" = true_sample_cace,
          "com_prop_sample" = com_prop_sample,
          "def_prop_sample" = def_prop_sample,
          "er1_sample" = er1_sample,
          "er0_sample" = er0_sample,
          "true_sample_ate" = true_sample_ate,
          "true_cace" = true_cace,
          "com_prop" = com_prop,
          "def_prop" = def_prop,
          "er1" = er1,
          "er0" = er0,
          "true_ate" = true_ate,
          "n_samples" = n_samples)
}
stopImplicitCluster()
return(sim_res)
}

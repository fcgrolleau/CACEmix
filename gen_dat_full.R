expit <- function(x) 1 / (1 + exp(-x))

gen_data <- function(n, p, p_bin = as.integer(p / 2), icpt = TRUE,
                     eigeinval = NULL, alpha = .5, sd_y = 1, np_noise = 0,
                     bin = TRUE, seed = 42, verbose = TRUE,
                     def = TRUE, er = FALSE) {
## Generate data from a CACE model with 4 classes of compliers, always takers

## Arguments
# n: number of observations
# p: number of covariates
# p_bin: number of binary covariates
# icpt: if TRUE, an intercept is added to the covariates
# alpha: randomization probability
# sd_y: standard deviation of the potential outcomes if they are continuous
# np_noise: non parametric noise
#   if ≠ 0, the parameters do not fully determine conditional expectations
# eigeinval: eigenvalues of the covariance matrix
# bin: if TRUE, the potential outcomes are binary, otherwise they are continuous
# seed: seed for reproducibility
# verbose: if TRUE, print the parameters of the model
# def: if TRUE, a defier class is added (i.e., monotonicity assumption does not hold) # nolint
# er: if TRUE, the exclusion restriction assumption holds

## Output
# sim_dat: simulated data frame

## Example
# sim_dat <- gen_data(n = 1000, p = 2, alpha = .5, bin = TRUE, seed = 42, def = TRUE, er = FALSE) # nolint

    set.seed(seed)

    # generate a random orthonormal matrix
    library(pracma)
    o_mat <- randortho(p, type = "orthonormal")

    # fill diagonal matrix with eigeinvalues
    d_mat <- matrix(0, nrow = p, ncol = p)
    eigeinval <- ifelse(is.null(eigeinval), seq(1, by = .2, length = p),
                                            eigeinval)
    diag(d_mat) <- eigeinval

    # generate a positive semi-definite covariance matrix
    sigma_mat <- o_mat %*% d_mat %*% t(o_mat)

    # sample from a multivariate normal distribution N(0, sigma_mat)
    library(mvtnorm)
    x_mat <- rmvnorm(n, mean = rep(0, nrow(sigma_mat)), sigma = sigma_mat)

    set.seed(seed + 1)
    # make covariates 1:p_bin binary
    x_mat[, 1:p_bin] <- ifelse(x_mat[, 1:p_bin] > 0, 1, 0)

    # make the remaining covariates positive (i.e., log-noramlly distributed)
    x_mat[, (p_bin + 1):p] <- exp(x_mat[, (p_bin + 1):p])

    if (icpt) {
    x_mat <- cbind(1, x_mat) # add intercept
    colnames(x_mat) <- paste0("X", 0:p)
    p <- p + 1
    } else {
    colnames(x_mat) <- paste0("X", 1:p)
    }

    set.seed(seed + 2)
    z <- rbinom(n, 1, alpha) # allocation variable

    # Parameters for the model of the probability of being in each class given x
    set.seed(seed + 3)
    delta_1 <- runif(p, -1, 1)
    delta_2 <- runif(p, -1, 1)
    delta_3 <- runif(p, -1, 1)
    if (def) {
        delta_4 <- runif(p, -1, 1)
    }

    # put parameters in a (4 x p) matrix

    if (def) {
        delta_mat <- as.matrix(rbind(delta_1, delta_2, delta_3, delta_4))
    } else {
        delta_mat <- as.matrix(rbind(delta_1, delta_2, delta_3))
    }

    colnames(delta_mat) <- colnames(x_mat)

    # parametric model of the probability of being in each class given x
    rho <- function(l, x, delta_mat) {
    # l: class of interest 1 complier, 2 always taker, 3 never taker, 4 defier
    # x: vector of length p
    # delta_mat: model's parameter in a (4 x p) matrix
    exp(delta_mat[l, ] %*% x) / sum(
        sapply(seq_len(nrow(delta_mat)), function(l)  exp(delta_mat[l, ] %*% x))
        ) # (this is a softmax function)
    }

    # for all x_i from the rho model predict the probability of c,a,n,d
    rhos_mat <- t(sapply(as.data.frame(t(x_mat)),
        function(x_i) rho(seq_len(nrow(delta_mat)), x_i, delta_mat)))
    colnames(rhos_mat) <- paste0("rho_", seq_len(nrow(delta_mat)))

    # sample an unique class from the probability of being in each class
    class_mat <- t(sapply(as.data.frame(t(rhos_mat)),
        function(x_i) rmultinom(n = 1, size = 1, prob = x_i)))
    cand <- apply(class_mat, 1, which.max)

    # from that unique class and the value of z (allocation at randomization)
    # predict the treatment actually taken
    ttt <- (cand == 1) * z + (cand == 2) + (cand == 4) * (1 - z)

    # Generate potential outcomes
    # Parameters for the potential outcome models
    set.seed(seed + 4)
    beta_c11 <- runif(p, -1, 1)
    beta_a11 <- runif(p, -1, 1)
    beta_a01 <- runif(p, -1, 1)
    beta_d01 <- runif(p, -1, 1)

    beta_c00 <- runif(p, -1, 1)
    beta_n10 <- runif(p, -1, 1)
    beta_n00 <- runif(p, -1, 1)
    beta_d10 <- runif(p, -1, 1)

    beta_mat <- as.matrix(rbind(beta_c11, beta_a11, beta_a01, beta_d01,
                                beta_c00, beta_n10, beta_n00, beta_d10))
    colnames(beta_mat) <- colnames(x_mat)

    # Compute dot products ± add non-parametric noise
    set.seed(seed + 5)
    dp_c11 <- x_mat %*% beta_c11 +  rexp(n, rate = 1 / np_noise)
    dp_a11 <- x_mat %*% beta_a11 +  rexp(n, rate = 1 / np_noise)
    dp_a01 <- x_mat %*% beta_a01 +  rexp(n, rate = 1 / np_noise)
    dp_d01 <- x_mat %*% beta_d01 +  rexp(n, rate = 1 / np_noise)

    dp_c00 <- x_mat %*% beta_c00 +  rexp(n, rate = 1 / np_noise)
    dp_n10 <- x_mat %*% beta_n10 +  rexp(n, rate = 1 / np_noise)
    dp_n00 <- x_mat %*% beta_n00 +  rexp(n, rate = 1 / np_noise)
    dp_d10 <- x_mat %*% beta_d10 +  rexp(n, rate = 1 / np_noise)

    # Generate elementary potential outcomes
    if (bin) {
    set.seed(seed + 6)
    y_c11 <- rbinom(n, 1, expit(dp_c11))
    y_a11 <- rbinom(n, 1, expit(dp_a11))
    y_a01 <- rbinom(n, 1, expit(dp_a01))
    y_d01 <- rbinom(n, 1, expit(dp_d01))

    y_c00 <- rbinom(n, 1, expit(dp_c00))
    y_n10 <- rbinom(n, 1, expit(dp_n10))
    y_n00 <- rbinom(n, 1, expit(dp_n00))
    y_d10 <- rbinom(n, 1, expit(dp_d10))
    } else {
    y_c11 <- rnorm(n, mean = dp_c11, sd = sd_y)
    y_a11 <- rnorm(n, mean = dp_a11, sd = sd_y)
    y_a01 <- rnorm(n, mean = dp_a01, sd = sd_y)
    y_d01 <- rnorm(n, mean = dp_d01, sd = sd_y)

    y_c00 <- rnorm(n, mean = dp_c00, sd = sd_y)
    y_n10 <- rnorm(n, mean = dp_n10, sd = sd_y)
    y_n00 <- rnorm(n, mean = dp_n00, sd = sd_y)
    y_d10 <- rnorm(n, mean = dp_d10, sd = sd_y)
    }

    # If exclusion restriction holds, we have:
    if (er) {
    y_a11 <- y_a01
    y_n00 <- y_n10
    }

    y_t1 <- (cand == 1) * y_c11 + (cand == 2) * z * y_a11 + (cand == 2) * (1 - z) * y_a01 + (cand == 4) * y_d01 # nolint
    y_t0 <- (cand == 1) * y_c00 + (cand == 3) * z * y_n10 + (cand == 3) * (1 - z) * y_n00 + (cand == 4) * y_d10 # nolint

    y <- ttt * y_t1 + (1 - ttt) * y_t0

    # Edit and print final simulated dataframe
    sim_dat <- data.frame(x_mat, rhos_mat, z, cand, ttt,
                          dp_c11, dp_a11, dp_a01, dp_d01,
                          dp_c00, dp_n10, dp_n00, dp_d10,
                          y_c11, y_a11, y_a01, y_d01,
                          y_c00, y_n10, y_n00, y_d10,
                          y_t1, y_t0, y) # nolint
    if (verbose) {
        cat("CAND parameters:\n")
        print(delta_mat)
        cat("\n Elementary potential outcomes parameters:\n")
        print(beta_mat)
    }
return(sim_dat)
}

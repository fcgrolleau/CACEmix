library(glmnet)
library(progress)

em_1_mon_fun <- function(dat, xvar_names, z_name, t_name,
                     maxit = 100, epsi = 1e-4,
                     bar = TRUE, verbose = TRUE) {

### EM algorithm to estimate conditional probabilities of c,a,n
### Here, assuming monotonicity (i.e., no defiers)

## Arguments
# dat: data frame
# xvar_names: names of the covariates
# z_name: name of the allocation variable
# t_name: name of the treatment variable
# maxit: maximum number of iterations
# epsi: stopping criteria
# bar: if TRUE, display a progress bar
# verbose: if TRUE, print the number of iterations

## Output
# res: list containing:
#      - (res$param_mat) a matrix of the estimated parameters
#      - (res$g_mat) the posterior probabilities

## Example
# res <- em_1_fun(dat, xvar_names = c("X0", "X1", "X2"),
#            z_name = "z", t_name = "ttt", maxit = 100, epsi = 1e-6,
#            bar = TRUE, verbose = TRUE)

x <- as.matrix(dat[xvar_names])
z <- dat[z_name]
ttt <- dat[t_name]

# initialize the prior probabilities associated with the nodes of the tree as
g_1 <- 1 / 3
g_2 <- 1 / 3
g_3 <- 1 / 3

# individual contributions to each expert’s likelihood
l_1 <- z * ttt + (1 - z) * (1 - ttt)
l_2 <- ttt
l_3 <- 1 - ttt

names(l_1) <- "C"
names(l_2) <- "A"
names(l_3) <- "N"

# initialize a vector of parameters
fitted_params <- rep(0, 3 * length(xvar_names))

# Create a progress bar
if (bar) {
    pb <- progress_bar$new(total = maxit)
}

# iterate until convergence
for (i in 1:maxit) {
old_params <- fitted_params

# E step:
# compute the posterior probabilities associated with the nodes of the tree
h_1 <- (g_1 * l_1) / (g_1 * l_1 + g_2 * l_2 + g_3 * l_3)
h_2 <- (g_2 * l_2) / (g_1 * l_1 + g_2 * l_2 + g_3 * l_3)
h_3 <- (g_3 * l_3) / (g_1 * l_1 + g_2 * l_2 + g_3 * l_3)

# M step
h_mat <- as.matrix(cbind(h_1, h_2, h_3))

fit <- tryCatch({
  glmnet(x, h_mat, family = "multinomial",
         lambda = 0, # no regularization
         alpha = 0, # no elastic net
         intercept = TRUE)
}, error = function(e) {
  warning(paste0("glmnet error at iteration ", i, " (monotonicity ie 3 classes)"))
  h_fix <- h_mat + matrix(runif(length(h_mat), 0, 1), nrow(h_mat), ncol(h_mat))
  h_fix <- t(apply(h_fix, 1, function(x) x/sum(x)))
  glmnet(x, h_fix, family = "multinomial",
         lambda = 0, # no regularization
         alpha = 0, # no elastic net
         intercept = TRUE)
})

g_mat <- drop(predict(fit, newx = x, type = "response"))
g_1 <- g_mat[, 1]
g_2 <- g_mat[, 2]
g_3 <- g_mat[, 3]

if (bar) {
    pb$tick()
}

param_mat <- sapply(1:3, function(i) coef(fit)[[i]]@x)

# check convergence on the squared frobenius norm
fitted_params <- as.vector(param_mat)

if (sum((fitted_params - old_params)^2) / length(fitted_params) < epsi) {
  if (verbose) {
    cat("\nEM no. 1 MON converged at iteration", i, "\n")
  }
  break
}
}
colnames(param_mat) <- c("C", "A", "N")
rownames(param_mat) <- xvar_names

res <- list(param_mat = t(param_mat), g_mat = g_mat)
return(res)
}

# load data simulation
#source("CACEmix/gen_dat.R")
#dat <- gen_data(n = 1000, p = 2, icpt = TRUE, sd_y = 1, eigeinval = NULL,
#         alpha = .5, bin = TRUE, seed = 42, verbose = TRUE)

###
#res <- em_1_mon_fun(dat, xvar_names = c("X0", "X1", "X2"),
#            z_name = "z", t_name = "ttt", maxit = 100, epsi = 1e-6)

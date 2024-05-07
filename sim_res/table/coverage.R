cover_fun <- function(res_sim) {
for (i in names(table(rownames(res_sim)))) {
    work_dat <-  res_sim[rownames(res_sim) == i, ]
    estimates <- work_dat[, "est"]
    #cis <- 2 * matrix(estimates, length(estimates), 2) - matrix(quantile(estimates, c(0.975, 0.025)), length(estimates), 2, byrow = TRUE)
    cis <- matrix(estimates, length(estimates), 2) + matrix(c(-1, 1) * qnorm(0.975) * sd(estimates), length(estimates), 2, byrow = TRUE)
    colnames(cis) <- c("lb", "ub")
    work_dat <- as.data.frame(cbind(work_dat, cis))
    cover <- mean(with(work_dat, true_cace > lb & true_cace < ub))
    len <- with(work_dat, ub - lb)[1]
    cat(i, "coverage =", cover, "length =", len, "\n")
}
cat("\n")
}

lapply(list(res_ori_n3, res_er_n3, res_mon_n3, res_mon_er_n3), cover_fun)

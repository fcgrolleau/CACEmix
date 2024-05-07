library(ggplot2)
library(latex2exp)
library(reshape2)
library(dplyr)

res_sim_fun_1 <- function(scenario, estimator) {

n_it <- length(unique(scenario[rownames(scenario) == estimator, "it"]))
n_samples <- scenario[rownames(scenario) == estimator, "n_samples"]

ate <- scenario[rownames(scenario) == estimator, "true_ate"]
ate_sample <- scenario[rownames(scenario) == estimator, "true_sample_ate"]
ate_ci <- ate + c(-1, 0, 1) * qnorm(.975) * sd(ate_sample)

true <- scenario[rownames(scenario) == estimator, "true_cace"]
true_sample <- scenario[rownames(scenario) == estimator, "true_sample_cace"]
true_ci <- true + c(-1, 0, 1) * qnorm(.975) * sd(true_sample)

est <- scenario[rownames(scenario) == estimator, "est"]

bias <- mean(est - true)
relative_bias <- bias / true[1]
variance <- mean((est - mean(est))^2)
mse <- bias^2 + variance

com_prop <- scenario[rownames(scenario) == estimator, "com_prop"][1] # nolint
com_prop_sample <- scenario[rownames(scenario) == estimator, "com_prop_sample"]  # nolint
com_prop_ci <- com_prop + c(-1, 0, 1) * qnorm(.975) * sd(com_prop_sample)

def_prop <- scenario[rownames(scenario) == estimator, "def_prop"][1] # nolint
def_prop_sample <- scenario[rownames(scenario) == estimator, "def_prop_sample"]  # nolint
def_prop_ci <- def_prop + c(-1, 0, 1) * qnorm(.975) * sd(def_prop_sample)

er0 <- scenario[rownames(scenario) == estimator, "er0"][1] # nolint
er0_sample <- scenario[rownames(scenario) == estimator, "er0_sample"] # nolint
er0_ci <- er0 + c(-1, 0, 1) * qnorm(.975) * sd(er0_sample)

er1 <- scenario[rownames(scenario) == estimator, "er1"][1] # nolint
er1_sample <- scenario[rownames(scenario) == estimator, "er1_sample"]  # nolint
er1_ci <- er1 + c(-1, 0, 1) * qnorm(.975) * sd(er1_sample)

res <- cbind(bias = bias,
                relative_bias = relative_bias,
                se = sqrt(variance),
                mse = mse,
                com_prop_lb = com_prop_ci[1],
                com_prop = com_prop_ci[2],
                com_prop_ub = com_prop_ci[3],
                def_prop_lb = def_prop_ci[1],
                def_prop = def_prop_ci[2],
                def_prop_ub = def_prop_ci[3],
                er0_lb = er0_ci[1],
                er0 = er0_ci[2],
                er0_ub = er0_ci[3],
                er1_lb = er1_ci[1],
                er1 = er1_ci[2],
                er1_ub = er1_ci[3],
                true_cace_lb = true_ci[1],
                true_cace = true_ci[2],
                true_cace_ub = true_ci[3],
                true_ate_lb = ate_ci[1],
                true_ate = ate_ci[2],
                true_ate_ub = ate_ci[3],
                n_samples = n_samples[1],
                n_it = n_it)

return(res)
}

res_sim_fun_2 <- function(scenario) {
temp <- lapply(c("no_mon_no_er_res",
                 "mon_no_er_res",
                 "no_mon_er_res",
                 "mon_er_res",
                 "wald_res",
                 "ps_matching"),
    function(x) res_sim_fun_1(scenario, x))

temp <- do.call(rbind, temp)

return(temp)
}

quick_fun <- function(scenario) {
    temp <- as.data.frame(res_sim_fun_2(get(scenario)))
    return(cbind(estimator = rownames(temp),
                 scenario = sub("_n.*", "", scenario),
                 misspecification = grepl("mis", scenario),
                 temp))
}

res_sim_fun_3 <- function(scenarios) {
temp <- lapply(scenarios, quick_fun)

do.call(rbind, temp)
}

plot_scenarios <- function(misspecification = FALSE, iv_res = TRUE) {

scenarios <- c(
 "res_er_n1",         "res_er_n1_mis",     "res_er_n2",
 "res_er_n2_mis",     "res_er_n3",         "res_er_n3_mis",
 "res_mon_er_n1",     "res_mon_er_n1_mis", "res_mon_er_n2",
 "res_mon_er_n2_mis", "res_mon_er_n3",     "res_mon_er_n3_mis",
 "res_mon_n1",        "res_mon_n1_mis",    "res_mon_n2",
 "res_mon_n2_mis",    "res_mon_n3",        "res_mon_n3_mis",
 "res_ori_n1",        "res_ori_n1_mis",    "res_ori_n2",
 "res_ori_n2_mis",    "res_ori_n3",        "res_ori_n3_mis")

res <- res_sim_fun_3(scenarios)
res$bias_sq <- res$bias^2

selection <- ifelse(rep(misspecification, length(res$misspecification)),
                  res$misspecification, !res$misspecification)

if (!iv_res) {
    selection <- selection & res$estimator != "wald_res" & res$estimator != "ps_matching" # nolint
}

proc_res <- res[selection,
                c("scenario", "estimator", "bias_sq", "mse", "n_samples")]

proc_res <- melt(proc_res,
            id.vars = setdiff(names(proc_res), c("bias_sq", "mse")),
            measure.vars = c("bias_sq", "mse"))

proc_res <- proc_res %>% arrange(desc(variable))

proc_res <- proc_res %>% mutate(scenario = recode(scenario,
                  "res_er" = "PI & ER", "res_mon" = "PI & MO", "res_mon_er" = "PI & MO & ER", "res_ori" = "PI"), # nolint
                        estimator = recode(estimator,
                  "mon_er_res" = "PI & MO & ER", "mon_no_er_res" = "PI & MO", "no_mon_er_res" = "PI & ER", "no_mon_no_er_res" = "PI"), # nolint
                  n_samples = paste0("n=", n_samples))

proc_res$n_samples <- factor(proc_res$n_samples, levels = unique(proc_res$n_samples)) # nolint

proc_res$value <- sqrt(proc_res$value)

gg_fig <- ggplot(proc_res, aes(fill = variable, y = value, x = estimator)) +
      facet_grid(scenario ~ n_samples) +
      scale_fill_manual("legend", values = c("mse" = "#00A1D5FF", "bias_sq" = "#00468BFF")) +  # nolint
      geom_bar(position = "identity", stat = "identity") +
      theme_light(base_size = 22) +
      theme(legend.position="none", strip.text = element_text(colour = 'black'), strip.background =element_rect(fill="gray90"), # nolint
            axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
            axis.title.y = element_text(face = "bold")) +
      xlab("") +
      ylab("Absolute bias / RMSE")

if (iv_res) {
gg_fig <- gg_fig + scale_x_discrete(labels =
            unname(TeX(c("$\\widehat{\\Delta}_{PI}$",
                         "$\\widehat{\\Delta}_{PI}^{ER}$",
                         "$\\widehat{\\Delta}_{PI,MO}$",
                         "$\\widehat{\\Delta}_{PI,MO}^{ER}$",
                         "$\\widehat{\\Delta}_{IV\\,matching}",
                         "$\\widehat{\\Delta}_{WALD}$"))))
} else {
gg_fig <- gg_fig + scale_x_discrete(labels =
            unname(TeX(c("$\\widehat{\\Delta}_{PI}$",
                         "$\\widehat{\\Delta}_{PI}^{ER}$",
                         "$\\widehat{\\Delta}_{PI,MO}$",
                         "$\\widehat{\\Delta}_{PI,MO}^{ER}$"))))
}
gg_fig
}

setwd("/Users/francois/Desktop/github repos/CACEmix/sim_res/reformated_from_2024_03_19_15_47_51_268704")
load("simulations_results.RData")

misspecification <- FALSE
p <- plot_scenarios(misspecification = misspecification, iv_res = TRUE)

ggsave(paste0("sim_res", ifelse(misspecification, "_mis", ""), ".pdf"),
            plot = p, device = "pdf", width = 10, height = 40 / 3)
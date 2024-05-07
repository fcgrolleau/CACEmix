library(dplyr)
library(ggplot2)
library(latex2exp)
source("/Users/francois/Desktop/github repos/CACEmix/sim_res/table/print_table.R")

setwd("/Users/francois/Desktop/github repos/CACEmix/sim_res/reformated_from_2024_03_19_15_47_51_268704")
load("simulations_results.RData") # nolint

misspecification  <- TRUE

log_vars2 <-prep_table(misspecification = misspecification) %>% filter(variable == "se") %>%  # nolint
  select(scenario, estimator, n_samples, value) %>%  # nolint
  mutate(value = log(value^2), n_samples= log(as.numeric(gsub("n=", "", as.character(n_samples)))), # nolint
         Estimator = estimator)

log_vars2$scenario <- factor(log_vars2$scenario)
log_vars2$Estimator <- factor(log_vars2$Estimator)

levels(log_vars2$scenario) <- c("Scenario 1:\nPrincipal Ignorability",
                                "Scenario 2:\nPrincipal Ignorability & Exclusion Restriction", # nolint
                                "Scenario 3:\nPrincipal Ignorability & Monotonicity", # nolint
                                "Scenario 4:\nPrincipal Ignorability & Monotonicity & Exclusion Restriction") # nolint

p2 <- ggplot(log_vars2, aes(x = n_samples, y = value, color = Estimator)) +
  facet_wrap("scenario") +
  geom_point(size = 3, alpha=.5) + geom_smooth(method='lm', se=FALSE, linetype="dotted", aes(color=Estimator)) + # nolint
  scale_color_discrete(labels = c("PI" = TeX("$\\widehat{\\Delta}_{PI}$"), # nolint
                                "PI & ER" = TeX("$\\widehat{\\Delta}_{PI}^{ER}$"), # nolint
                                "PI & MO" = TeX("$\\widehat{\\Delta}_{PI,MO}$"), # nolint
                                "PI & MO & ER" = TeX("$\\widehat{\\Delta}_{PI,MO}^{ER}$"), # nolint
                                "ps_matching" = TeX("$\\widehat{\\Delta}_{IV\\,matching}$"), # nolint
                                "wald_res" = TeX("$\\widehat{\\Delta}_{WALD}$"))) + # nolint
  labs(color = "Estimator") +
  xlab(TeX("ln $n$")) +
  ylab(TeX("ln var $(\\hat{\\Delta})$")) +
  coord_cartesian(xlim = c(7.5, 9.5), ylim = c(-8.5, -4.5)) +
  theme_minimal() + theme(
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 12),
    axis.title = element_text(size = 25, face = "bold"),
    legend.background = element_rect(fill = "white"),
    legend.key.size = unit(.9, "cm"),
    legend.position = c(.9, .3),
    legend.text.align = 0)
p2

# Save the plot to a square PDF

ggsave(paste0("convergence_plot", ifelse(misspecification, "_mis", ""), ".pdf"),
      plot = p2, device = "pdf", width = 10, height = 10)
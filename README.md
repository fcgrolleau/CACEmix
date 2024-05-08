# CACEmix
 
<!-- badges: start -->
[![license](https://img.shields.io/badge/license-MIT-blue)](https://github.com/fcgrolleau/Dynamic-RRT/blob/main/LICENSE)
[![R badge](https://img.shields.io/badge/Build%20with-%20R,%20♥%20and%20Python-blue)](https://rstudio.github.io/reticulate/index.html)
[![R 4.2.1](https://img.shields.io/badge/R-4.3.3-blue.svg)](https://www.r-project.org) 
<!-- badges: end -->

This repository reproduces results from the paper *Estimating Complier Average Causal Effects with Mixtures of Experts* [<a href="https://arxiv.org/pdf/2405.02779">*arXiv*</a>].

### Authors
This repository is written and maintained by François Grolleau (grolleau@stanford.edu).

### Reproducibility

`em_part1.R` implements the EM Algorithm C.1 from the paper (step 1.)

`em_part1_mon.R` implements the EM Algorithm C.12 from the paper (step 1 with monotonicity.)

`em_part2_bin_er.R` implements the EM Algorithm C.7 from the paper (step 2 with exclusion restriction.)

`em_part2_bin_test.R` implements the EM Algorithm C.2 from the paper (step 2.)

`full_fun_er_fix2.R` implements the EM Algorithm C.7 from the paper (step 2 with exclusion restriction.)

`gen_dat_full.R` data generating process from the simulations.

`plot_sim.R` plots figures 1 and 4 from simulation results.

`simulation_funs.r` functions to prepare the parallelized simulation.

`simulations.R` code to run the simulations and save the results.

`sim_res/table/convergence.R` reproduces the convergence plots given figures 2 and 5.

Miscellananeous: `full_fun_er_fix.R`, `full_fun_er_fix_eff_wrong.R`

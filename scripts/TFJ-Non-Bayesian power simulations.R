library(ggplot2)
library(tictoc)
library(mvtnorm)
library(foreach)
library(doParallel)
library(future.apply)

source("../scripts/functions.R")

make_subj_cov = function(SD) {
  matrix(
    c(SD^2, rep(0, 15)), 
    nrow = 4,
    dimnames = list(
      c("Intercept", "TS.vs.CNTL", "MT.vs.CNTL", "ST.vs.CNTL"), 
      c("Intercept", "TS.vs.CNTL", "MT.vs.CNTL", "ST.vs.CNTL")
    )
  )
  
}

get_significance = function(d) {
  with(d, 
       rbind(
         round(prop.table(table(TS.vs.CNTL.z > 1.96)), 3) * 100,
         round(prop.table(table(MT.vs.CNTL.z > 1.96)), 3) * 100,
         round(prop.table(table(ST.vs.CNTL.z > 1.96)), 3) * 100,
         round(prop.table(table(TS.vs.MT.z > 1.96)), 3) * 100,
         round(prop.table(table(MT.vs.ST.z > 1.96)), 3) * 100
       )[,2])
}

get_significant_convergence_failure = function(d) {
  with(d, 
       rbind(
         round(prop.table(table(TS.vs.CNTL.z > 1.96 & ConvergenceFailure.treat)), 3) * 100,
         round(prop.table(table(MT.vs.CNTL.z > 1.96 & ConvergenceFailure.treat)), 3) * 100,
         round(prop.table(table(ST.vs.CNTL.z > 1.96 & ConvergenceFailure.treat)), 3) * 100,
         round(prop.table(table(TS.vs.MT.z > 1.96 & ConvergenceFailure.diff)), 3) * 100,
         round(prop.table(table(MT.vs.ST.z > 1.96 & ConvergenceFailure.diff)), 3) * 100
       )[,2])
}

get_significance_and_convergence = function(d) {
  paste0(
    get_significance(d),
    "\% (",
    get_significant_convergence_failure(d),
    "\%)"
  )
}


# Run simulations if they have not already been stored.
if (!file.exists("../models/powersims.RData")) {
  # Treatment-coded data generation
  e = cbind("Intercept" = 1.71, "TS.vs.CNTL" = .82, "MT.vs.CNTL" = .51, "ST.vs.CNTL" = .30)

  # Make variance-covariance matrices of random effects
  subj_SD = .65
  subj.ranef.covar = make_subj_cov(subj_SD)
  item.sd = c(1.13, .44, .14, .42)
  item.cor = matrix(
    c(1, .07, .09, -.24, .07, 1, .4, .63, 
      .09, 0.4, 1, .27, -.24, .63, .27, 1), 
    nrow = length(e),
    dimnames = list(colnames(e), colnames(e))
  ) 
  item.ranef.covar = MBESS::cor2cov(cor.mat = item.cor, sd = item.sd)
  
  my.simulator <- make.data.generator(
    estimates = e, 
    n.subj.perCondition = 80, 
    n.obs = 51, 
    subject.ranef.covar = subj.ranef.covar,
    item.ranef.covar = item.ranef.covar
  )
  
  my.BB08.simulator <- make.data.generator(
    estimates = e, 
    n.subj.perCondition = c(10, 40, 10, 10), 
    n.obs = 51, 
    subject.ranef.covar = subj.ranef.covar,
    item.ranef.covar = item.ranef.covar
  )
  
  my.half.simulator <- make.data.generator(
    estimates = e / 2, 
    n.subj.perCondition = 80, 
    n.obs = 51, 
    subject.ranef.covar = subj.ranef.covar,
    item.ranef.covar = item.ranef.covar
  )
  
  # Don't change order of data generator creation
  item.sd = item.sd * 2
  item.ranef.covar = MBESS::cor2cov(cor.mat = item.cor, sd = item.sd)
  subject.ranef.covar = make_subj_cov(subj_SD * 2)
  
  my.2var.simulator <- make.data.generator(
    estimates = e, 
    n.subj.perCondition = 80, 
    n.obs = 51, 
    subject.ranef.covar = subj.ranef.covar,
    item.ranef.covar = item.ranef.covar
  )
  
  num.sims = 1000
  n.workers = 12
  plan(multisession, workers = n.workers)

  tic()
  d.sim.small = future_replicate(n = num.sims, workers = n.workers,
                                 expr = fit.models(my.BB08.simulator()), simplify = "matrix")
  toc()
  d.sim.small = 
    data.frame(t(d.sim.small)) %>% 
    as_tibble() %>% 
    mutate_all(.funs = unlist)
  
  tic()
  d.sim = future_replicate(n = num.sims, workers = n.workers,
                           expr = fit.models(my.simulator()), simplify = "matrix")
  toc()
  d.sim = 
    data.frame(t(d.sim)) %>% 
    as_tibble() %>% 
    mutate_all(.funs = unlist)
  
  tic()
  d.sim.half = future_replicate(n = num.sims, workers = n.workers,
                                expr = fit.models(my.half.simulator()), simplify = "matrix")
  toc()
  d.sim.half = 
    data.frame(t(d.sim.half)) %>% 
    as_tibble() %>% 
    mutate_all(.funs = unlist)
  
  tic()
  d.sim.2var = future_replicate(n = num.sims, workers = n.workers,
                                expr = fit.models(my.2var.simulator()), simplify = "matrix")
  toc()
  d.sim.2var = 
    data.frame(t(d.sim.2var)) %>% 
    as_tibble() %>% 
    mutate_all(.funs = unlist)
  
  # Save power simulations
  save(d.sim, d.sim.small, d.sim.half, d.sim.2var, file = "../models/powersims.RData", compress = T)
  
  # Go back to sequential processing
  plan(sequential)
} else {
  cat("\nLoading existing file with power simulations.\nTo re-run simulations, delete powersims.RData.\n\n")
  load("../models/powersims.RData")
}

# Power = proportion of significant effects in expected direction (for which model also converged)
library(knitr)
k = kable(
  matrix(cbind(
    get_significance(d.sim.small),
    get_significance(d.sim),
    get_significance(d.sim.2var),
    get_significance(d.sim.half)
  ), nrow = 3, ncol = 5, 
  dimnames = list(
    c("TS vs. CNTL", "MT vs. CNTL", "ST vs. CNTL", "TS vs. MT", "MT vs. ST"),
    col.names = c("BB08", "Exp 1", "Exp 1 (2 * \sigma)", "Exp 1 (.5 \beta)")
  )),
  format = "latex"
)
k



## ----------------------------------------------------------------------------
#
# PLOTS
#
## ----------------------------------------------------------------------------
myGplot.defaults(type = "paper")
ggplot(d.sim %>% as_tibble() %>% mutate_all(.funs = unlist), 
       aes(x = TS.vs.MT.z, y = MT.vs.ST.z)) +
  geom_vline(xintercept = c(-1.96, 1.96), color = "red", linetype = 2) +
  geom_hline(yintercept = c(-1.96, 1.96), color = "red", linetype = 2) +
  geom_rect(xmin = -1.96,
            xmax = 1.96,
            ymin = -1.96,
            ymax = 1.96, fill = "red", alpha = .01) +
  geom_point(alpha = .5) +
  geom_density_2d() +
  theme_bw() +
  scale_x_continuous("TS vs. MT\n(z-value)",
                     limits = c(-15,15)) +
  scale_y_continuous("MT vs. ST\n(z-value)",
                     limits = c(-10,10))
ggsave("figures/power.pdf", height = 4, width = 4)


ggplot(d.sim %>% as_tibble() %>% mutate_all(.funs = unlist), 
       aes(x = TS.vs.MT.estimate, y = MT.vs.ST.estimate)) +
  geom_point(aes(shape = ifelse(abs(TS.vs.MT.z) > 1.96 & abs(MT.vs.ST.z) > 1.96, "yes", "no"),
                 alpha = ifelse(abs(TS.vs.MT.z) > 1.96 & abs(MT.vs.ST.z) > 1.96, 1, .7))) +
  geom_point(x = mean(d.sim$TS.vs.MT.estimate), 
             y = mean(d.sim$MT.vs.ST.estimate),
             color = "blue",
             size = 2
  ) +
  geom_density_2d() +
  geom_vline(xintercept = 0, color = "black", linetype = 2) +
  geom_hline(yintercept = 0, color = "black", linetype = 2) + 
  theme_bw() +
  scale_shape_discrete("Condition significant?") +
  guides(alpha = "none")


d.sim.old = d.sim


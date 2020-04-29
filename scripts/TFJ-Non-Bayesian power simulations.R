library(ggplot2)
library(tictoc)
library(mvtnorm)
library(foreach)
library(doParallel)

source("../scripts/functions.R")

e = cbind("Intercept" = 2.07, "TS.vs.MT" = .33, "MT.vs.ST" = .16, "ST.vs.CNTL" = .33)
item.sd = c(1.1, .26, .22, .4)
item.cor = matrix(
   c(1, .17, .2, 0, .17, 1, -.29, .54, .2, -.29, 1, -.53, 0, .54, -.53, 1), 
   nrow = length(e),
   dimnames = list(colnames(e), colnames(e))
)
item.ranef.covar = MBESS::cor2cov(cor.mat = item.cor, sd = item.sd)

my.simulator <- make.data.generator(
  estimates = e, 
  n.subj.perCondition = 80, 
  n.obs = 51, # Include all five blocks of 7 stims during Test block
  subject.ranef.covar = matrix(
    c(.58^2, rep(0, length(e)**2-1)), 
    nrow = length(e),
    dimnames = list(colnames(e), colnames(e))
  ),
  item.ranef.covar = item.ranef.covar
)

my.small.simulator <- make.data.generator(
  estimates = e, 
  n.subj.perCondition = 10, 
  n.obs = 51, # Include all five blocks of 7 stims during Test block
  subject.ranef.covar = matrix(
    c(.58^2, rep(0, length(e)**2-1)), 
    nrow = length(e),
    dimnames = list(colnames(e), colnames(e))
  ),
  item.ranef.covar = item.ranef.covar
)

my.half.simulator <- make.data.generator(
  estimates = cbind("Intercept" = 2.07, "TS.vs.MT" = .165, "MT.vs.ST" = .08, "ST.vs.CNTL" = .165), 
  n.subj.perCondition = 80, 
  n.obs = 51, # Include all five blocks of 7 stims during Test block
  subject.ranef.covar = matrix(
    c(.58^2, rep(0, 15)), 
    nrow = 4,
    dimnames = list(
      c("Intercept", "TS.vs.MT", "MT.vs.ST", "ST.vs.CNTL"), 
      c("Intercept", "TS.vs.MT", "MT.vs.ST", "ST.vs.CNTL")
    )
  ),
  item.ranef.covar = item.ranef.covar
)

item.sd = (c(1.1, .26, .22, .4)^2*2)^.5
item.cor = matrix(
  c(1, .17, .2, 0, .17, 1, -.29, .54, .2, -.29, 1, -.53, 0, .54, -.53, 1), 
  nrow = length(e),
  dimnames = list(colnames(e), colnames(e))
)
item.ranef.covar = MBESS::cor2cov(cor.mat = item.cor, sd = item.sd)
my.2var.simulator <- make.data.generator(
  estimates = cbind("Intercept" = 2.07, "TS.vs.MT" = .33, "MT.vs.ST" = .16, "ST.vs.CNTL" = .33), 
  n.subj.perCondition = 80, 
  n.obs = 51, # Include all five blocks of 7 stims during Test block
  subject.ranef.covar = matrix(
    c((.58^2)*2, rep(0, 15)), 
    nrow = 4,
    dimnames = list(
      c("Intercept", "TS.vs.MT", "MT.vs.ST", "ST.vs.CNTL"), 
      c("Intercept", "TS.vs.MT", "MT.vs.ST", "ST.vs.CNTL")
    )
  ),
  item.ranef.covar = item.ranef.covar
)

#library(doSNOW)
# registerDoSNOW(cl)
# pb <- txtProgressBar(max = num.sims, style = 3)
# progress <- function(n) setTxtProgressBar(pb, n)
# opts <- list(progress = progress)

# Initiate cluster
# cl<-parallel::makeCluster(10)
# doParallel::registerDoParallel(cl)
# num.sims = 1000
# 
# tic()
# d.sim.small = foreach(n = 1:num.sims,
#         .verbose = T,
#         .combine = rbind,
#         .final = as.data.frame,
#         .packages = c("mvtnorm")) %dopar% fit.models(my.small.simulator())
#           
# toc()
# parallel::stopCluster(cl)

library(future.apply)
plan(multiprocess)
d.sim = future_replicate(n = num.sims, workers = 10,
                expr = fit.models(my.simulator()), simplify = "matrix")
d.sim = data.frame(t(d.sim))

d.sim.half = future_replicate(n = num.sims, workers = 10,
                         expr = fit.models(my.half.simulator()), simplify = "matrix")
d.sim.half = data.frame(t(d.sim.half))

d.sim.2var = future_replicate(n = num.sims, workers = 10,
                              expr = fit.models(my.2var.simulator()), simplify = "matrix")
d.sim.2var = data.frame(t(d.sim.2var))
saveRDS(list(d.sim, d.sim.small, d.sim.half, d.sim.2var), "../models/powersims.RDS", compress = T)


l = readRDS("../models/powersims.RDS")
d.sim = l[[1]]
d.sim.small = l[[2]]
d.sim.half = l[[3]]
d.sim.2var = l[[4]]

summary(d.sim)
summary(d.sim.small)
summary(d.sim.half)
summary(d.sim.2var)

# Power = proportion of "true" (significant effects)
library(knitr)
k = kable(
  matrix(cbind(
    rbind(
      round(prop.table(table(d.sim.small$TS.vs.MT.z > 1.96)), 3) * 100,
      round(prop.table(table(d.sim.small$MT.vs.ST.z > 1.96)), 3) * 100,
      round(prop.table(table(d.sim.small$ST.vs.CNTL.z > 1.96)), 3) * 100
    )[,2],
    rbind(
      round(prop.table(table(d.sim$TS.vs.MT.z > 1.96)), 3) * 100,
      round(prop.table(table(d.sim$MT.vs.ST.z > 1.96)), 3) * 100,
      round(prop.table(table(d.sim$ST.vs.CNTL.z > 1.96)), 3) * 100
    )[,2],
    rbind(
      round(prop.table(table(d.sim.half$TS.vs.MT.z > 1.96)), 3) * 100,
      round(prop.table(table(d.sim.half$MT.vs.ST.z > 1.96)), 3) * 100,
      round(prop.table(table(d.sim.half$ST.vs.CNTL.z > 1.96)), 3) * 100
    )[,2],
    rbind(
      round(prop.table(table(d.sim.2var$TS.vs.MT.z > 1.96)), 3) * 100,
      round(prop.table(table(d.sim.2var$MT.vs.ST.z > 1.96)), 3) * 100,
      round(prop.table(table(d.sim.2var$ST.vs.CNTL.z > 1.96)), 3) * 100
    )[,2]
  ), nrow = 3, ncol = 4, 
  dimnames = list(
    c("TS vs. MT", "MT vs. ST", "ST vs. CNTL"),
    col.names = c("BB08", "Exp 1", "Exp 1 (.5 \beta)", "Exp 1 (2* sigma)")
  )),
  format = "latex"
)
## ----------------------------------------------------------------------------
#
# PLOTS
#
## ----------------------------------------------------------------------------
myGplot.defaults(type = "paper")
ggplot(d.sim, 
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


ggplot(d.sim, 
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


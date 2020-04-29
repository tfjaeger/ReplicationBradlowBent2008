## -------------------------------------------------------------------------------
# If you use this or other scripts for our replication analyses, we would appreciate 
# if you cite the article as the sources. Thank you.
## -------------------------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(brms)
library(bridgesampling)

setwd("~/Box Sync/_Papers - Box/Xie, Liu, Jaeger - Foreign-accented speech/scripts")

## ----------------------------------------------------------------------------------
# Functions
## ----------------------------------------------------------------------------------
source("../scripts/functions.R")

data.generator <- function(
  estimates = cbind("Intercept" = 0, "TS" = 0), 
  n.subj.perCondition = 80,
  n.subj = n.subj.perCondition * 2, # Total number of subjects (for between subject design)
  n.obs = 80,                        # Number of items per subject
  subject.ranef.covar = matrix(
    c(.4, rep(0, length(estimates)**2-1)),
    nrow = length(estimates),
    dimnames = list(colnames(estimates), colnames(estimates))
  ),
  item.ranef.covar = matrix(
    c(.4, 0, 0, .4),
    nrow = length(estimates),
    dimnames = list(colnames(estimates), colnames(estimates))
  )
) {
  require(mvtnorm)
  require(magrittr)
  require(dplyr)
  
  n = n.obs * n.subj 
  d = expand.grid(Subject = 1:n.subj,
                  Item = 1:n.obs,
                  Intercept = 1)
  d = d[order(d$Subject),]
  d$Cond = factor(
    sort(rep(
      c("Control", "Talker-specific"),
      n / 2)),
    levels = c("Talker-specific", "Control"))
  contrasts(d$Cond) = cbind("TS" = c(1, 0))
  
  d %<>%
    mutate(
      TS = case_when(
        Cond == "Talker-specific" ~ 1,
        T ~ 0
      )
    )
  
  generate.data <- function() {
    if (all(dim(item.ranef.covar) == 1)) {
      item.adjustment <- data.frame(sapply(1:n.obs, function(x) { rnorm(n = 1, mean = 0, eval(item.ranef.covar**.5)) } ))
    } else {
      item.adjustment <- as.data.frame(t(sapply(1:n.obs, function(x) t(rmvnorm(n = 1, sigma = eval(item.ranef.covar))))))
    }
    names(item.adjustment) <- colnames(eval(item.ranef.covar))
    item.adjustment$Item <- 1:n.obs
    
    if (all(dim(subject.ranef.covar) == 1)) {
      subject.adjustment <- data.frame(sapply(1:n.subj, function(x) { rnorm(n = 1, mean = 0, subject.ranef.covar**.5) } ))
    } else{
      subject.adjustment <- as.data.frame(t(sapply(1:n.subj, function(x) t(rmvnorm(n = 1, sigma = subject.ranef.covar)))))
    }
    names(subject.adjustment) <- colnames(subject.ranef.covar)
    subject.adjustment$Subject <- 1:n.subj
    
    d$y <- rowSums(d[,colnames(estimates)] * (
      matrix(
        rep(as.numeric(estimates), n),
        byrow = T, 
        nrow = n
      ) +
        item.adjustment[d$Item, colnames(estimates)] + 
        subject.adjustment[d$Subject, colnames(estimates)]
    ))
    d$Item = factor(d$Item)
    d$Subject = factor(d$Subject)
    
    d$Outcome = sapply(d$y, FUN = function(x) { rbinom(1,1,plogis(x))})
    
    return(d)
  }
  
  return(generate.data)
}



fit.replication.models <- function(data1, data2, 
                                   model.orig, model.rep.full, model.rep.reduced, 
                                   niter = 11000, warmup = 1000, chains = 8, cores = chains, 
                                   save_models = F, print_priors = T,
                                   n.bs.rep = 1) {
  require(brms)
  require(bridgesampling)
  
  if (save_models) 
    message("NB: Models will be saved to disk. This can take up a lot of space. (To change this, set save_models = FALSE.)\n")

  message(paste0("Sampling with ", niter, " samples per chain for ", chains, " chains (warmup=", warmup, ").\n"))
  message(paste0("Using ", cores, " cores.\n"))
    
  cat("Fitting model to original data ...\n")  
  b.orig <- update(
    model.orig,
    newdata = data1, 
    iter = niter, warmup = warmup, chains = chains, cores = cores,
    refresh = 0, silent =T,
    recompile = FALSE,
    file = if (save_models) paste0("../models/tmp/m_orig_full_", Sys.time()) else NULL
  )
  
  # NB: This has been carefully checked to create output in exactly the right order for the my.stanvars
  #     IF THE STANVARS ABOVE CHANGE, THIS PIECE OF CODE MIGHT ALSO NEED TO CHANGE!
  # NB: Currently, the script removes the correlation information since brm only accept the lkj prior for
  #     correlation matrices.
  p = posterior_samples(b.orig) %>%
      dplyr::select(-matches("^(r|lp|cor)_.*")) %<>%
      summarise_all(
        .funs = c("mean", "sd")
      ) %>%
    gather() 
    
  my.rep.stanvars = 
    stanvar(p[1, "value"], name = "mean_intercept") +
    stanvar(p[2, "value"], name = "mean_condition") +
    stanvar(p[3, "value"], name = "mean_intercept_byItem") +
    stanvar(p[4, "value"], name = "mean_condition_byItem") +
    stanvar(p[5, "value"], name = "mean_intercept_bySubject") +
    stanvar(p[6, "value"], name = "sd_intercept") +
    stanvar(p[7, "value"], name = "sd_condition") +
    stanvar(p[8, "value"], name = "sd_intercept_byItem") +
    stanvar(p[9, "value"], name = "sd_condition_byItem") +
    stanvar(p[10, "value"], name = "sd_intercept_bySubject")
  
  if (print_priors) {
    cat("Priors for replication data:\n")  
    print(p)
  }
  
  cat("Fitting model to replication data ...\n")  
  b.rep <- update(
    model.rep.full,
    newdata = data2, 
    stanvars = my.rep.stanvars,
    iter = niter, warmup = warmup, chains = chains, cores = cores,
    save_all_pars = T,
    refresh = 0, silent =T,
    recompile = FALSE,
    file = if (save_models) paste0("../models/tmp/m_rep_full_", Sys.time()) else NULL
  )
  
  cat("Fitting reduced model to replication data ...\n")  
  b.rep_reduced <- update(
    model.rep.reduced,
    newdata = data2, 
    stanvars = my.rep.stanvars,
    iter = niter, warmup = warmup, chains = chains, cores = cores,
    save_all_pars = T,
    refresh = 0, silent =T,
    recompile = FALSE,
    file = if (save_models) paste0("../models/tmp/m_rep_reduced_", Sys.time()) else NULL
  )
  
  simulation = NA
  
  cat("Bridgesampling full model ...\n")
  if (file.exists("../models/tmp/bs.rep_full.rds")) {
    bs.rep = readRDS("../models/tmp/bs.rep_full.rds")
  } else {
    bs.rep = bridge_sampler(samples = b.rep, cores = 8, method = "warp3", maxiter = 1000, repetitions = n.bs.rep)
    saveRDS(bs.rep, file = "../models/tmp/bs.rep_full.rds")
  }

  cat("Bridgesampling reduced model ...\n")
  if (file.exists("../models/tmp/bs.rep_reduced.rds")) {
    bs.rep_reduced = readRDS("../models/tmp/bs.rep_reduced.rds")
  } else {
    bs.rep_reduced = bridge_sampler(samples = b.rep_reduced, cores = 8, method = "warp3", maxiter = 1000, repetitions = n.bs.rep)
    saveRDS(bs.rep_reduced, file = "../models/tmp/bs.rep_reduced.rds")
  }

  cat("Calculating replication BF and storing information about this analysis ...\n")
  s.orig = summary(b.orig)$fixed
  s.rep = summary(b.rep)$fixed
  simulation <- data.frame(
    Int.1.estimate = s.orig[1,1],
    Int.1.lower = s.orig[1,3],
    Int.1.upper = s.orig[1,4],
    Cond.1.estimate = s.orig[2,1],
    Cond.1.lower = s.orig[2,3],
    Cond.1.upper = s.orig[2,4],
    Cond.1.hypothesis.p = hypothesis(b.orig, "CondTS > 0")$hypothesis$Post.Prob,
    Int.2.estimate = s.rep[1,1],
    Int.2.lower = s.rep[1,3],
    Int.2.upper = s.rep[1,4],
    Cond.2.estimate = s.rep[2,1],
    Cond.2.lower = s.rep[2,3],
    Cond.2.upper = s.rep[2,4],
    Cond.2.hypothesis.p = hypothesis(b.rep, "CondTS > 0")$hypothesis$Post.Prob,
    BF = bf(bs.rep, bs.rep_reduced)$bf
  )

  # store simulation results
  if (file.exists("../models/tmp/simulation.rds"))
    simulation = rbind(readRDS("../models/tmp/simulation.rds"), simulation)

  saveRDS(simulation, "../models/tmp/simulation.rds")

  # delete all temporary model files
  file.remove(list.files(path = "../models/tmp", pattern = "b.*.rds", full.names = T))
  rm(b.orig, b.rep, b.rep_reduced, bs.rep, bs.rep_reduced)
  
  return(simulation)
}
  

## ----------------------------------------------------------------------------------
# Setting up data generators for original and replication data
## ----------------------------------------------------------------------------------
nsubj = 80
nitem = 51
estimates = cbind("Intercept" = NA, "TS" = NA)
my.data.orig = data.generator(
  estimates = cbind("Intercept" = 1.66, "TS" = .49),
  n.subj.perCondition = nsubj,
  n.obs = nitem,
  subject.ranef.covar = matrix(
    c(.58, rep(0, length(estimates)**2-1)),
    nrow = length(estimates),
    dimnames = list(colnames(estimates), colnames(estimates))
  ),
  item.ranef.covar = matrix(
    c(1.1, 0, 0, .13),
    nrow = length(estimates),
    dimnames = list(colnames(estimates), colnames(estimates))
  )
)

my.data.rep = data.generator(
  estimates = cbind("Intercept" = 1.66, "TS" = .49),
  n.subj.perCondition = nsubj,
  n.obs = nitem,
  subject.ranef.covar = matrix(
    c(.58, rep(0, length(estimates)**2-1)),
    nrow = length(estimates),
    dimnames = list(colnames(estimates), colnames(estimates))
  ),
  item.ranef.covar = matrix(
    c(1.1, 0, 0, .13),
    nrow = length(estimates),
    dimnames = list(colnames(estimates), colnames(estimates))
  )
)


## ----------------------------------------------------------------------------------
# Setting up compilation of models
## ----------------------------------------------------------------------------------
formula.full = formula(Outcome ~ 1 + Cond + (1 | Subject) + (1 + Cond | Item))
formula.reduced = formula(Outcome ~ 1 + (1 | Subject) + (1 + Cond | Item))

# Setting up the values for priors
my.stanvars = 
  stanvar(0, name = "mean_intercept") +
  stanvar(0, name = "mean_condition") + 
  stanvar(0, name = "mean_intercept_byItem") +
  stanvar(0, name = "mean_condition_byItem") +
  stanvar(0, name = "mean_intercept_bySubject") +
  stanvar(1, name = "sd_intercept") +
  stanvar(1, name = "sd_condition") +
  stanvar(1, name = "sd_intercept_byItem") +
  stanvar(1, name = "sd_condition_byItem") +
  stanvar(1, name = "sd_intercept_bySubject") 
  
# Setting up priors
# for original model
my.priors = c(
  prior(normal(0, 10), class = b),
  prior(cauchy(0, 2), class = sd),
  prior(lkj(1), class = cor))

# for reduced replication model
my.reducedpriors = c(
  prior(normal(mean_intercept, sd_intercept), class = Intercept),
  prior(normal(mean_intercept_bySubject, sd_intercept_bySubject), class = sd, coef = Intercept, group = Subject),
  prior(normal(mean_intercept_byItem, sd_intercept_byItem), class = sd, coef = Intercept, group = Item),
  prior(normal(mean_condition_byItem, sd_condition_byItem), class = sd, coef = CondTS, group = Item),
  prior(lkj(1), class = cor)
)

# Compiling and setting up DSO for reduced model
model.rep.reduced = brm(
  formula = formula.reduced,  
  family = bernoulli(),
  data = my.data.rep(),
  # Setting up priors with stanvars so that they specific values can be manipulated within the loop
  prior = my.reducedpriors,
  stanvars = my.stanvars,
  iter = 2, chains = 1,
  save_dso = T, save_model = "../models/tmp/m_sim_rep_reducedmodel.stan", file = "../models/tmp/m_sim_rep_reducedmodel"
)

# Compiling and setting up DSO for full model
model.rep.full = brm(
  formula = formula.full,  
  family = bernoulli(),
  data = my.data.rep(),
  # Setting up priors with stanvars so that they specific values can be manipulated within the loop
  prior = my.reducedpriors + prior(normal(mean_condition, sd_condition), class = b, coef = CondTS),
  stanvars = my.stanvars,
  iter = 2, chains = 1,
  save_dso = T, save_model = "../models/tmp/m_sim_rep_fullmodel.stan", file = "../models/tmp/m_sim_rep_fullmodel"
)

# Compiling and setting up DSO for original model
model.orig = brm(
  formula = formula.full,  
  family = bernoulli(),
  data = my.data.orig(),
  # Using standard uninformative weakly regularizing priors
  prior = my.priors,
  iter = 2, chains = 1,
  save_dso = T, save_model = "../models/tmp/m_sim_orig_fullmodel.stan", file = "../models/tmp/m_sim_orig_fullmodel"
)


## ----------------------------------------------------------------------------------
# Running loop
## ----------------------------------------------------------------------------------
num.sims = 5
tictoc::tic()
d.sim.orig_rep = plyr::rdply(.n = num.sims, 
                             fit.replication.models(
                               data1 = my.data.orig(), 
                               data2 = my.data.rep(), 
                               model.orig = model.orig, 
                               model.rep.full = model.rep.full, 
                               model.rep.reduced = model.rep.reduced,
                               niter = 11000, warmup = 1000),
                             .progress = "tk")
tictoc::toc()
saveRDS(d.sim.orig_rep, "../models/replication_simulations_orig and rep-MTvsCNTL_same intercept.RDS", compress = T)



my.test = function(data1, data2) {
  data1 %>% 
    group_by(Cond) %>%
    summarise(meanByCondition1 = qlogis(mean(Outcome))) %>%
    left_join(
      data2 %>% 
        group_by(Cond) %>%
        summarise(meanByCondition2 = qlogis(mean(Outcome))),
      by = "Cond"
    )
}

d.test = plyr::rdply(.n = 10000,
            my.test(my.data.orig(), my.data.rep())
)

d.test %>%
  filter(Cond == "Control") %>%
  summarise(DIFF = mean(meanByCondition1 - meanByCondition2))


## -------------------------------------------------------------------------------
# This replication analysis was also incorporated into the R Markdown report that 
# now constitutes the supplementary information. The present file contains the 
# original scripts used to conduct the replication analysis, in case they prove 
# useful to other researchers. If you develop replication analyses from these 
# scripts, we would appreciate if you cite the article as the sources.
## -------------------------------------------------------------------------------
library(dplyr)
library(tibble)
library(magrittr)
library(tidyr)
library(brms)
library(bridgesampling)
library(shinystan)
library(bayesplot)
library(actuar)
library(fitdistrplus)

source("../scripts/functions.R")
shape = function(x, ...) { 
  if (sum(1/x > 0) == length(x)) {
    fitdistrplus::fitdist(x, "invgamma")$estimate["shape"]
  } else 0 
}
scale = function(x, ...) { 
  if (sum(1/x > 0) == length(x)) {
    fitdistrplus::fitdist(x, "invgamma")$estimate["scale"]
  } else 0 
}


# load data
d.exp1ab.training = readRDS("../data/Exp1ab.training postRMD.RDS")
d.exp1ab.test = readRDS("../data/Exp1ab.test postRMD.RDS")

myGplot.defaults(type = "paper")
d.exp1ab.test %>%
  filter(Experiment == "1b") %>%
  group_by(WorkerID, Sentence) %>%
  distinct(.keep_all = TRUE) %>%
  group_by(WorkerID, Condition2) %>%
  dplyr::summarise(PropKeywordsCorrect = mean(PropKeywordsCorrect)) %>%
  ggplot(aes(x = Condition2, y = PropKeywordsCorrect, colour = Condition2, fill = Condition2)) +
  stat_summary(fun.y = "mean", geom = "point", alpha = 1, size= 1.1, stroke=3) +
  geom_dotplot(binaxis = "y", binwidth = 0.015, stackdir = "center", alpha = .3, position = position_dodge(.9)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, linetype=1, size =1)+ 
  xlab("Exposure condition") +
  scale_y_continuous("Proportion of keywords correctly transcribed") +
  scale_colour_manual("Exposure Condition", values = setNames(c("#00BFC4", "#F8766D", "#7CAE00", "#C77CFF"), c("Talker-specific", "Control", "Single talker", "Multi-talker")), guide=FALSE) +
  scale_fill_manual("Exposure Condition", values = setNames(c("#00BFC4", "#F8766D", "#7CAE00", "#C77CFF"), c("Talker-specific", "Control", "Single talker", "Multi-talker")), guide=FALSE) +
  coord_cartesian(ylim = c(0.3,1.02)) +
  theme(legend.position = "none", panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Step 1: set priors and contrasts
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
my.priors = c(
  prior(normal(0, 10), class = b),
  prior(cauchy(0, 2), class = sd),
  prior(lkj(1), class = cor)
)

defineContrasts = function(dataset) {
  require(dplyr)
  require(magrittr)
  require(MASS)
  
  dataset %<>%
    mutate(
      Cond.treat = factor(Condition2, levels = c("Talker-specific", "Multi-talker", "Single talker", "Control")),
      Cond.treatTS.vs.CNTL = case_when(
        Cond.treat == "Talker-specific" ~ 1,
        Cond.treat == "Multi-talker" ~ 0,
        Cond.treat == "Single talker" ~ 0,
        T ~ 0
      ),
      Cond.treatMT.vs.CNTL = case_when(
        Cond.treat == "Talker-specific" ~ 0,
        Cond.treat == "Multi-talker" ~ 1,
        Cond.treat == "Single talker" ~ 0,
        T ~ 0
      ),
      Cond.treatST.vs.CNTL = case_when(
        Cond.treat == "Talker-specific" ~ 0,
        Cond.treat == "Multi-talker" ~ 0,
        Cond.treat == "Single talker" ~ 1,
        T ~ 0
      ),
      Cond.diff = factor(Condition2, levels = c("Talker-specific", "Multi-talker", "Single talker", "Control")),
      Cond.diffTS.vs.MT = case_when(
        Cond.diff == "Talker-specific" ~ 3/4,
        Cond.diff == "Multi-talker" ~ -1/4,
        Cond.diff == "Single talker" ~ -1/4,
        T ~ -1/4
      ),
      Cond.diffMT.vs.ST = case_when(
        Cond.diff == "Talker-specific" ~ 1/2,
        Cond.diff == "Multi-talker" ~ 1/2,
        Cond.diff == "Single talker" ~ -1/2,
        T ~ -1/2
      ),
      Cond.diffST.vs.CNTL = case_when(
        Cond.diff == "Talker-specific" ~ 1/4,
        Cond.diff == "Multi-talker" ~ 1/4,
        Cond.diff == "Single talker" ~ 1/4,
        T ~ -3/4
      )
    )

  # set sliding treaterence coding
  contrasts(dataset$Cond.treat) = cbind("TS.vs.CNTL" = c(1, 0, 0, 0),
                                        "MT.vs.CNTL" = c(0, 1, 0, 0),
                                        "ST.vs.CNTL" = c(0, 0, 1, 0))
  # set sliding difference coding
  contrasts(dataset$Cond.diff) = cbind("TS.vs.MT" = c(3/4, -1/4, -1/4, -1/4),
                                       "MT.vs.ST" = c(1/2, 1/2, -1/2, -1/2),
                                       "ST.vs.CNTL" = c(1/4, 1/4, 1/4, -3/4))
  
  return(dataset) 
}

d.exp1a.test = d.exp1ab.test %>% filter(Experiment == "1a") %>% defineContrasts()
d.exp1b.test = d.exp1ab.test %>% filter(Experiment == "1b") %>% defineContrasts()

# # Test newly created numerical variables that capture the 4-way condition contrast
# contrasts(d.exp1a.test$Cond.diff)
# d.exp1a.test %>%
#   dplyr::select(Cond.diff, Cond.diffTS.vs.MT, Cond.diffMT.vs.ST, Cond.diffST.vs.CNTL) %>%
#   distinct() %>% 
#   arrange(Cond.diff)
# 
# contrasts(d.exp1b.test$Cond.diff)
# d.exp1b.test %>%
#   dplyr::select(c(Cond.diff, Cond.diffTS.vs.MT, Cond.diffMT.vs.ST, Cond.diffST.vs.CNTL)) %>%
#   distinct() %>% 
#   arrange(Cond.diff)

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Step 2: Run model on original data. But with numerically coded variables instead of contrast
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
m.exp1a.test.diff <- brm(IsCorrect ~ Cond.diffTS.vs.MT + Cond.diffMT.vs.ST + Cond.diffST.vs.CNTL + Intercept.ranef.1ab + 
                           (1 + Cond.diffTS.vs.MT + Cond.diffMT.vs.ST + Cond.diffST.vs.CNTL | Sentence) + (1 | WorkerID),  
                         data = d.exp1a.test, family = bernoulli, iter = 11000, warmup = 1000, chains = 8, cores = 8,
                         prior = my.priors, 
                         save_model = "../models/m.exp1a.test.diff.full.numerically_coded_conditions.stan", 
                         file = "../models/exp1a.Test.diff.numerical_coded_conditions"
)

# Let's look at this model
stancode(m.exp1a.test.diff)

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Step 3: Derive priors from original data.
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
model = m.exp1a.test.diff
p = posterior_samples(model)


p %<>%
  dplyr::select(-starts_with("r_")) %>%
  dplyr::select(-starts_with("lp_")) %>%
  summarise_all(
    .funs = c("mean", "mode", "sd", "se") 
  ) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "coef") %>%
  mutate(
    property = gsub("^.*_([a-z]+)$", "\\1", coef),
    coef = gsub("^(.*)_[a-z]+$", "\\1", coef)
  ) %>%
  spread(property, V1) %>%
  mutate(
    class = gsub("^([a-z]+)_.*$", "\\1", coef),
    coef = gsub("^[a-z]+_(.*)$", "\\1", coef),
    group = ifelse(gsub("^([A-Za-z]+)__.*$", "\\1", coef) %in% c("Sentence", "WorkerID"), gsub("^([A-Za-z]+)__.*$", "\\1", coef), ""),
    coef = gsub("^[A-Za-z]+__(.*)$", "\\1", coef),
    prior = case_when(
      class == "b" ~ paste0("normal(", mean, ",", sd, ")"),
      class == "sd" ~ paste0("normal(", mean, ",", sd, ")"),
      #      class == "sd" ~ paste0("inv_gamma(", shape, ",", scale, ")"), # I checked in Stan reference manual that parameters for inv_gamma are indeed first shape and then scale
    )
  ) %>% 
  filter(class != "cor") %>%                          # For now removing correlation priors
  dplyr::select(prior, class, coef, group) %>%
  mutate(
    class = ifelse(coef == "Intercept" & class == "b", "Intercept", class),
    coef = ifelse(class == "Intercept", "", coef)
  )

my.newpriors = c(
  set_prior(p[1,"prior"], class = p[1,"class"], coef = p[1,"coef"], group = p[1,"group"]),
  set_prior(p[2,"prior"], class = p[2,"class"], coef = p[2,"coef"], group = p[2,"group"]),
  set_prior(p[3,"prior"], class = p[3,"class"], coef = p[3,"coef"], group = p[3,"group"]),
  set_prior(p[4,"prior"], class = p[4,"class"], coef = p[4,"coef"], group = p[4,"group"]),
  set_prior(p[5,"prior"], class = p[5,"class"], coef = p[5,"coef"], group = p[5,"group"]),
  set_prior(p[6,"prior"], class = p[6,"class"], coef = p[6,"coef"], group = p[6,"group"]),
  set_prior(p[7,"prior"], class = p[7,"class"], coef = p[7,"coef"], group = p[7,"group"]),
  set_prior(p[8,"prior"], class = p[8,"class"], coef = p[8,"coef"], group = p[8,"group"]),
  set_prior(p[9,"prior"], class = p[9,"class"], coef = p[9,"coef"], group = p[9,"group"]),
  set_prior(p[10,"prior"], class = p[10,"class"], coef = p[10,"coef"], group = p[10,"group"]),
  prior(lkj(1), class = cor)
)
my.newpriors.woTS.vs.MT = my.newpriors[-3,]
my.newpriors.woMT.vs.ST = my.newpriors[-1,]
my.newpriors.woST.vs.CNTL = my.newpriors[-2,]


## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Step 4: Run same numerically-coded model on new data, but with priors based on original data  
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
m.exp1b.test.diff.full <- brm(IsCorrect ~ Cond.diffTS.vs.MT + Cond.diffMT.vs.ST + Cond.diffST.vs.CNTL + Intercept.ranef.1ab + 
                                (1 + Cond.diffTS.vs.MT + Cond.diffMT.vs.ST + Cond.diffST.vs.CNTL | Sentence) + (1 | WorkerID),  
                              data = d.exp1b.test, family = bernoulli, iter = 11000, warmup = 1000, chains = 8, cores = 8,
                              prior = my.newpriors, save_all_pars = T, 
                              save_model = "../models/m.exp1b.test.diff.full.numerically_coded_conditions.stan", 
                              file = "../models/exp1b.Test.diff.numerical_coded_conditions"
)


## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Step 4: Run all three 'nullmodels' -- i.e. the same numerically-coded model on new data, but 
#         with one parameter set to 0 (removed). Again use priors based on original data. 
#         
#         NB: We're keeping the random effects for the null effect in the model.
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
m.exp1b.test.diff.woTS.vs.MT <- brm(IsCorrect ~ Cond.diffMT.vs.ST + Cond.diffST.vs.CNTL + Intercept.ranef.1ab + 
                                      (1 + Cond.diffTS.vs.MT + Cond.diffMT.vs.ST + Cond.diffST.vs.CNTL | Sentence) + (1 | WorkerID),  
                                    data = d.exp1b.test, family = bernoulli, iter = 11000, warmup = 1000, chains = 8, cores = 8,
                                    prior = my.newpriors.woTS.vs.MT, save_all_pars = T, 
                                    save_model = "../models/m.exp1b.test.diff.woTS.vs.MT.numerically_coded_conditions.stan",
                                    file = "../models/exp1b.Test.diff.woTS.vs.MT.numerical_coded_conditions"
)
m.exp1b.test.diff.woMT.vs.ST <- brm(IsCorrect ~ Cond.diffTS.vs.MT + Cond.diffST.vs.CNTL + Intercept.ranef.1ab + 
                                      (1 + Cond.diffTS.vs.MT + Cond.diffMT.vs.ST + Cond.diffST.vs.CNTL | Sentence) + (1 | WorkerID),  
                                    data = d.exp1b.test, family = bernoulli, iter = 11000, warmup = 1000, chains = 8, cores = 8,
                                    prior = my.newpriors.woMT.vs.ST, save_all_pars = T, 
                                    save_model = "../models/m.exp1b.test.diff.woMT.vs.ST.numerically_coded_conditions.stan",
                                    file = "../models/exp1b.Test.diff.woMT.vs.ST.numerical_coded_conditions"
)
m.exp1b.test.diff.woST.vs.CNTL <- brm(IsCorrect ~ Cond.diffTS.vs.MT + Cond.diffMT.vs.ST + Intercept.ranef.1ab + 
                                        (1 + Cond.diffTS.vs.MT + Cond.diffMT.vs.ST + Cond.diffST.vs.CNTL | Sentence) + (1 | WorkerID),  
                                      data = d.exp1b.test, family = bernoulli, iter = 11000, warmup = 1000, chains = 8, cores = 8,
                                      prior = my.newpriors.woST.vs.CNTL, save_all_pars = T, 
                                      save_model = "../models/m.exp1b.test.diff.woST.vs.CNTL.numerically_coded_conditions.stan",
                                      file = "../models/exp1b.Test.diff.woST.vs.CNTL.numerical_coded_conditions"
)

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Step 5: Run the bridge sampler on all new models fit to the replication data
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
b1 = bridge_sampler(samples = m.exp1b.test.diff.full, cores = 8, method = "warp3", maxiter = 1000)
saveRDS(b1, file = "../models/bridge_exp1b.Test.diff.numerical_coded_conditions.rds")
b2.woTS.vs.MT = bridge_sampler(samples = m.exp1b.test.diff.woTS.vs.MT, cores = 8, method = "warp3", maxiter = 1000)
saveRDS(b2.woTS.vs.MT, file = "../models/bridge_exp1b.Test.diff.woTS.vs.MT.numerical_coded_conditions.rds")
b2.woMT.vs.ST = bridge_sampler(samples = m.exp1b.test.diff.woMT.vs.ST, cores = 8, method = "warp3", maxiter = 1000)
saveRDS(b2.woMT.vs.ST, file = "../models/bridge_exp1b.Test.diff.woMT.vs.ST.numerical_coded_conditions.rds")
b2.woST.vs.CNTL = bridge_sampler(samples = m.exp1b.test.diff.woST.vs.CNTL, cores = 8, method = "warp3", maxiter = 1000)
saveRDS(b2.woST.vs.CNTL, file = "../models/bridge_exp1b.Test.diff.woST.vs.CNTL.numerical_coded_conditions.rds")
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Step 6: Calculate all the relevant Bayes Factors
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
b1 = readRDS(file = "../models/bridge_exp1b.Test.diff.numerical_coded_conditions.rds")
b2.woTS.vs.MT = readRDS(file = "../models/bridge_exp1b.Test.diff.woTS.vs.MT.numerical_coded_conditions.rds")
b2.woMT.vs.ST = readRDS(file = "../models/bridge_exp1b.Test.diff.woMT.vs.ST.numerical_coded_conditions.rds")
b2.woST.vs.CNTL = readRDS(file = "../models/bridge_exp1b.Test.diff.woST.vs.CNTL.numerical_coded_conditions.rds")

bf(b1, b2.woTS.vs.MT)
bf(b1, b2.woMT.vs.ST)
bf(b1, b2.woST.vs.CNTL)


# Use the previous data to establish a prior. You may want to construct a normal prior centered at their estimate, with a standard deviation equal to their SE estimate. Do this for every coefficient.
# Construct two models, WITHOUT using the ~ syntax. You need to use the target += normal_lpdf(…) syntax. (This is because we will need to compute the marginal likelihood, with all normalizing constants present, and the ~ syntax drops those constants). The first model has the full model, with priors defined in step 1. The second model has the reduced model, with the coefficient of interest implicitly set to zero (as in, it isn’t present). Everything else is the same as in the first model.
# Estimate both models using the replication data.


## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# AND NOW FOR TREATMENT-CODED MODELS
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Contrasts were already set above

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Step 2: Run model on original data. But with numerically coded variables instead of contrast
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
m.exp1a.test.treat <- brm(IsCorrect ~ Cond.treatTS.vs.CNTL + Cond.treatMT.vs.CNTL + Cond.treatST.vs.CNTL + Intercept.ranef.1ab + 
                           (1 + Cond.treatTS.vs.CNTL + Cond.treatMT.vs.CNTL + Cond.treatST.vs.CNTL | Sentence) + (1 | WorkerID),  
                         data = d.exp1a.test, family = bernoulli, iter = 11000, warmup = 1000, chains = 8, cores = 8,
                         prior = my.priors, 
                         save_model = "../models/m.exp1a.test.treat.full.numerically_coded_conditions.stan", 
                         file = "../models/exp1a.Test.treat.numerical_coded_conditions"
)

# Let's look at this model
stancode(m.exp1a.test.treat)

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Step 3: Derive priors from original data.
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
model = m.exp1a.test.treat
p = posterior_samples(model)


p %<>%
  dplyr::select(-starts_with("r_")) %>%
  dplyr::select(-starts_with("lp_")) %>%
  summarise_all(
    .funs = c("mean", "mode", "sd", "se") 
  ) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "coef") %>%
  mutate(
    property = gsub("^.*_([a-z]+)$", "\\1", coef),
    coef = gsub("^(.*)_[a-z]+$", "\\1", coef)
  ) %>%
  spread(property, V1) %>%
  mutate(
    class = gsub("^([a-z]+)_.*$", "\\1", coef),
    coef = gsub("^[a-z]+_(.*)$", "\\1", coef),
    group = ifelse(gsub("^([A-Za-z]+)__.*$", "\\1", coef) %in% c("Sentence", "WorkerID"), gsub("^([A-Za-z]+)__.*$", "\\1", coef), ""),
    coef = gsub("^[A-Za-z]+__(.*)$", "\\1", coef),
    prior = case_when(
      class == "b" ~ paste0("normal(", mean, ",", sd, ")"),
      class == "sd" ~ paste0("normal(", mean, ",", sd, ")"),
      #      class == "sd" ~ paste0("inv_gamma(", shape, ",", scale, ")"), # I checked in Stan reference manual that parameters for inv_gamma are indeed first shape and then scale
    )
  ) %>% 
  filter(class != "cor") %>%                          # For now removing correlation priors
  dplyr::select(prior, class, coef, group) %>%
  mutate(
    class = ifelse(coef == "Intercept" & class == "b", "Intercept", class),
    coef = ifelse(class == "Intercept", "", coef)
  )

my.newpriors = c(
  set_prior(p[1,"prior"], class = p[1,"class"], coef = p[1,"coef"], group = p[1,"group"]),
  set_prior(p[2,"prior"], class = p[2,"class"], coef = p[2,"coef"], group = p[2,"group"]),
  set_prior(p[3,"prior"], class = p[3,"class"], coef = p[3,"coef"], group = p[3,"group"]),
  set_prior(p[4,"prior"], class = p[4,"class"], coef = p[4,"coef"], group = p[4,"group"]),
  set_prior(p[5,"prior"], class = p[5,"class"], coef = p[5,"coef"], group = p[5,"group"]),
  set_prior(p[6,"prior"], class = p[6,"class"], coef = p[6,"coef"], group = p[6,"group"]),
  set_prior(p[7,"prior"], class = p[7,"class"], coef = p[7,"coef"], group = p[7,"group"]),
  set_prior(p[8,"prior"], class = p[8,"class"], coef = p[8,"coef"], group = p[8,"group"]),
  set_prior(p[9,"prior"], class = p[9,"class"], coef = p[9,"coef"], group = p[9,"group"]),
  set_prior(p[10,"prior"], class = p[10,"class"], coef = p[10,"coef"], group = p[10,"group"]),
  prior(lkj(1), class = cor)
)
my.newpriors.woTS.vs.CNTL = my.newpriors[-3,]
my.newpriors.woMT.vs.CNTL = my.newpriors[-1,]
my.newpriors.woST.vs.CNTL = my.newpriors[-2,]


## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Step 4: Run same numerically-coded model on new data, but with priors based on original data  
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
m.exp1b.test.treat.full <- brm(IsCorrect ~ Cond.treatTS.vs.CNTL + Cond.treatMT.vs.CNTL + Cond.treatST.vs.CNTL + Intercept.ranef.1ab + 
                                (1 + Cond.treatTS.vs.CNTL + Cond.treatMT.vs.CNTL + Cond.treatST.vs.CNTL | Sentence) + (1 | WorkerID),  
                              data = d.exp1b.test, family = bernoulli, iter = 11000, warmup = 1000, chains = 8, cores = 8,
                              prior = my.newpriors, save_all_pars = T, 
                              save_model = "../models/m.exp1b.test.treat.full.numerically_coded_conditions.stan", 
                              file = "../models/exp1b.Test.treat.numerical_coded_conditions"
)


## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Step 4: Run all three 'nullmodels' -- i.e. the same numerically-coded model on new data, but 
#         with one parameter set to 0 (removed). Again use priors based on original data. 
#         
#         NB: We're keeping the random effects for the null effect in the model.
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
m.exp1b.test.treat.woTS.vs.CNTL <- brm(IsCorrect ~ Cond.treatMT.vs.CNTL + Cond.treatST.vs.CNTL + Intercept.ranef.1ab + 
                                      (1 + Cond.treatTS.vs.CNTL + Cond.treatMT.vs.CNTL + Cond.treatST.vs.CNTL | Sentence) + (1 | WorkerID),  
                                    data = d.exp1b.test, family = bernoulli, iter = 11000, warmup = 1000, chains = 8, cores = 8,
                                    prior = my.newpriors.woTS.vs.CNTL, save_all_pars = T, 
                                    save_model = "../models/m.exp1b.test.treat.woTS.vs.CNTL.numerically_coded_conditions.stan",
                                    file = "../models/exp1b.Test.treat.woTS.vs.CNTL.numerical_coded_conditions"
)
m.exp1b.test.treat.woMT.vs.CNTL <- brm(IsCorrect ~ Cond.treatTS.vs.CNTL + Cond.treatST.vs.CNTL + Intercept.ranef.1ab + 
                                      (1 + Cond.treatTS.vs.CNTL + Cond.treatMT.vs.CNTL + Cond.treatST.vs.CNTL | Sentence) + (1 | WorkerID),  
                                    data = d.exp1b.test, family = bernoulli, iter = 11000, warmup = 1000, chains = 8, cores = 8,
                                    prior = my.newpriors.woMT.vs.CNTL, save_all_pars = T, 
                                    save_model = "../models/m.exp1b.test.treat.woMT.vs.CNTL.numerically_coded_conditions.stan",
                                    file = "../models/exp1b.Test.treat.woMT.vs.CNTL.numerical_coded_conditions"
)
m.exp1b.test.treat.woST.vs.CNTL <- brm(IsCorrect ~ Cond.treatTS.vs.CNTL + Cond.treatMT.vs.CNTL + Intercept.ranef.1ab + 
                                        (1 + Cond.treatTS.vs.CNTL + Cond.treatMT.vs.CNTL + Cond.treatST.vs.CNTL | Sentence) + (1 | WorkerID),  
                                      data = d.exp1b.test, family = bernoulli, iter = 11000, warmup = 1000, chains = 8, cores = 8,
                                      prior = my.newpriors.woST.vs.CNTL, save_all_pars = T, 
                                      save_model = "../models/m.exp1b.test.treat.woST.vs.CNTL.numerically_coded_conditions.stan",
                                      file = "../models/exp1b.Test.treat.woST.vs.CNTL.numerical_coded_conditions"
)

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Step 5: Run the bridge sampler on all new models fit to the replication data
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
b1 = bridge_sampler(samples = m.exp1b.test.treat.full, cores = 8, method = "warp3", maxiter = 1000)
saveRDS(b1, file = "../models/bridge_exp1b.Test.treat.numerical_coded_conditions.rds")
b2.woTS.vs.CNTL = bridge_sampler(samples = m.exp1b.test.treat.woTS.vs.CNTL, cores = 8, method = "warp3", maxiter = 1000)
saveRDS(b2.woTS.vs.CNTL, file = "../models/bridge_exp1b.Test.treat.woTS.vs.CNTL.numerical_coded_conditions.rds")
b2.woMT.vs.CNTL = bridge_sampler(samples = m.exp1b.test.treat.woMT.vs.CNTL, cores = 8, method = "warp3", maxiter = 1000)
saveRDS(b2.woMT.vs.CNTL, file = "../models/bridge_exp1b.Test.treat.woMT.vs.CNTL.numerical_coded_conditions.rds")
b2.woST.vs.CNTL = bridge_sampler(samples = m.exp1b.test.treat.woST.vs.CNTL, cores = 8, method = "warp3", maxiter = 1000)
saveRDS(b2.woST.vs.CNTL, file = "../models/bridge_exp1b.Test.treat.woST.vs.CNTL.numerical_coded_conditions.rds")
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Step 6: Calculate all the relevant Bayes Factors
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bf(b1, b2.woTS.vs.CNTL)
bf(b1, b2.woMT.vs.CNTL)
bf(b1, b2.woST.vs.CNTL)










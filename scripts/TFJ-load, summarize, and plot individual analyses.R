# This script provides examples of how to load models saved my the R markdown script, and how to 
# summarize or visualize those models. This script was NOT used to generate the analyses reported 
# in the paper, but is rather meant to help other researchers interested in interacting with the 
# our data. 

library(dplyr)
library(tibble)
library(magrittr)
library(tidyr)
library(brms)
library(shinystan)
library(knitr)
source("../scripts/functions.R")

# load data
d.exp1ab.training = readRDS("../data/Exp1ab.training postRMD.RDS")
d.exp1ab.test = readRDS("../data/Exp1ab.test postRMD.RDS")

# load models
setwd("../models/")

# Experiment 1a
m.exp1a.test.treat = readRDS("exp1a.Test.treat.RDS")
m.exp1a.test.diff = readRDS("exp1a.Test.diff.RDS")
m.exp1a.test.treat.td = readRDS("exp1a.Test.treat.Talker-dependent.RDS")
m.exp1a.test.diff.td = readRDS("exp1a.Test.diff.Talker-dependent.RDS")

# Experiment 1b
m.exp1b.test.treat = readRDS("exp1b.Test.treat.RDS")
m.exp1b.test.diff = readRDS("exp1b.Test.diff.RDS")
m.exp1b.test.treat.td = readRDS("exp1b.Test.treat.Talker-dependent.RDS")
m.exp1b.test.diff.td = readRDS("exp1b.Test.diff.Talker-dependent.RDS")

launch_shinystan(m.exp1a.test.diff)

## -------------------------------------------------------------------------------
# By condition means (for comparison to BB08)
## -------------------------------------------------------------------------------
rau = function(s, n = 16) {
  AU = asin(sqrt(s / (n+1))) + asin(sqrt((s+1) / (n+1)))
  RAU = 146/pi *AU - 23
  
  return(RAU)
}

d.exp1ab.test %>%
  group_by(Experiment, Condition2, WorkerID) %>%
  dplyr::summarise(
    meanAccuracy = mean(IsCorrect),
    meanRAU = rau(sum(IsCorrect), length(IsCorrect))
  ) %>%
  group_by(Experiment, Condition2) %>%
  dplyr::summarise(
    meanAccuracy = mean(meanAccuracy),
    meanRAU = mean(meanRAU)
  ) %>%
  mutate(
    meanRAUacc = rau(meanAccuracy * 80, 80)
  )

## -------------------------------------------------------------------------------
# Example table of fit
## -------------------------------------------------------------------------------
# Summary of Bayesian mixed-logistic regression with sliding difference coded exposure conditions
s = summary(m.exp1a.test.diff)$fixed
colnames(s) = c("Est.", "SE", "$CI_{lower}$", "$CI_{upper}$", "Eff. samples", "$\\hat{R}$")
rownames(s) = c("Intercept", "TS vs. MT", "MT vs. ST", "ST vs. CNTL", "{\\em a priori} performance")
kable(s, format = "latex", digits = c(2, 3, 2, 2, 1, 2), 
      caption = "Experiment 1a")

s = summary(m.exp1b.test.diff)$fixed
colnames(s) = c("Est.", "SE", "$CI_{lower}$", "$CI_{upper}$", "Eff. samples", "$\\hat{R}$")
rownames(s) = c("Intercept", "TS vs. MT", "MT vs. ST", "ST vs. CNTL", "{\\em a priori} performance")
kable(s, format = "latex", digits = c(2, 3, 2, 2, 1, 2), 
      caption = "Experiment 1b")



xhypotheses(
  list(
    hypothesis(m.exp1a.test.diff, "Cond.diffTS.vs.MT + Cond.diffMT.vs.ST + Cond.diffST.vs.CNTL > 0", class = "b"),
    hypothesis(m.exp1a.test.diff, "Cond.diffMT.vs.ST + Cond.diffST.vs.CNTL > 0", class = "b"),
    hypothesis(m.exp1a.test.diff, "Cond.diffST.vs.CNTL > 0", class = "b"),
    hypothesis(m.exp1a.test.diff, "Cond.diffMT.vs.ST > 0", class = "b")
  ),
  labels = c("Adaptation: TS vs. CNTL", "Question 1: MT vs. CNTL", "Question 2: ST vs. CNTL", "Question 3: MT vs. ST"), 
  format = "latex"
)

xhypotheses(
  list(
    hypothesis(m.exp1b.test.diff, "Cond.diffTS.vs.MT + Cond.diffMT.vs.ST + Cond.diffST.vs.CNTL > 0", class = "b"),
    hypothesis(m.exp1b.test.diff, "Cond.diffMT.vs.ST + Cond.diffST.vs.CNTL > 0", class = "b"),
    hypothesis(m.exp1b.test.diff, "Cond.diffST.vs.CNTL > 0", class = "b"),
    hypothesis(m.exp1b.test.diff, "Cond.diffMT.vs.ST > 0", class = "b")
  ),
  labels = c("Adaptation: TS vs. CNTL", "Question 1: MT vs. CNTL", "Question 2: ST vs. CNTL", "Question 3: MT vs. ST"), 
  format = "latex"
)

## -------------------------------------------------------------------------------
# Example visualization of fit
## -------------------------------------------------------------------------------

library(tidybayes)
d.exp1ab.test %>%
  filter(Experiment == "1a") %>%
  data_grid(Cond.diff, .model = m.exp1a.test.diff) %>%
  # by default all random effects are considered in prediction
  # choose scale = "linear" for log-odds prediction, "response" for proportions
  add_fitted_draws(m.exp1a.test.diff, scale = "linear") %>% 
  ggplot(aes(x = .value, y = Cond.diff)) +
  geom_halfeyeh(.width = c(.5, .95)) + 
  scale_x_continuous("Predicted log-odds correct") +
  scale_y_discrete("Exposure condition") 

m.exp1a.Test.diff %>%
  gather_draws(b_Cond.diffTS.vs.MT, b_Cond.diffMT.vs.ST, b_Cond.diffST.vs.CNTL) %>%
  ggplot(aes(y = .variable, x = .value)) +
  geom_halfeyeh(color = "blue") + 
  scale_x_continuous(expression(paste(hat(beta)," (in log-odds)"))) +
  scale_y_discrete("",
                   labels = c("Single talker vs.\nControl",
                              "Multi-talker vs.\nSingle talker", 
                              "Talker-specific vs.\nMulti-talker", 
                              "Intercept"
                   )
  ) + 
  geom_vline(xintercept = 0, linetype = 2) 
  
myGplot.defaults(type = "paper")
m.exp1a.test.treat %>% 
  gather_draws(b_Cond.treatTS.vs.CNTL, b_Cond.treatMT.vs.CNTL) %>%
  rbind(
    m.exp1a.test.diff %>%
      gather_draws(b_Cond.diffMT.vs.ST, b_Cond.diffST.vs.CNTL)
  ) %>%
  ggplot(aes(y = .variable, x = .value)) +
  #  geom_rect(xmin = -1, xmax = 2, ymin = .8, ymax = 1.8, fill = "darkgray", alpha = .3, inherit.aes = FALSE) +
  geom_halfeyeh(color = "blue", fill = "darkgray", relative_scale = .8) + 
  scale_x_continuous(
    expression(paste(hat(beta)," (in log-odds)")),
    limits = c(-.5, 1.4)
  ) +
  scale_y_discrete("",
                   labels = c(
                     expression(paste(bold("Question 3   "), "Multi-talker vs.\nSingle talker")),
                     expression(paste(bold("Question 2   "), "Single talker vs.\nControl")),
                     expression(paste(bold("Question 1   "), "Multi-talker vs.\nControl")),
                     expression(paste(bold("Adaptation   "), "Talker-specific vs.\nControl"))
                   )
  ) +
  coord_cartesian(ylim = c(.99, 4.5)) +
  geom_vline(xintercept = 0, linetype = 2) +
  ggtitle("Experiment 1a") + 
  theme(axis.text.y = element_text(hjust=0), panel.grid.major.y = element_blank()) 







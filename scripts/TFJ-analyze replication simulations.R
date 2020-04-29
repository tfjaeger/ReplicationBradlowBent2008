library(tidyverse)
source("../scripts/functions.R")
simulation = readRDS("../models/replication_simulations_orig and rep-MTvsCNTL_same intercept.rds")

length(which(simulation$BF > 20)) / nrow(simulation)
length(which(simulation$BF > 1)) / nrow(simulation)

library(ggplot2)

# density of BF values
myGplot.defaults(type = "paper")
xlimits = c(10^-2, 10^7)
p = ggplot(simulation,
           aes(
             x = BF)
) +
  geom_density(n = 2 ^ 10) + theme_bw() +
  scale_x_log10(
    expression(paste("BF"["H"["rep"]], ""["H"["0"]])),
    limits = xlimits
  ) + coord_cartesian(xlim = xlimits, ylim = c(0,1))

d <- ggplot_build(p)$data[[1]]
p + geom_area(data = subset(d, x >= 0), aes(x=10^x, y=y), fill="gray") +
  geom_vline(xintercept = 3, linetype = 2, color = "black") +
  geom_vline(xintercept = 20, linetype = 2, color = "black") +
  geom_vline(xintercept = 150, linetype = 2, color = "black") +
  geom_text(label = paste0(round(length(which(simulation$BF >= 3)) / length(simulation$BF) * 100, 1),"% positive"), 
            angle = 90, x = log10(6), y = 0.005, color = "black", hjust = 0, vjust = 0) + 
  geom_text(label = paste0(round(length(which(simulation$BF >= 20)) / length(simulation$BF) * 100, 1),"% strong"),
            angle = 90, x = log10(40), y = 0.005, color = "black", hjust = 0, vjust = 0) +
  geom_text(label = paste0(round(length(which(simulation$BF >= 150)) / length(simulation$BF) * 100, 1),"% very strong"), 
            angle = 90, x = log10(300), y = 0.005, color = "black", hjust = 0, vjust = 0) +
  geom_text(label = paste0("BFs supporting replication: ", round(length(which(simulation$BF > 1)) / length(simulation$BF) * 100, 1), "%"),
            angle = 0, x = log10(min(xlimits)), y = 1, hjust = 0, color = "darkgray") 
#ggsave("../figures/Power of simulation for ST effect.pdf", width = 6, height = 3)
#ggsave("../figures/Power of simulation for MT effect.pdf", width = 6, height = 3)
ggsave("../figures/Power of simulation for MT effect_same intercept.pdf", width = 6, height = 3)


simulation = 
  rbind(
    readRDS("../models/replication_simulations_orig and rep-STvsCNTL.RDS") %>%
      mutate(Comparison = "ST vs. CNTL"),
    readRDS("../models/replication_simulations_orig and rep-MTvsCNTL.RDS") %>%
      mutate(Comparison = "MT vs. CNTL")
  )

ggplot(simulation,
       aes(
         x = Int.2.estimate,
         y = BF,
         color = Comparison)
) +
  geom_point() +
  geom_hline(yintercept = 1, linetype = 1, color = "gray") +
  geom_hline(yintercept = 3, linetype = 2, color = "gray") +
  geom_hline(yintercept = 20, linetype = 2, color = "gray") +
  geom_hline(yintercept = 150, linetype = 2, color = "gray") +
  geom_smooth() +
  scale_x_continuous("Intercept estimate of replication\n(in log-odds)") +
  scale_y_log10("Replication BF") + 
  scale_color_manual(values = setNames(c("#7CAE00", "#C77CFF"), c("ST vs. CNTL", "MT vs. CNTL")))
ggsave("../figures/Power by intercept of replication (ST and MT).pdf", width = 6, height = 3)

ggplot(simulation,
       aes(
         x = Int.2.estimate - Int.1.estimate,
         y = Cond.2.estimate,
         z = log10(BF),
         size = ((Cond.1.estimate - Cond.1.lower) * 2),
         color = log10(BF))
) +
  geom_density2d(color = "gray") +
  geom_point() +
  scale_color_gradient2(high = "green", mid = "blue", low="red")  + theme_bw() + 
  scale_x_continuous("Replication intercept - original intercept") +
  scale_y_continuous("Effect size estimate (replication)") +
  scale_size_continuous("Width of 95% CI\n(replication)", trans = "reverse")
ggsave("../figures/Power of simulation for ST and MT effect by intercept difference between replication and original.pdf")

# # Some sanity checks
# my.data.rep = simple.data.generator(
#   estimates = cbind("Intercept" = 2.00, "TS" = .49),
#   n.subj.perCondition = nsubj,
#   n.obs = nitem
# ) 
# 
# my.data.rep() %>%
#   filter(Cond == "Control") %>%
#   summarise(meanOutcome = qlogis(mean(Outcome)), meanY = mean(y))


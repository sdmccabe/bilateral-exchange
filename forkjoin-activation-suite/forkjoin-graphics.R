#forkjoin analysis
library(readr)
library(dplyr)
library(ggplot2)
library(scales)

#d <- read_csv("~/Dropbox/data-pull/forkjoin-run-3-3.csv")
d <- read_csv("~/c/exchange/results/new/forkjoin_results.csv")
# Drop some bad data
d$activation.method <- factor(d$activation.method, levels = c("-1", "0", "1", "2", "3", "4", "5"), labels = c("Fixed", "Random", "Uniform", "Poisson-1", "Poisson-Poor", "Poisson-Rich", "Poisson-Middle"))
d <- dplyr::filter(d, activation.method != "Poisson-1") # re-do this
d <- dplyr::mutate(d, run = (row_number() - 1) %/% 100)


interactions.plot <- ggplot(d, aes(x=time, y=interactions, group = run, color = factor(activation.method))) + 
  #geom_point(alpha = 0.01) + 
  stat_summary(fun.y = "mean", geom = "line", size = 1, alpha = 1, aes(group = factor(activation.method))) + 
  scale_color_discrete(name="Activation Method") +
  theme_classic(base_family = "serif", base_size = 12) +
  xlab("Time") +
  ylab("Number of Interactions") +
  scale_y_continuous(labels = comma) +
  ggtitle("Number of Interactions over Time")

wealthmin.plot <- ggplot(d, aes(x=time, y=currentwealth.min, group = run, color = factor(activation.method))) + 
  #geom_point(alpha = 0.01) + 
  stat_summary(fun.y = "mean", geom = "line", size = 1, alpha = 1, aes(group = factor(activation.method))) + 
  scale_color_discrete(name="Activation Method") +
  theme_classic(base_family = "serif", base_size = 12) +
  xlab("Time") +
  ylab("Minimum Agent Wealth") +
  scale_y_continuous(labels = comma) +
  ggtitle("Minimum Agent Wealth over Time")

wealthmax.plot <- ggplot(d, aes(x=time, y=currentwealth.max, group = run, color = factor(activation.method))) + 
  #geom_point(alpha = 0.01) + 
  stat_summary(fun.y = "mean", geom = "line", size = 1, alpha = 1, aes(group = factor(activation.method))) + 
  scale_color_discrete(name="Activation Method") +
  theme_classic(base_family = "serif", base_size = 12) +
  xlab("Time") +
  ylab("Maximum Agent Wealth") +
  scale_y_continuous(labels = comma)+ 
  ggtitle("Max Agent Wealth over Time")

wealthavg.plot <- ggplot(d, aes(x=time, y=currentwealth.avg, group = run, color = factor(activation.method))) + 
  #geom_point(alpha = 0.01) + 
  stat_summary(fun.y = "mean", geom = "line", size = 1, alpha = 1, aes(group = factor(activation.method))) + 
  scale_color_discrete(name="Activation Method") +
  theme_classic(base_family = "serif", base_size = 12) +
  xlab("Time") +
  ylab("Average Agent Wealth") +
  scale_y_continuous(labels = comma) +
  ggtitle("Mean Agent Wealth over Time")

wealthsd.plot <- ggplot(d, aes(x=time, y=currentwealth.sd, group = run, color = factor(activation.method))) + 
  #geom_point(alpha = 0.01) + 
  stat_summary(fun.y = "mean", geom = "line", size = 1, alpha = 1, aes(group = factor(activation.method))) + 
  scale_color_discrete(name="Activation Method") +
  theme_classic(base_family = "serif", base_size = 12) +
  xlab("Time") +
  ylab("S. D. of Agent Wealth") +
  scale_y_continuous(labels = comma) +
  ggtitle("Standard Deviation of Agent Wealth over Time")

utilitymin.plot <- ggplot(d, aes(x=time, y=utility.min, group = run, color = factor(activation.method))) + 
  #geom_point(alpha = 0.01) + 
  stat_summary(fun.y = "mean", geom = "line", size = 1, alpha = 1, aes(group = factor(activation.method))) + 
  scale_color_discrete(name="Activation Method") +
  theme_classic(base_family = "serif", base_size = 12) +
  xlab("Time") +
  ylab("Minimum Agent Utility") +
  scale_y_continuous(labels = comma) +
  ggtitle("Minimum Agent Utility over Time")

utilitymax.plot <- ggplot(d, aes(x=time, y=utility.max, group = run, color = factor(activation.method))) + 
  #geom_point(alpha = 0.01) + 
  stat_summary(fun.y = "mean", geom = "line", size = 1, alpha = 1, aes(group = factor(activation.method))) + 
  scale_color_discrete(name="Activation Method") +
  theme_classic(base_family = "serif", base_size = 12) +
  xlab("Time") +
  ylab("Maximum Agent Utility") +
  scale_y_continuous(labels = comma) + 
  ggtitle("Maximum Agent Utility over Time")

utilityavg.plot <- ggplot(d, aes(x=time, y=utility.avg, group = run, color = factor(activation.method))) + 
  #geom_point(alpha = 0.01) + 
  stat_summary(fun.y = "mean", geom = "line", size = 1, alpha = 1, aes(group = factor(activation.method))) + 
  scale_color_discrete(name="Activation Method") +
  theme_classic(base_family = "serif", base_size = 12) +
  xlab("Time") +
  ylab("Average Agent Utility") +
  scale_y_continuous(labels = comma) + 
  ggtitle("Average Agent Utility over Time")

utilitysd.plot <- ggplot(d, aes(x=time, y=utility.sd, group = run, color = factor(activation.method))) + 
  #geom_point(alpha = 0.01) + 
  stat_summary(fun.y = "mean", geom = "line", size = 1, alpha = 1, aes(group = factor(activation.method))) + 
  scale_color_discrete(name="Activation Method") +
  theme_classic(base_family = "serif", base_size = 12) +
  xlab("Time") +
  ylab("S. D. of Agent Utility") +
  scale_y_continuous(labels = comma) + 
  ggtitle("Standard Deviation of Agent Utility Over Time")

l2sdmrs.plot <- ggplot(d, aes(x=time, y=L2.sd.MRS, group = run, color = factor(activation.method))) + 
  #geom_point(alpha = 0.01) + 
  stat_summary(fun.y = "mean", geom = "line", size = 1, alpha = 1, aes(group = factor(activation.method))) + 
  scale_color_discrete(name="Activation Method") +
  theme_classic(base_family = "serif", base_size = 12) +
  xlab("Time") +
  ylab("L2 S. D. of Agent MRSes") +
  scale_y_continuous(labels = comma) + 
  ggtitle("Standard Deviation of Agent MRSes over Time")

maxsdmrs.plot <- ggplot(d, aes(x=time, y=max.sd.MRS, group = run, color = factor(activation.method))) + 
  #geom_point(alpha = 0.01) + 
  stat_summary(fun.y = "mean", geom = "line", size = 1, alpha = 1, aes(group = factor(activation.method))) + 
  scale_color_discrete(name="Activation Method") +
  theme_classic(base_family = "serif", base_size = 12) +
  xlab("Time") +
  ylab("Max S. D. of Agent MRSes") +
  scale_y_continuous(labels = comma) +
  ggtitle("Standard Deviation of Agent MRSes over Time")


setwd("~/Dropbox/Thesis/images/")
ggsave("forkjoin-interactions.pdf", interactions.plot, width = 6.2, height = 3, units = "in")
ggsave("forkjoin-wealth-max.pdf", wealthmax.plot, width = 6.2, height = 3, units = "in")
ggsave("forkjoin-wealth-min.pdf", wealthmin.plot, width = 6.2, height = 3, units = "in")
ggsave("forkjoin-wealth-sd.pdf", wealthsd.plot, width = 6.2, height = 3, units = "in")
ggsave("forkjoin-wealth-avg.pdf", wealthavg.plot, width = 6.2, height = 3, units = "in")
ggsave("forkjoin-utility-max.pdf", utilitymax.plot, width = 6.2, height = 3, units = "in")
ggsave("forkjoin-utlity-min.pdf", utilitymin.plot, width = 6.2, height = 3, units = "in")
ggsave("forkjoin-utility-sd.pdf", utilitysd.plot, width = 6.2, height = 3, units = "in")
ggsave("forkjoin-utility-avg.pdf", utilityavg.plot, width = 6.2, height = 3, units = "in")
ggsave("forkjoin-l2-sd-mrs.pdf", l2sdmrs.plot, width = 6.2, height = 3, units = "in")
ggsave("forkjoin-max-sd-mrs.pdf", maxsdmrs.plot, width = 6.2, height = 3, units = "in")

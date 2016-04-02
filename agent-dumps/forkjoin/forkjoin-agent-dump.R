library(readr)
library(dplyr)
library(ggplot2)
library(ggthemes)

fixed<-read_csv("~/c/exchange/agent-dumps/forkjoin/fixedactivations.csv",col_names=F)
fixed <- dplyr::rename(fixed, id=X1, activations=X2, trades=X3)
fixed$type = "Fixed"
random<-read_csv("~/c/exchange/agent-dumps/forkjoin/randomactivations.csv",col_names=F)
random <- dplyr::rename(random, id=X1, activations=X2, trades=X3)
random$type = "Random"
uniform<-read_csv("~/c/exchange/agent-dumps/forkjoin/uniformactivations.csv",col_names=F)
uniform <- dplyr::rename(uniform, id=X1, activations=X2, trades=X3)
uniform$type = "Uniform"
# poisson1<-read_csv("~/c/exchange/agent-dumps/poisson1activations.csv",col_names=F)
# poisson1 <- dplyr::rename(poisson1, id=X1, activations=X2, trades=X3)
# poisson1$type = "Poisson-1"
poisson2<-read_csv("~/c/exchange/agent-dumps/forkjoin/poisson2activations.csv",col_names=F)
poisson2 <- dplyr::rename(poisson2, id=X1, activations=X2, trades=X3)
poisson2$type = "Poisson-Poor"
poisson3<-read_csv("~/c/exchange/agent-dumps/forkjoin/poisson3activations.csv",col_names=F)
poisson3 <- dplyr::rename(poisson3, id=X1, activations=X2, trades=X3)
poisson3$type = "Poisson-Rich"
poisson4<-read_csv("~/c/exchange/agent-dumps/forkjoin/poisson4activations.csv",col_names=F)
poisson4 <- dplyr::rename(poisson4, id=X1, activations=X2, trades=X3)
poisson4$type = "Poisson-Middle"


dat <- rbind(fixed,random,uniform,poisson2,poisson3,poisson4)
dat$type = factor(dat$type, levels = c("Fixed","Random","Uniform","Poisson-Poor","Poisson-Rich","Poisson-Middle"))
dat$activations[dat$activations == 0] <- -0.9999999999 # pseudo-count

ggplot(dat, aes(x=activations, group = type)) + 
  geom_histogram(data = subset(dat,type=="Fixed"), fill = "red", alpha = 0.5, bins= 30) + 
  geom_histogram(data = subset(dat,type=="Random"), fill = "blue", alpha = 0.5, bins= 30) + 
  geom_histogram(data = subset(dat,type=="Uniform"), fill = "green", alpha = 0.5, bins= 30) + 
  geom_histogram(data = subset(dat,type=="Poisson-1"), fill = "yellow", alpha = 0.5, bins= 30) + 
  geom_histogram(data = subset(dat,type=="Poisson-2"), fill = "orange", alpha = 0.5, bins= 30) + 
  geom_histogram(data = subset(dat,type=="Poisson-3"), fill = "brown", alpha = 0.5, bins= 30) + 
  geom_histogram(data = subset(dat,type=="Poisson-4"), fill = "black", alpha = 0.5, bins= 30) +
  scale_x_log10() + 
  theme_classic(base_family = "serif", base_size = 12)

p <- ggplot(dplyr::filter(dat, as.numeric(type) < 4), aes(x=activations)) + 
  geom_histogram() +
  facet_wrap(~type,ncol = 3) + 
  #scale_x_log10(na.value = 0.000001) + 
  #annotation_logticks(sides="b", size = 0.1) +
  theme_classic(base_family = "serif", base_size = 12) +
  xlab ("Number of agent activations") + 
  ylab("Count") + 
  ggtitle("Distribution of Agent Activations")

p2 <- ggplot(dplyr::filter(dat, as.numeric(type) >= 4), aes(x=activations)) + 
  geom_histogram() +
  facet_wrap(~type, ncol = 3) + 
  scale_x_log10(na.value = 0.000001) + 
  annotation_logticks(sides="b", size = 0.1) +
  theme_classic(base_family = "serif", base_size = 12) +
  xlab ("Number of agent activations") + 
  ylab("Count")

p3 <- ggplot(dplyr::filter(dat, as.numeric(type) < 4), aes(x=trades)) + 
  geom_histogram() +
  facet_wrap(~type,ncol = 3) + 
  #scale_x_log10(na.value = 0.000001) + 
  #annotation_logticks(sides="b", size = 0.1) +
  theme_classic(base_family = "serif", base_size = 12) +
  xlab ("Number of agent trades") + 
  ylab("Count") +
  ggtitle("Distribution of Agent Trades")

p4 <- ggplot(dplyr::filter(dat, as.numeric(type) >= 4), aes(x=trades)) + 
  geom_histogram() +
facet_wrap(~type, ncol = 3) + 
  scale_x_log10(na.value = 0.000001) + 
  annotation_logticks(sides="b", size = 0.1) +
  theme_classic(base_family = "serif", base_size = 12) +
  xlab ("Number of agent trades") + 
  ylab("Count")

setwd("~/Dropbox/Thesis/images")
ggsave("forkjoin-agent-dump-histogram.pdf", p, width = 6.2, height = 3, units = "in")
ggsave("forkjoin-agent-dump-histogram2.pdf", p2, width = 6.2, height = 3, units = "in")
ggsave("forkjoin-agent-dump-histogram3.pdf", p3, width = 6.2, height = 3, units = "in")
ggsave("forkjoin-agent-dump-histogram4.pdf", p4, width = 6.2, height = 3, units = "in")

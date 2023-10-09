library(ggplot2)
library(viridis)
library(stats)
library(signs)
library(NSM3)

theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank())

#Power for one-way ANOVA and Kruskal-Wallis test for a normally assumed distribution.

Fpower_normalsample <- function(effect, grpSize, numTrial) {
 
  rejRule <- rep(0, 2)
  
  for (i in 1:numTrial) {
    obs <- c()
    for (j in 1:length(grpSize)) {
      obs <- c(obs, effect[j] + rnorm(grpSize[j]))
    }
    group <- rep(1:length(grpSize), grpSize)
    simData <- data.frame(obs, group)
    
    rejRule[1] <- rejRule[1] + (summary(aov(obs~as.factor(group), data = simData))[[1]][["Pr(>F)"]][1] < 0.05)
    
    rejRule[2] <- rejRule[2] + (kruskal.test(obs~as.factor(group), data = simData)$p.value < 0.05)
    
  }
  
  return(rejRule/numTrial)
}

#Initialize a matrix for empirical powers
powerF_normal <- matrix(rep(0, 21*2), nrow = 21)
#Compute power for several effect sizes between -2 and 2.
for (i in 0:20) {
  powerF_normal[i + 1, ] <- Fpower_normalsample((i-10)/5*c(0, 1, 2, 3), c(5, 5, 5, 5), 10000)
}
#Put the results in a data frame.
powerF_normal <- as.data.frame(powerF_normal)
#Name your column

 colnames(powerF_normal) <- c("ANOVA", "Kruskal_Wallis")
#Create a list of the effect sizes.
powerF_normal$effect <- seq(-2, 2, by = 0.2)
#Create a plot to summarize the empirical powers for the one-way ANOVA test.
ggplot(powerF_normal, aes(effect)) +
  geom_line(aes(y = ANOVA, color = "One-Way ANOVA")) +
   geom_line(aes(y = Kruskal_Wallis, color = "Kruskal-Wallis")) +
  geom_point(aes(y = ANOVA, shape = "One-Way ANOVA", color = "One-Way ANOVA")) +
   geom_point(aes(y = Kruskal_Wallis, shape = "Kruskal-Wallis", color = "Kruskal-Wallis")) + 
  labs(x = "Value of Scaling factor",
       y = "Empirical Power",
       title = expression(paste(Empirical~Power~of~One, "-", Way~ANOVA~and~the~Kruskal, "-", Wallis~Test~Against~Effect~Size, 
                                ", " , Normal~Error))) +
  scale_x_continuous(breaks = seq(-2, 2, by = 0.4), expand = c(0, 0.1), 
                     labels = signs_format(accuracy = .1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), expand = c(0, 0), limits = c(-0.025, 1.025)) +
  scale_shape_manual(name = "Test", breaks = c("One-Way ANOVA", "Kruskal-Wallis"),
                     values = c("One-Way ANOVA" = 0, "Kruskal-Wallis" = 1)) +
  scale_color_manual(name = "Test", breaks = c("One-Way ANOVA", "Kruskal-Wallis"),
                     values = c("One-Way ANOVA" = viridis(1, begin = 0.75), 
                                "Kruskal-Wallis" = viridis(1, begin = 0.25)))

#####################################################################################################################

#Power for one-way ANOVA and Kruskal-Wallis test for a cauchy assumed distribution.

Fpower_CauchySample <- function(effect, grpSize, numTrial) {
  
  rejRule <- rep(0, 2)
  
  for (i in 1:numTrial) {
    obs <- c()
    for (j in 1:length(grpSize)) {
      obs <- c(obs, effect[j] + rcauchy(grpSize[j]))
    }
    group <- rep(1:length(grpSize), grpSize)
    simData <- data.frame(obs, group)
    
    rejRule[1] <- rejRule[1] + (summary(aov(obs~as.factor(group), data = simData))[[1]][["Pr(>F)"]][1] < 0.05)
    
    rejRule[2] <- rejRule[2] + (kruskal.test(obs~as.factor(group), data = simData)$p.value < 0.05)
    
  }
  
  return(rejRule/numTrial)
}



#Initialize a matrix for empirical powers
powerF_cauchy <- matrix(rep(0, 41*2), nrow = 41)
#Compute power for several effect sizes between -2 and 2.
for (i in 0:40) {
  powerF_cauchy[i + 1, ] <- Fpower_CauchySample((i-20)/2.5*c(0, 1, 2, 3), c(5, 5, 5, 5), 10000)
}
#Put the results in a data frame.
powerF_cauchy <- as.data.frame(powerF_cauchy)
#Name your column
colnames(powerF_cauchy) <- c("ANOVA", "Kruskal_Wallis")
#Create a list of the effect sizes
powerF_cauchy$effect <- seq(-8, 8, by = 0.4)
#Create a plot to summarize the empirical powers for the one-way ANOVA test.
ggplot(powerF_cauchy, aes(effect)) +
  geom_line(aes(y = ANOVA, color = "One-Way ANOVA")) +
   geom_line(aes(y = Kruskal_Wallis, color = "Kruskal-Wallis")) +
  geom_point(aes(y = ANOVA, shape = "One-Way ANOVA", color = "One-Way ANOVA")) +
   geom_point(aes(y = Kruskal_Wallis, shape = "Kruskal-Wallis", color = "Kruskal-Wallis")) + 
  labs(x = "Value of Scaling Factor",
       y = "Empirical Power",
       title = expression(paste(Empirical~Power~of~One, "-", Way~ANOVA~and~the~Kruskal, "-", Wallis~Test~Against~Effect~Size, 
                                ", " , Cauchy~Error))) +
  scale_x_continuous(breaks = seq(-8, 8, by = 1.6), expand = c(0, 0.4), 
                     labels = signs_format(accuracy = .1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), expand = c(0, 0), limits = c(-0.025, 1.025)) +
  scale_shape_manual(name = "Test", breaks = c("One-Way ANOVA", "Kruskal-Wallis"),
                     values = c("One-Way ANOVA" = 0, "Kruskal-Wallis" = 1)) +
  scale_color_manual(name = "Test", breaks = c("One-Way ANOVA", "Kruskal-Wallis"),
                     values = c("One-Way ANOVA" = viridis(1, begin = 0.75), 
                                "Kruskal-Wallis" = viridis(1, begin = 0.25)))

##############################################################################################################

#Power for one-way ANOVA and Kruskal-Wallis test for a larger sample size.

Fpower_LargeSample <- function(effect, grpSize, numTrial) {
  
  rejRule <- rep(0, 2)
  
  for (i in 1:numTrial) {
    obs <- c()
    for (j in 1:length(grpSize)) {
      obs <- c(obs, effect[j] + rnorm(grpSize[j]))
    }
    group <- rep(1:length(grpSize), grpSize)
    simData <- data.frame(obs, group)
    
    rejRule[1] <- rejRule[1] + (summary(aov(obs~as.factor(group), data = simData))[[1]][["Pr(>F)"]][1] < 0.05)
    
    rejRule[2] <- rejRule[2] + (kruskal.test(obs~as.factor(group), data = simData)$p.value < 0.05)
    
  }
  
  return(rejRule/numTrial)
}

#Initialize a matrix for empirical powers
powerF_Large <- matrix(rep(0, 21*2), nrow = 21)
#Compute power for several effect sizes between -2 and 2.
for (i in 0:20) {
  powerF_Large[i + 1, ] <- Fpower_LargeSample((i-10)/5*c(0, 1, 2, 3), c(10, 10, 10, 10), 10000)
}
#Put the results in a data frame.
powerF_Large <- as.data.frame(powerF_Large)
#Name your column
colnames(powerF_Large) <- c("ANOVA", "Kruskal_Wallis")
#Create a list of the effect sizes.
powerF_Large$effect <- seq(-2, 2, by = 0.2)
#Create a plot to summarize the empirical powers for the one-way ANOVA test.
ggplot(powerF_Large, aes(effect)) +
geom_line(aes(y = ANOVA, color = "One-Way ANOVA")) +
 geom_line(aes(y = Kruskal_Wallis, color = "Kruskal-Wallis")) +
  geom_point(aes(y = ANOVA, shape = "One-Way ANOVA", color = "One-Way ANOVA")) +
 geom_point(aes(y = Kruskal_Wallis, shape = "Kruskal-Wallis", color = "Kruskal-Wallis")) + 
  labs(x = "Value of Scaling Factor",
       y = "Empirical Power",
       title = expression(paste(Empirical~Power~of~One, "-", Way~ANOVA~and~the~Kruskal, "-", Wallis~Test~Against~Effect~Size, 
                                ",",~Large~Samples))) +
  scale_x_continuous(breaks = seq(-2, 2, by = 0.2), expand = c(0, 0.1), 
                     labels = signs_format(accuracy = .1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), expand = c(0, 0), limits = c(-0.025, 1.025)) +
  scale_shape_manual(name = "Test", breaks = c("One-Way ANOVA", "Kruskal-Wallis"),
                     values = c("One-Way ANOVA" = 0, "Kruskal-Wallis" = 1)) +
  scale_color_manual(name = "Test", breaks = c("One-Way ANOVA", "Kruskal-Wallis"),
                     values = c("One-Way ANOVA" = viridis(1, begin = 0.75), 
                                "Kruskal-Wallis" = viridis(1, begin = 0.25)))

##############################################################################################################

#Power for one-way ANOVA and Kruskal-Wallis test for non-constant variance.


Fpower_nonconstantVar <- function(effect, grpSize, numTrial) {
  
  rejRule <- rep(0, 2)
  
  for (i in 1:numTrial) {
    obs <- c()
    for (j in 1:length(grpSize)) {
      obs <- c(obs, effect[j] + rnorm(grpSize[j], sd=sqrt(j)))
    }
    group <- rep(1:length(grpSize), grpSize)
    simData <- data.frame(obs, group)
    
    rejRule[1] <- rejRule[1] + (summary(aov(obs~as.factor(group), data = simData))[[1]][["Pr(>F)"]][1] < 0.05)
    
    rejRule[2] <- rejRule[2] + (kruskal.test(obs~as.factor(group), data = simData)$p.value < 0.05)
    
  }
  
  return(rejRule/numTrial)
}


 
#Initialize a matrix for empirical powers
powerF_nonConstantvar <- matrix(rep(0, 21*2), nrow = 21)
#Compute power for several effect sizes between -2 and 2.
for (i in 0:20) {
  powerF_nonConstantvar[i + 1, ] <- Fpower_nonconstantVar((i-10)/5*c(0, 1, 2, 3), c(5, 5, 5, 5), 10000)
}
#Put the results in a data frame.
powerF_nonConstantvar <- as.data.frame(powerF_nonConstantvar)
#Name your column
colnames(powerF_nonConstantvar) <- c("ANOVA", "Kruskal_Wallis")
#Create a list of the effect sizes.
powerF_nonConstantvar$effect <- seq(-2, 2, by = 0.2)
#Create a plot to summarize the empirical powers for the one-way ANOVA test.
ggplot(powerF_nonConstantvar, aes(effect)) +
  geom_line(aes(y = ANOVA, color = "One-Way ANOVA")) +
  geom_line(aes(y = Kruskal_Wallis, color = "Kruskal-Wallis")) +
  geom_point(aes(y = ANOVA, shape = "One-Way ANOVA", color = "One-Way ANOVA")) +
   geom_point(aes(y = Kruskal_Wallis, shape = "Kruskal-Wallis", color = "Kruskal-Wallis")) + 
  labs(x = "Value of Scaling Factor",
       y = "Empirical Power",
       title = expression(paste(Empirical~Power~of~One, "-", Way~ANOVA~and~the~Kruskal, "-", Wallis~Test~Against~Effect~Size, 
                                ", " , non-cons~Var))) +
  scale_x_continuous(breaks = seq(-2, 2, by = 0.2), expand = c(0, 0.1), 
                     labels = signs_format(accuracy = .1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), expand = c(0, 0), limits = c(-0.025, 1.025)) +
  scale_shape_manual(name = "Test", breaks = c("One-Way ANOVA", "Kruskal-Wallis"),
                     values = c("One-Way ANOVA" = 0, "Kruskal-Wallis" = 1)) +
  scale_color_manual(name = "Test", breaks = c("One-Way ANOVA", "Kruskal-Wallis"),
                     values = c("One-Way ANOVA" = viridis(1, begin = 0.75), 
                                "Kruskal-Wallis" = viridis(1, begin = 0.25)))

##############################################################################################################

#Power for one-way ANOVA and Kruskal-Wallis test for Unbalanced studies.


Fpower_unbalanced <- function(effect, grpSize, numTrial) {
  
  rejRule <- rep(0, 2)
  
  for (i in 1:numTrial) {
    obs <- c()
    for (j in 1:length(grpSize)) {
      obs <- c(obs, effect[j] + rnorm(grpSize[j]))
    }
    group <- rep(1:length(grpSize), grpSize)
    simData <- data.frame(obs, group)
    
    rejRule[1] <- rejRule[1] + (summary(aov(obs~as.factor(group), data = simData))[[1]][["Pr(>F)"]][1] < 0.05)
    
    rejRule[2] <- rejRule[2] + (kruskal.test(obs~as.factor(group), data = simData)$p.value < 0.05)
    
  }
  
  return(rejRule/numTrial)
}



#Initialize a matrix for empirical powers
powerF_unbal <- matrix(rep(0, 21*2), nrow = 21)
#Compute power for several effect sizes between -2 and 2.
for (i in 0:20) {
  powerF_unbal[i + 1, ] <- Fpower_unbalanced((i-10)/5*c(0, 1, 2, 3), c(2, 3, 4, 5), 10000)
}
#Put the results in a data frame.
powerF_unbal <- as.data.frame(powerF_unbal)
#Name your column
colnames(powerF_unbal) <- c("ANOVA", "Kruskal_Wallis")
#Create a list of the effect sizes.
powerF_unbal$effect <- seq(-2, 2, by = 0.2)
#Create a plot to summarize the empirical powers for the one-way ANOVA test.
ggplot(powerF_unbal, aes(effect)) +
  geom_line(aes(y = ANOVA, color = "One-Way ANOVA")) +
  geom_line(aes(y = Kruskal_Wallis, color = "Kruskal-Wallis")) +
  geom_point(aes(y = ANOVA, shape = "One-Way ANOVA", color = "One-Way ANOVA")) +
  geom_point(aes(y = Kruskal_Wallis, shape = "Kruskal-Wallis", color = "Kruskal-Wallis")) + 
  labs(x = "Value of Scaling Factor",
       y = "Empirical Power",
       title = expression(paste(Empirical~Power~of~One, "-", Way~ANOVA~and~the~Kruskal, "-", Wallis~Test~Against~Effect~Size, 
                                ", " , unbalanced))) +
  scale_x_continuous(breaks = seq(-2, 2, by = 0.2), expand = c(0, 0.1), 
                     labels = signs_format(accuracy = .1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), expand = c(0, 0), limits = c(-0.025, 1.025)) +
  scale_shape_manual(name = "Test", breaks = c("One-Way ANOVA", "Kruskal-Wallis"),
                     values = c("One-Way ANOVA" = 0, "Kruskal-Wallis" = 1)) +
  scale_color_manual(name = "Test", breaks = c("One-Way ANOVA", "Kruskal-Wallis"),
                     values = c("One-Way ANOVA" = viridis(1, begin = 0.75), 
                                "Kruskal-Wallis" = viridis(1, begin = 0.25)))

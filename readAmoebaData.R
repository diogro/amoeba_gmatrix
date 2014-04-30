library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(lme4)
library(MCMCglmm)
library(gridExtra)

dicty_Phen = read.csv("./dicty_phenotypes.csv")
dicty_Phen_std = read.csv("./dicty_phenotypes_standardized.csv")
dicty_Phen_std$Strain = paste("X", dicty_Phen_std$Strain, sep = '')
dicty_Phen$Strain = paste("X", dicty_Phen$Strain, sep = '')
names(dicty_Phen_std)[[4]] = 'variable'
names(dicty_Phen)[[4]] = 'variable'

dicty_Phen_std$ID = NULL
dicty_Phen$ID = NULL
dicty_Phen_std = dicty_Phen_std[complete.cases(dicty_Phen_std),]
dicty_Phen = dicty_Phen[complete.cases(dicty_Phen),]

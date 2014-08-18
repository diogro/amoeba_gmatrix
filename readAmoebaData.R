library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(lme4)
library(MCMCglmm)
library(gridExtra)
library(gtools)

dicty_Phen = read.csv("./data/dicty\ phenotypes_new.csv", as.is = T)
#dicty_Phen_std = read.csv("./dicty_phenotypes_standardized.csv")
#dicty_Phen_std$strain = paste("X", dicty_Phen_std$strain, sep = '')
dicty_Phen$strain = paste("X", dicty_Phen$strain, sep = '')
#names(dicty_Phen_std)[[4]] = 'variable'
names(dicty_Phen)[[3]] = 'variable'

dicty_Phen$variable = gsub('succ', 'succes', dicty_Phen$variable)
dicty_Phen$variable = gsub('TSC', 'tsc', dicty_Phen$variable)

dicty_Phen = filter(dicty_Phen, (value > 0 & variable == 'viab') | variable != 'viab')
dicty_Phen$value[dicty_Phen$variable == 'viab'] = logit(dicty_Phen$value[dicty_Phen$variable == 'viab'])
#dicty_Phen$value[dicty_Phen$variable == 'length'] = log10(dicty_Phen$value[dicty_Phen$variable == 'length'])

#dicty_Phen = dicty_Phen[!(dicty_Phen$variable == 'tsc' & dicty_Phen$value > 12),]

#par(mfrow= c(2, 2))
#qqnorm(dicty_Phen$value[dicty_Phen$variable == 'viab'], main = 'Viability')
#qqline(dicty_Phen$value[dicty_Phen$variable == 'viab'])
#qqnorm(dicty_Phen$value[dicty_Phen$variable == 'tsc'], main = 'Spore number')
#qqline(dicty_Phen$value[dicty_Phen$variable == 'tsc'])
#qqnorm(dicty_Phen$value[dicty_Phen$variable == 'size'], main = 'Spore size')
#qqline(dicty_Phen$value[dicty_Phen$variable == 'size'])
#qqnorm(dicty_Phen$value[dicty_Phen$variable == 'succes'], main = 'Success')
#qqline(dicty_Phen$value[dicty_Phen$variable == 'succes'])

dicty_Phen_std = dicty_Phen
dicty_Phen_std$value[dicty_Phen_std$variable == 'succes'] = scale(dicty_Phen_std$value[dicty_Phen_std$variable == 'succes'])
dicty_Phen_std$value[dicty_Phen_std$variable == 'size'] = scale(dicty_Phen_std$value[dicty_Phen_std$variable == 'size'])
dicty_Phen_std$value[dicty_Phen_std$variable == 'tsc'] = scale(dicty_Phen_std$value[dicty_Phen_std$variable == 'tsc'])
dicty_Phen_std$value[dicty_Phen_std$variable == 'viab'] = scale(dicty_Phen_std$value[dicty_Phen_std$variable == 'viab'])

ggplot(dicty_Phen_std, aes(value, valueA, color  = variable)) + geom_point()

dicty_Phen$variable = factor(dicty_Phen$variable)
dicty_Phen$strain = factor(dicty_Phen$strain)

dicty_Phen_std$variable = factor(dicty_Phen_std$variable)
dicty_Phen_std$strain = factor(dicty_Phen_std$strain)

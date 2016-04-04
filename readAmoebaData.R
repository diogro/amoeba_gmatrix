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
#dicty_Phen$strain = paste("X", dicty_Phen$strain, sep = '')
dicty_Phen$strain = as.character(dicty_Phen$strain)
#names(dicty_Phen_std)[[4]] = 'variable'
names(dicty_Phen)[[3]] = 'variable'

dicty_Phen$variable = gsub('succ', 'succes', dicty_Phen$variable)
dicty_Phen$variable = gsub('TSC', 'tsc', dicty_Phen$variable)

#dicty_Phen$value[dicty_Phen$variable == 'tsc'] = dicty_Phen$value[dicty_Phen$variable == 'tsc'] + 6.9
#dicty_Phen$value[dicty_Phen$variable == 'size'] = dicty_Phen$value[dicty_Phen$variable == 'size'] + log10(1000/40.535)

#par(mfrow= c(2, 2))
#qqnorm(dicty_Phen$value[dicty_Phen$variable == 'viab'], main = 'Viability')
#qqline(dicty_Phen$value[dicty_Phen$variable == 'viab'])
#qqnorm(dicty_Phen$value[dicty_Phen$variable == 'tsc'], main = 'Spore number')
#qqline(dicty_Phen$value[dicty_Phen$variable == 'tsc'])
#qqnorm(dicty_Phen$value[dicty_Phen$variable == 'size'], main = 'Spore size')
#qqline(dicty_Phen$value[dicty_Phen$variable == 'size'])
#qqnorm(dicty_Phen$value[dicty_Phen$variable == 'succes'], main = 'Proportion in chimera')
#qqline(dicty_Phen$value[dicty_Phen$variable == 'succes'])

variable_names <- list('viab' = 'Viability',
                       'tsc'= 'Spore number',
                       'size' = 'Spore size',
                       'succes' = 'Proportion in chimera')

dicty_Phen$plot_variable <- unlist(variable_names[dicty_Phen$variable])

boxplots <- ggplot(dicty_Phen, aes(strain, value)) + geom_boxplot() + facet_wrap(~plot_variable, scale = 'free') + theme_bw() + labs(x = 'Strain ID (See Supplementary Table)', y = 'Trait value (see legend)')+ theme(panel.margin = unit(2, "lines")) 
mask = dicty_Phen$variable == 'viab' & dicty_Phen$strain %in% c('80.1', '87.1')
violinplots <- ggplot(dicty_Phen, aes(strain, value))  + geom_violin(data = dicty_Phen[!mask,]) + geom_point() + facet_wrap(~plot_variable, scale = 'free') + theme_bw() + labs(x = 'Strain', y = 'Scaled trait value')+ theme(panel.margin = unit(2, "lines")) 
ggsave('./figures/dicty_boxplots.tiff', boxplots, height = 10, width = 20, dpi = 500)
ggsave('./figures/dicty_violinplots.tiff', violinplots, height = 10, width = 20, dpi = 500)

dicty_Phen_std = dicty_Phen
dicty_Phen_std$value[dicty_Phen_std$variable == 'succes'] = scale(dicty_Phen_std$value[dicty_Phen_std$variable == 'succes'])
dicty_Phen_std$value[dicty_Phen_std$variable == 'size']   = scale(dicty_Phen_std$value[dicty_Phen_std$variable == 'size'])
dicty_Phen_std$value[dicty_Phen_std$variable == 'tsc']    = scale(dicty_Phen_std$value[dicty_Phen_std$variable == 'tsc'])
dicty_Phen_std$value[dicty_Phen_std$variable == 'viab']   = scale(dicty_Phen_std$value[dicty_Phen_std$variable == 'viab'])

ggplot(dicty_Phen_std, aes(value, valueA, color  = variable)) + geom_point() + facet_wrap(~variable)

dicty_Phen$variable = factor(dicty_Phen$variable)
dicty_Phen$strain = factor(dicty_Phen$strain)

dicty_Phen_std$variable = factor(dicty_Phen_std$variable)
dicty_Phen_std$strain = factor(dicty_Phen_std$strain)

succes_plot = ggplot(filter(dicty_Phen, variable == 'succes'), aes(value+51)) + geom_histogram() + theme_classic() + labs(x = 'Proportion in chimera')

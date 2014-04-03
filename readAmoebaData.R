library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(lme4)
library(MCMCglmm)
library(gridExtra)

dicty_Phen = read.csv("./dicty_phenotypes.csv")
dicty_Phen_std = read.csv("./dicty_phenotypes_standardized.csv")
relat_matrix = as.matrix(read.csv("./relatedness_matrices_one_locus.csv",  header = TRUE)[-1])
rownames(relat_matrix) = colnames(relat_matrix)
relat_matrix_allLoci = as.matrix(read.csv("./relatedness_matrices_all_loci.csv",  header = TRUE)[-1])
rownames(relat_matrix_allLoci) = colnames(relat_matrix_allLoci)
dicty_Phen_std$Strain = paste("X", dicty_Phen_std$Strain, sep = '')
dicty_Phen$Strain = paste("X", dicty_Phen$Strain, sep = '')
names(dicty_Phen_std)[[4]] = 'variable'
names(dicty_Phen)[[4]] = 'variable'

dicty_Phen_std$ID = NULL
dicty_Phen_std = dicty_Phen_std[complete.cases(dicty_Phen_std),]

ggplot(dicty_Phen_std, aes(variable, value, group = interaction(Strain, variable), color=Strain)) + geom_boxplot() + theme_classic()

model = lmer(value ~ variable + (0 + variable|Strain), data = dicty_Phen_std_c)
G_lme4 = VarCorr(model)[[1]]
rownames(G_lme4) = colnames(G_lme4) = gsub('variable', '', rownames(G_lme4))
dimnames(attr(G_lme4, 'correlation')) = dimnames(G_lme4)
names(attr(G_lme4, 'stddev')) = rownames(G_lme4)

num_traits = length(unique(dicty_Phen_std$variable))
num_strains = length(unique(dicty_Phen_std$Strain))
prior = list(R = list(V = 1, n = 0.002),
             G = list(G1 = list(V = diag(num_traits) * 0.02, n = num_traits+1)))
mcmc_model = MCMCglmm(value ~ variable - 1,
                      random = ~us(variable):Strain,
                      data = dicty_Phen_std,
                      prior = prior,
                      pl = TRUE,
                      family = "gaussian")
summary(mcmc_model)
Gs = array(mcmc_model$VCV[,1:(num_traits*num_traits)], dim = c(1000, num_traits, num_traits))
G_mcmc = apply(Gs, 2:3, mean)
G_mcmc_conf = apply(Gs, 2:3, quantile, c(0.025, 0.975))
G_mcmc_conf = aperm(G_mcmc_conf, c(2, 3, 1))
containsZero = function(x) ifelse(0 > x[1] & 0 < x[2], TRUE, FALSE)
significant = !aaply(G_mcmc_conf, 1:2, containsZero)
dimnames(significant) = dimnames(G_mcmc) = dimnames(G_mcmc_conf)[1:2] = dimnames(G_lme4)
HPDinterval(mcmc_model$Sol)

sim_strains = adply(1:1000, 1, function(index) mvtnorm::rmvnorm(1, mcmc_model$Sol[index,], Gs[index,,]))
names(sim_strains) = gsub('variable', '', names(sim_strains))
cast.phen = dcast(dicty_Phen_std, Strain~variable, function(x) mean(x, na.rm = T))
plot11 = ggplot(sim_strains, aes(length, succes, group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic() + geom_point(data = cast.phen, aes(length , succes , group = 1), color = 'red')
plot12 = ggplot(sim_strains, aes(length, tsc   , group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic() + geom_point(data = cast.phen, aes(length , tsc    , group = 1), color = 'red')
plot13 = ggplot(sim_strains, aes(length, viab  , group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic() + geom_point(data = cast.phen, aes(length , viab   , group = 1), color = 'red')
plot21 = ggplot(sim_strains, aes(tsc   , succes, group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic() + geom_point(data = cast.phen, aes(tsc    , succes , group = 1), color = 'red')
plot22 = ggplot(sim_strains, aes(tsc   , viab  , group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic() + geom_point(data = cast.phen, aes(tsc    , viab   , group = 1), color = 'red')
plot23 = ggplot(sim_strains, aes(viab  , succes, group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic() + geom_point(data = cast.phen, aes(viab   , succes , group = 1), color = 'red')
png("~/Desktop/sim.png", heigh = 720, width = 1080)
grid.arrange(plot11, plot12, plot13, plot21, plot22, plot23, ncol = 3)
dev.off()

plot11 = ggplot(cast.phen, aes(length , succes , group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic()
plot12 = ggplot(cast.phen, aes(length , tsc    , group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic()
plot13 = ggplot(cast.phen, aes(length , viab   , group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic()
plot21 = ggplot(cast.phen, aes(tsc    , succes , group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic()
plot22 = ggplot(cast.phen, aes(tsc    , viab   , group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic()
plot23 = ggplot(cast.phen, aes(viab   , succes , group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic()
png("~/Desktop/real.png", heigh = 720, width = 1080)
grid.arrange(plot11, plot12, plot13, plot21, plot22, plot23, ncol = 3)
dev.off()

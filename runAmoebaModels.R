source('./readAmoebaData.R')

ggplot(dicty_Phen_std, aes(variable, value, group = interaction(Strain, variable), color=Strain)) + geom_boxplot() + theme_classic()

cast.phen = dcast(dicty_Phen_std, Strain~variable, function(x) mean(x, na.rm = T))
plot11 = ggplot(cast.phen, aes(length , succes , group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic()
plot12 = ggplot(cast.phen, aes(length , tsc    , group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic()
plot13 = ggplot(cast.phen, aes(length , viab   , group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic()
plot21 = ggplot(cast.phen, aes(tsc    , succes , group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic()
plot22 = ggplot(cast.phen, aes(tsc    , viab   , group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic()
plot23 = ggplot(cast.phen, aes(viab   , succes , group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic()
png("~/Desktop/real.png", heigh = 720, width = 1080)
grid.arrange(plot11, plot12, plot13, plot21, plot22, plot23, ncol = 3)
dev.off()

##
# lme4 traditional mixed model with clonal design
##

model = lmer(value ~ variable + (0 + variable|Strain), data = dicty_Phen_std)
G_lme4 = VarCorr(model)[[1]]
rownames(G_lme4) = colnames(G_lme4) = gsub('variable', '', rownames(G_lme4))
dimnames(attr(G_lme4, 'correlation')) = dimnames(G_lme4)
names(attr(G_lme4, 'stddev')) = rownames(G_lme4)

##
# MCMC mixed model with gaussian priors and clonal design
##

num_traits = length(unique(dicty_Phen_std$variable))
num_strains = length(unique(dicty_Phen_std$Strain))
prior = list(R = list(V = 1, n = 0.002),
             G = list(G1 = list(V = diag(num_traits) * 0.02, n = num_traits+1)))
mcmc_model = MCMCglmm(value ~ variable - 1,
                      random = ~us(variable):Strain,
                      data = dicty_Phen_std,
                      prior = prior,
                      verbose = FALSE,
                      family = "gaussian")
summary(mcmc_model)
Gs = array(mcmc_model$VCV[,1:(num_traits*num_traits)], dim = c(1000, num_traits, num_traits))
G_mcmc = apply(Gs, 2:3, mean)
G_mcmc_conf = apply(Gs, 2:3, quantile, c(0.025, 0.975))
G_mcmc_conf = aperm(G_mcmc_conf, c(2, 3, 1))
containsZero = function(x) ifelse(0 > x[1] & 0 < x[2], TRUE, FALSE)
significant = !aaply(G_mcmc_conf, 1:2, containsZero)
dimnames(significant) = dimnames(G_mcmc) = dimnames(G_mcmc_conf)[1:2] = dimnames(G_lme4)
sim_strains = adply(1:1000, 1, function(index) mvtnorm::rmvnorm(1, mcmc_model$Sol[index,], Gs[index,,]))
names(sim_strains) = gsub('variable', '', names(sim_strains))

plot11 = ggplot(sim_strains, aes(length, succes, group = 1)) + geom_point(alpha = 0.3) + geom_smooth(method="lm", color = 'black') + theme_classic()
plot12 = ggplot(sim_strains, aes(length, tsc   , group = 1)) + geom_point(alpha = 0.3) + geom_smooth(method="lm", color = 'black') + theme_classic()
plot13 = ggplot(sim_strains, aes(length, viab  , group = 1)) + geom_point(alpha = 0.3) + geom_smooth(method="lm", color = 'black') + theme_classic()
plot21 = ggplot(sim_strains, aes(tsc   , succes, group = 1)) + geom_point(alpha = 0.3) + geom_smooth(method="lm", color = 'black') + theme_classic()
plot22 = ggplot(sim_strains, aes(tsc   , viab  , group = 1)) + geom_point(alpha = 0.3) + geom_smooth(method="lm", color = 'black') + theme_classic()
plot23 = ggplot(sim_strains, aes(viab  , succes, group = 1)) + geom_point(alpha = 0.3) + geom_smooth(method="lm", color = 'black') + theme_classic()
plot11 = plot11 + geom_point(data = cast.phen, aes(length, succes, group = 1), color = 'red')
plot12 = plot12 + geom_point(data = cast.phen, aes(length, tsc   , group = 1), color = 'red')
plot13 = plot13 + geom_point(data = cast.phen, aes(length, viab  , group = 1), color = 'red')
plot21 = plot21 + geom_point(data = cast.phen, aes(tsc   , succes, group = 1), color = 'red')
plot22 = plot22 + geom_point(data = cast.phen, aes(tsc   , viab  , group = 1), color = 'red')
plot23 = plot23 + geom_point(data = cast.phen, aes(viab  , succes, group = 1), color = 'red')
png("~/Desktop/sim.png", heigh = 720, width = 1080)
grid.arrange(plot11, plot12, plot13, plot21, plot22, plot23, ncol = 3)
dev.off()

##
# MCMC mixed model with gaussian priors and relatedness matrix
##

relat_matrix = list()
relat_matrix[['one_locus']] = as.matrix(read.csv("./relatedness_matrices_one_locus.csv",  header = TRUE)[-1])
rownames(relat_matrix[['one_locus']]) = colnames(relat_matrix[['one_locus']])
relat_matrix[['all_locus']] = as.matrix(read.csv("./relatedness_matrices_all_loci.csv",  header = TRUE)[-1])
rownames(relat_matrix[['all_locus']]) = colnames(relat_matrix[['one_locus']])
relat_matrix[['distances']] = as.matrix(read.csv("./relatedness_matrices_distances.csv",  header = TRUE)[-1])
diag(relat_matrix[['distances']]) = 1
rownames(relat_matrix[['distances']]) = colnames(relat_matrix[['one_locus']])
#relat_matrix[['allelecor']] = as.matrix(read.csv("./relatedness_matrices_allelic_correlation.csv",  header = TRUE)[-1])
#rownames(relat_matrix[['allelecor']]) = colnames(relat_matrix[['one_locus']])

library(ape)
tree = llply(relat_matrix, triangMtd)
rooted = root(tree[[3]], "X34.1", resolve.root = TRUE)
rooted <- compute.brlen(rooted, 1)
#plot(rooted)
strain_relat = inverseA (rooted, scale = FALSE)
rownames(strain_relat$Ainv) = strain_relat$node.names
num_traits = length(unique(dicty_Phen_std$variable))
num_strains = length(unique(dicty_Phen_std$Strain))
prior = list(R = list(V = 1, n = 0.002),
             G = list(G1 = list(V = diag(num_traits) * 0.02, n = num_traits+1)))
mcmc_r_model = MCMCglmm(value ~ variable - 1,
                      random = ~us(variable):Strain,
                      data = dicty_Phen_std,
                      prior = prior,
                      ginverse = list(Strain = strain_relat$Ainv),
                      verbose = FALSE,
                      family = "gaussian")
summary(mcmc_r_model)
Gs_r = array(mcmc_r_model$VCV[,1:(num_traits*num_traits)], dim = c(1000, num_traits, num_traits))
G_mcmc_r = apply(Gs_r, 2:3, mean)
G_mcmc_r_conf = apply(Gs_r, 2:3, quantile, c(0.025, 0.975))
G_mcmc_r_conf = aperm(G_mcmc_r_conf, c(2, 3, 1))
significant_r = !aaply(G_mcmc_r_conf, 1:2, containsZero)
dimnames(significant_r) = dimnames(G_mcmc_r) = dimnames(G_mcmc_r_conf)[1:2] = dimnames(G_lme4)
sim_strains_r = adply(1:1000, 1, function(index) mvtnorm::rmvnorm(1, mcmc_r_model$Sol[index,], Gs_r[index,,]))
names(sim_strains_r) = gsub('variable', '', names(sim_strains_r))

plot11 = ggplot(sim_strains_r, aes(length, succes, group = 1)) + geom_point(color = 'blue', alpha = 0.3) + geom_smooth(method="lm") + theme_classic()
plot12 = ggplot(sim_strains_r, aes(length, tsc   , group = 1)) + geom_point(color = 'blue', alpha = 0.3) + geom_smooth(method="lm") + theme_classic()
plot13 = ggplot(sim_strains_r, aes(length, viab  , group = 1)) + geom_point(color = 'blue', alpha = 0.3) + geom_smooth(method="lm") + theme_classic()
plot21 = ggplot(sim_strains_r, aes(tsc   , succes, group = 1)) + geom_point(color = 'blue', alpha = 0.3) + geom_smooth(method="lm") + theme_classic()
plot22 = ggplot(sim_strains_r, aes(tsc   , viab  , group = 1)) + geom_point(color = 'blue', alpha = 0.3) + geom_smooth(method="lm") + theme_classic()
plot23 = ggplot(sim_strains_r, aes(viab  , succes, group = 1)) + geom_point(color = 'blue', alpha = 0.3) + geom_smooth(method="lm") + theme_classic()
plot11 = plot11 + geom_point(data = cast.phen, aes(length, succes, group = 1), color = 'red')
plot12 = plot12 + geom_point(data = cast.phen, aes(length, tsc   , group = 1), color = 'red')
plot13 = plot13 + geom_point(data = cast.phen, aes(length, viab  , group = 1), color = 'red')
plot21 = plot21 + geom_point(data = cast.phen, aes(tsc   , succes, group = 1), color = 'red')
plot22 = plot22 + geom_point(data = cast.phen, aes(tsc   , viab  , group = 1), color = 'red')
plot23 = plot23 + geom_point(data = cast.phen, aes(viab  , succes, group = 1), color = 'red')
png("~/Desktop/sim_relatedness.png", heigh = 720, width = 1080)
grid.arrange(plot11, plot12, plot13, plot21, plot22, plot23, ncol = 3)
dev.off()

##
# All results in one plot
##

plot11 = ggplot(sim_strains, aes(length, succes, group = 1)) + geom_point(alpha = 0.3) + geom_smooth(method="lm", color = 'black') + theme_classic()
plot12 = ggplot(sim_strains, aes(length, tsc   , group = 1)) + geom_point(alpha = 0.3) + geom_smooth(method="lm", color = 'black') + theme_classic()
plot13 = ggplot(sim_strains, aes(length, viab  , group = 1)) + geom_point(alpha = 0.3) + geom_smooth(method="lm", color = 'black') + theme_classic()
plot21 = ggplot(sim_strains, aes(tsc   , succes, group = 1)) + geom_point(alpha = 0.3) + geom_smooth(method="lm", color = 'black') + theme_classic()
plot22 = ggplot(sim_strains, aes(tsc   , viab  , group = 1)) + geom_point(alpha = 0.3) + geom_smooth(method="lm", color = 'black') + theme_classic()
plot23 = ggplot(sim_strains, aes(viab  , succes, group = 1)) + geom_point(alpha = 0.3) + geom_smooth(method="lm", color = 'black') + theme_classic()
plot11 = plot11 + geom_point(data = sim_strains_r, aes(length, succes, group = 1), color = 'blue', alpha = 0.3) + geom_smooth(data = sim_strains_r, method="lm")
plot12 = plot12 + geom_point(data = sim_strains_r, aes(length, tsc   , group = 1), color = 'blue', alpha = 0.3) + geom_smooth(data = sim_strains_r, method="lm")
plot13 = plot13 + geom_point(data = sim_strains_r, aes(length, viab  , group = 1), color = 'blue', alpha = 0.3) + geom_smooth(data = sim_strains_r, method="lm")
plot21 = plot21 + geom_point(data = sim_strains_r, aes(tsc   , succes, group = 1), color = 'blue', alpha = 0.3) + geom_smooth(data = sim_strains_r, method="lm")
plot22 = plot22 + geom_point(data = sim_strains_r, aes(tsc   , viab  , group = 1), color = 'blue', alpha = 0.3) + geom_smooth(data = sim_strains_r, method="lm")
plot23 = plot23 + geom_point(data = sim_strains_r, aes(viab  , succes, group = 1), color = 'blue', alpha = 0.3) + geom_smooth(data = sim_strains_r, method="lm")
plot11 = plot11 + geom_point(data = cast.phen, aes(length, succes, group = 1), color = 'red')
plot12 = plot12 + geom_point(data = cast.phen, aes(length, tsc   , group = 1), color = 'red')
plot13 = plot13 + geom_point(data = cast.phen, aes(length, viab  , group = 1), color = 'red')
plot21 = plot21 + geom_point(data = cast.phen, aes(tsc   , succes, group = 1), color = 'red')
plot22 = plot22 + geom_point(data = cast.phen, aes(tsc   , viab  , group = 1), color = 'red')
plot23 = plot23 + geom_point(data = cast.phen, aes(viab  , succes, group = 1), color = 'red')
png("~/Desktop/sim_full.png", heigh = 720, width = 1080)
grid.arrange(plot11, plot12, plot13, plot21, plot22, plot23, ncol = 3)
dev.off()

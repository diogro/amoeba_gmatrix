library(ape)
library(cluster)
library(Morphometrics)
library(corrgram)

source('./readAmoebaData.R')

find_CI = function(x, prob = 0.95){
    n = length(x)
    xs = sort(x)
    nint = floor(prob*n)
    lowest_int = abs(xs[n] - xs[1])
    #print(lowest_int)
    for(i in 1:(n-nint)){
        current_int = abs(xs[i] - xs[i+nint])
        if(current_int <= lowest_int){
            lowest_int = current_int
            pos = i
        }
    }
    return(c(xs[pos], xs[pos+nint]))
}

ggplot(dicty_Phen, aes(x = Strain, y = value, group = Strain)) +
geom_boxplot() + theme_classic(base_size = 20) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
facet_wrap(~variable, ncol = 4, scale = "free_y")

cast_phen = dcast(dicty_Phen_std, Strain~variable, function(x) mean(x, na.rm = T))
cast_phen_orig = dcast(dicty_Phen, Strain~variable, function(x) mean(x, na.rm = T))
plot11 = ggplot(cast_phen, aes(length, succes, group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic(base_size = 20)
plot12 = ggplot(cast_phen, aes(length, tsc   , group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic(base_size = 20)
plot13 = ggplot(cast_phen, aes(length, viab  , group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic(base_size = 20)
plot21 = ggplot(cast_phen, aes(tsc   , succes, group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic(base_size = 20)
plot22 = ggplot(cast_phen, aes(tsc   , viab  , group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic(base_size = 20)
plot23 = ggplot(cast_phen, aes(viab  , succes, group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic(base_size = 20)
#png("./figures/real.png", heigh = 720, width = 1080)
#grid.arrange(plot11, plot12, plot13, plot21, plot22, plot23, ncol = 3)
#dev.off()

##
# lme4 traditional mixed model with clonal design
##

model = lmer(value ~ variable + (0 + variable|Strain), data = dicty_Phen_std)
G_lme4 = VarCorr(model)[[1]]
rownames(G_lme4) = colnames(G_lme4) = gsub('variable', '', rownames(G_lme4))
dimnames(attr(G_lme4, 'correlation')) = dimnames(G_lme4)
names(attr(G_lme4, 'stddev')) = rownames(G_lme4)
#corrgram(G_lme4)
MatrixCompare(G_lme4, G_mcmc)

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
corr_Gs = aaply(Gs, 1, cov2cor)
G_mcmc = apply(Gs, 2:3, mean)
corr_G = apply(corr_Gs, 2:3, mean)
G_mcmc_conf = apply(Gs, 2:3, find_CI)
G_mcmc_conf = aperm(G_mcmc_conf, c(2, 3, 1))
corr_G_mcmc_conf = apply(corr_Gs, 2:3, find_CI)
corr_G_mcmc_conf = aperm(corr_G_mcmc_conf, c(2, 3, 1))
containsZero = function(x) ifelse(0 > x[1] & 0 < x[2], TRUE, FALSE)
significant = !aaply(G_mcmc_conf, 1:2, containsZero)
dimnames(significant) = dimnames(G_mcmc) = dimnames(corr_G) = dimnames(corr_G_mcmc_conf)[1:2] = dimnames(G_mcmc_conf)[1:2] = dimnames(G_lme4)
sim_strains = adply(1:1000, 1, function(index) mvtnorm::rmvnorm(1, mcmc_model$Sol[index,], Gs[index,,]))
names(sim_strains) = gsub('variable', '', names(sim_strains))
herit = summary(mcmc_model)$Gcovariances[c(1, 6, 11, 16),1:3]
rownames(herit) = c("spore size", "success", "total spore count", "viability")
write.table(herit, "./heritabilities.csv")
corr_g = round(corr_G, 3)
upper_g = round(corr_G_mcmc_conf[,,2], 3)
lower_g = round(corr_G_mcmc_conf[,,1], 3)
out_G = rbind(corr_g, lower_g, upper_g)
#write.table(out_G, "./G_correlation.csv")
#corrgram(corr_G)

plot11 = ggplot(sim_strains, aes(length, succes, group = 1)) + geom_point(alpha = 0.3) + geom_smooth(method="lm", color = 'black') + theme_classic(base_size = 20)
plot12 = ggplot(sim_strains, aes(length, tsc   , group = 1)) + geom_point(alpha = 0.3) + geom_smooth(method="lm", color = 'black') + theme_classic(base_size = 20)
plot13 = ggplot(sim_strains, aes(length, viab  , group = 1)) + geom_point(alpha = 0.3) + geom_smooth(method="lm", color = 'black') + theme_classic(base_size = 20)
plot21 = ggplot(sim_strains, aes(tsc   , succes, group = 1)) + geom_point(alpha = 0.3) + geom_smooth(method="lm", color = 'black') + theme_classic(base_size = 20)
plot22 = ggplot(sim_strains, aes(tsc   , viab  , group = 1)) + geom_point(alpha = 0.3) + geom_smooth(method="lm", color = 'black') + theme_classic(base_size = 20)
plot23 = ggplot(sim_strains, aes(viab  , succes, group = 1)) + geom_point(alpha = 0.3) + geom_smooth(method="lm", color = 'black') + theme_classic(base_size = 20)
plot11 = plot11 + geom_point(data = cast_phen, aes(length, succes, group = 1), color = 'red', size = 3) + labs(x = 'Spore size', y = 'Success')
plot12 = plot12 + geom_point(data = cast_phen, aes(length, tsc   , group = 1), color = 'red', size = 3) + labs(x = 'Spore size', y = 'Spore number')
plot13 = plot13 + geom_point(data = cast_phen, aes(length, viab  , group = 1), color = 'red', size = 3) + labs(x = 'Spore size', y = 'Viability')
plot21 = plot21 + geom_point(data = cast_phen, aes(tsc   , succes, group = 1), color = 'red', size = 3) + labs(x = 'Spore number', y = 'Success')
plot22 = plot22 + geom_point(data = cast_phen, aes(tsc   , viab  , group = 1), color = 'red', size = 3) + labs(x = 'Spore number', y = 'Viability')
plot23 = plot23 + geom_point(data = cast_phen, aes(viab  , succes, group = 1), color = 'red', size = 3) + labs(x = 'Viability', y = 'Success')
tiff("./figures/simulated.tiff", heigh = 720, width = 1080)
grid.arrange(plot11, plot12, plot13, plot21, plot22, plot23, ncol = 3)
dev.off()

##
# Fitness model
##

mcmcVar <- function(){
    mcmc_model = MCMCglmm(value ~ variable - 1,
                          random = ~idh(variable):Strain,
                          data = dicty_Phen,
                          prior = prior,
                          verbose = FALSE,
                          family = "gaussian")
    vars = aaply(mcmc_model$VCV[,1:4], 1, function(x) outer(sqrt(x), sqrt(x)))
    cors = aaply(Gs, 1, cov2cor)
    covars = cors * vars
    return(list(means = mcmc_model$Sol[,1:4], vars = covars))
}
mean_vars = mcmcVar()
sim_phens = adply(1:1000, 1, function(index) mvtnorm::rmvnorm(1, mean_vars[[1]][index,], mean_vars[[2]][index,,]))
names(sim_phens) = gsub('variable', '', names(sim_phens))
sim_phens = mutate(sim_phens, succes = (succes+50)/100)
sim_phens = mutate(sim_phens, viab = inv.logit(viab))
sim_phens = mutate(sim_phens, fitness = succes*viab)

cast_phen_orig = mutate(cast_phen_orig, succes = (succes+50)/100)
cast_phen_orig = mutate(cast_phen_orig, viab = inv.logit(viab))
cast_phen_orig = mutate(cast_phen_orig, fitness = succes*viab)

suc_tsc_plot = ggplot(sim_strains, aes(tsc, succes, group = 1)) + geom_point(alpha = 0.3) + geom_smooth(method="lm", color = 'black') + theme_classic(base_size = 20)
via_tsc_plot = ggplot(sim_strains, aes(tsc, viab  , group = 1)) + geom_point(alpha = 0.3) + geom_smooth(method="lm", color = 'black') + theme_classic(base_size = 20)
suc_tsc_plot = suc_tsc_plot + geom_point(data = cast_phen, aes(tsc, succes, group = 1), size = 3, color = 'red') + labs(x = 'Spore number', y = 'Success')
via_tsc_plot = via_tsc_plot + geom_point(data = cast_phen, aes(tsc, viab  , group = 1), size = 3, color = 'red') + labs(x = 'Spore number', y = 'Viability')
fit_tsc_plot = ggplot(sim_phens , aes(tsc , fitness , group = 1)) + geom_point(alpha = 0.3)  + theme_classic(base_size = 20)
fit_tsc_plot = fit_tsc_plot + labs(x = 'Spore number', y = 'Fitness') +
#geom_smooth(color = 'red', span = 0.1, method = loess) +
#geom_smooth(color = 'blue', span = 0.5, method = loess) +
#geom_smooth(color = 'blue', span = 0.75, method=loess) + 
geom_smooth(color = 'black', method='lm', formula = y ~ poly(x, 2))
tiff("./figures/fitness_tsc.tiff", heigh = 500, width = 1080)
grid.arrange(suc_tsc_plot, fit_tsc_plot, via_tsc_plot, ncol = 3)
dev.off()

suc_length_plot = ggplot(sim_strains, aes(length, succes, group = 1)) + geom_point(alpha = 0.3) + geom_smooth(method="lm", color = 'black') + theme_classic(base_size = 20)
via_length_plot = ggplot(sim_strains, aes(length, viab  , group = 1)) + geom_point(alpha = 0.3) + geom_smooth(method="lm", color = 'black') + theme_classic(base_size = 20)
suc_length_plot = suc_length_plot + geom_point(data = cast_phen, aes(length, succes, group = 1), size=3, color = 'red')  + labs(x = 'Spore size', y = 'Success')
via_length_plot = via_length_plot + geom_point(data = cast_phen, aes(length, viab  , group = 1), size=3, color = 'red') + labs(x = 'Spore size', y = 'Viability')
fit_length_plot = ggplot(sim_phens , aes(length , fitness , group = 1)) + geom_point(alpha = 0.3) + theme_classic(base_size = 20)
fit_length_plot = fit_length_plot + labs(x = 'Spore size', y = 'Fitness') +
#geom_smooth(color = 'red', span = 0.1, method = loess) +
#geom_smooth(color = 'blue', span = 0.5, method = loess) +
#geom_smooth(color = 'blue', span = 0.75, method=loess) + 
geom_smooth(color = 'black', method='lm', formula = y ~ poly(x, 2))
tiff("./figures/fitness_spore_size.tiff", heigh = 500, width = 1080)
grid.arrange(suc_length_plot, fit_length_plot, via_length_plot, ncol = 3)
dev.off()

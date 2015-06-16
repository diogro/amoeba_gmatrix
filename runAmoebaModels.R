if(!require(ape)) { install.packages("ape"); library(ape) }
if(!require(cluster)) { install.packages("cluster"); library(cluster) }
if(!require(Morphometrics)) {
    if(!require(devtools))
        instal.packages("devtools")
    devtools::install_github("lem-usp/evolqg"); library(evolqg) }
if(!require(corrgram)) { install.packages("corrgram"); library(corrgram) }

source('./readAmoebaData.R')

#Function for searching smallest high probability intervals

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

# Diagnostic plot for observed data

ggplot(dicty_Phen, aes(x = strain, y = value, group = strain)) +
geom_boxplot() + theme_classic(base_size = 20) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
facet_wrap(~variable, ncol = 4, scale = "free_y")

# Plot of observed means per lineage and linear lines

cast_phen = dcast(dicty_Phen_std, strain~variable, function(x) mean(x, na.rm = T))
cast_phen_orig = dcast(dicty_Phen, strain~variable, function(x) mean(x, na.rm = T))
plot11 = ggplot(cast_phen, aes(size, succes, group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic(base_size = 20)
plot12 = ggplot(cast_phen, aes(size, tsc   , group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic(base_size = 20)
plot13 = ggplot(cast_phen, aes(size, viab  , group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic(base_size = 20)
plot21 = ggplot(cast_phen, aes(tsc   , succes, group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic(base_size = 20)
plot22 = ggplot(cast_phen, aes(tsc   , viab  , group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic(base_size = 20)
plot23 = ggplot(cast_phen, aes(viab  , succes, group = 1)) + geom_point() + geom_smooth(method="lm") + theme_classic(base_size = 20)
#png("./figures/real.png", heigh = 720, width = 1080)
#grid.arrange(plot11, plot12, plot13, plot21, plot22, plot23, ncol = 3)
#dev.off()

##
# lme4 traditional mixed model with clonal design
##

model = lmer(value ~ variable + (0 + variable|strain), data = dicty_Phen_std)
G_lme4 = VarCorr(model)[[1]]
rownames(G_lme4) = colnames(G_lme4) = gsub('variable', '', rownames(G_lme4))
dimnames(attr(G_lme4, 'correlation')) = dimnames(G_lme4)
names(attr(G_lme4, 'stddev')) = rownames(G_lme4)
#corrgram(G_lme4)
#MatrixCompare(G_lme4, G_mcmc)

##
# MCMC mixed model with gaussian priors and clonal design
##

num_traits = length(unique(dicty_Phen_std$variable))
num_strains = length(unique(dicty_Phen_std$strain))
prior = list(R = list(V = 1, n = 0.002),
             G = list(G1 = list(V = diag(num_traits) * 0.02, n = num_traits+1)))
mcmc_model = MCMCglmm(value ~ variable - 1,
                      random = ~us(variable):strain,
                      data = dicty_Phen_std,
                      prior = prior,
                      verbose = TRUE,
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
significant = !aaply(corr_G_mcmc_conf, 1:2, containsZero)
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
write.table(out_G, "./G_correlation.csv")
#corrgram(corr_G)

##
# Plotting scaled simulated lineages, fitted linear line and observed means
##

simulated_plot <- function(x, y){
    ggplot(sim_strains, aes_string(x, y, group = 1)) +
    geom_point(alpha = 0.3, size = 5) + geom_smooth(method="lm", color = 'black') +
    geom_point(data = cast_phen, aes_string(x, y, group = 1), color = 'red', size = 7) +
    theme_classic(base_size = 70) +
    scale_x_continuous(limits = c(-2.2, 2.2)) + scale_y_continuous(limits = c(-2.2, 2.2)) +
    theme(axis.title.x=element_text(vjust=-1),
          axis.title.y=element_text(vjust=2),
          plot.margin = unit(c(1,1,1,1), "cm"),
          axis.line = element_line(size=2, color = "black"))
}
plot11 = simulated_plot("size", "tsc") + labs(x = 'Spore size', y = 'Spore number')
plot12 = simulated_plot("size", "viab") + labs(x = 'Spore size', y = 'Viability')
plot13 = simulated_plot("tsc" , "viab") + labs(x = 'Spore number', y = 'Viability')
plot21 = simulated_plot("tsc" , "succes") + labs(x = 'Spore number', y = 'Proportion in chimera')
plot22 = simulated_plot("size", "succes") + labs(x = 'Spore size', y = 'Proportion in chimera')
plot23 = simulated_plot("viab", "succes") + labs(x = 'Viability', y = 'Proportion in chimera')
tiff("./figures/simulated.tiff", heigh = 2000, width = 3000)
grid.arrange(plot11, plot12, plot13, plot21, plot22, plot23, ncol = 3)
grid.text("A.", x = 0.012, y = 0.98, gp = gpar(fontsize = 70))
grid.text("D.", x = 0.012, y = 0.48, gp = gpar(fontsize = 70))
grid.text("B.", x = 0.341, y = 0.98, gp = gpar(fontsize = 70))
grid.text("E.", x = 0.341, y = 0.48, gp = gpar(fontsize = 70))
grid.text("C.", x = 0.674, y = 0.98, gp = gpar(fontsize = 70))
grid.text("F.", x = 0.674, y = 0.48, gp = gpar(fontsize = 70))
dev.off()

##
# Fitness model
##

# MCMC mixed model with gaussian priors and clonal design and original variances included

mcmcVar <- function(){
    mcmc_model = MCMCglmm(value ~ variable - 1,
                          random = ~idh(variable):strain,
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
sim_phens = mutate(sim_phens, succes         = (succes+50)/100)
sim_phens = mutate(sim_phens, viab           = inv.logit(viab))
sim_phens = mutate(sim_phens, tsc_scaled     = inv.logit(tsc))
sim_phens = mutate(sim_phens, fitness        = succes*viab)
sim_phens = mutate(sim_phens, clonal_fitness = tsc_scaled*viab)

cast_phen_orig = mutate(cast_phen_orig, succes_01      = (succes+50)/100)
cast_phen_orig = mutate(cast_phen_orig, viab_unscaled  = inv.logit(viab))
cast_phen_orig = mutate(cast_phen_orig, fitness        = succes_01*viab)
cast_phen_orig = mutate(cast_phen_orig, tsc_scaled     = inv.logit(tsc))
cast_phen_orig = mutate(cast_phen_orig, clonal_fitness = tsc_scaled*viab)

fitness_plot <- function(x, y = "fitness"){
    ggplot(sim_phens , aes_string(x , y, group = 1)) + 
    geom_point(alpha = 0.3, size = 5) +
    geom_smooth(color = 'black', method='lm', formula = y ~ poly(x, 2), size = 2) + scale_x_continuous(limits = c(-1, 1)) +
    theme_classic(base_size = 80) +
    theme(axis.title.x=element_text(vjust=-1.3),
          axis.title.y=element_text(vjust=3),
          plot.margin = unit(c(1,1,1,1), "cm"),
          axis.line = element_line(size=2, color = "black"))
}
fit_tsc_plot = fitness_plot("tsc") + labs(x = 'Spore number', y = 'Realized social fitness')
fit_size_plot = fitness_plot("size") + labs(x = 'Spore size', y = 'Realized social fitness')
tiff("./figures/fitness.tiff", heigh = 2000, width = 2500)
grid.arrange(fit_size_plot, fit_tsc_plot, ncol = 2)
grid.text("A.", x = 0.04, y = 0.98, gp = gpar(fontsize = 80))
grid.text("B.", x = 0.54, y = 0.98, gp = gpar(fontsize = 80))
dev.off()

suc_tsc_plot = ggplot(sim_strains, aes(tsc, size, group = 1)) + geom_point(alpha = 0.3) + geom_smooth(method="lm", color = 'black') + theme_classic(base_size = 20)
via_tsc_plot = ggplot(sim_strains, aes(tsc, viab  , group = 1)) + geom_point(alpha = 0.3) + geom_smooth(method="lm", color = 'black') + theme_classic(base_size = 20)
suc_tsc_plot = suc_tsc_plot + geom_point(data = cast_phen, aes(tsc, succes, group = 1), size = 3, color = 'red') + labs(x = 'Spore number', y = 'Proportion in chimera') +
scale_x_continuous(limits = c(-1.5, 1.5)) + scale_y_continuous(limits = c(-2.2, 2.2))
via_tsc_plot = via_tsc_plot + geom_point(data = cast_phen, aes(tsc, viab  , group = 1), size = 3, color = 'red') + labs(x = 'Spore number', y = 'Viability') +
scale_x_continuous(limits = c(-1.5, 1.5)) + scale_y_continuous(limits = c(-2.2, 2.2))
fit_tsc_plot = ggplot(sim_phens, aes(tsc , clonal_fitness , group = 1)) + geom_point(alpha = 0.3)  + theme_classic(base_size = 20)
fit_tsc_plot = fit_tsc_plot + labs(x = 'Spore number', y = 'Clonal fitness') +
geom_point(data = cast_phen_orig, aes(tsc   , clonal_fitness  , group = 1), color = 'red', size = 3) +
geom_smooth(color = 'black', method='lm', formula = y ~ poly(x, 2)) + scale_x_continuous(limits = c(-1, 1))
png("./figures/clonal_fitness_tsc.png", heigh = 500, width = 1080)
grid.arrange(suc_tsc_plot, fit_tsc_plot, via_tsc_plot, ncol = 3)
dev.off()
tiff("./figures/clonal_fitness_tsc.tiff", heigh = 500, width = 1080)
grid.arrange(suc_tsc_plot, fit_tsc_plot, via_tsc_plot, ncol = 3)
dev.off()

suc_size_plot = ggplot(sim_strains, aes(size, tsc, group = 1)) +
geom_point(alpha = 0.3) + geom_smooth(method="lm", color = 'black') +
theme_classic(base_size = 20)
via_size_plot = ggplot(sim_strains, aes(size, viab , group = 1)) +
geom_point(alpha = 0.3) + geom_smooth(method="lm", color = 'black') +
theme_classic(base_size = 20)
suc_size_plot = suc_size_plot + geom_point(data = cast_phen, aes(size, succes, group = 1), size=3, color = 'red') +
labs(x = 'Spore size', y = 'Spore number') +
scale_x_continuous(limits = c(-2.2, 2.2)) + scale_y_continuous(limits = c(-2.2, 2.2))
via_size_plot = via_size_plot + geom_point(data = cast_phen, aes(size, viab , group = 1), size=3, color = 'red') +
labs(x = 'Spore size', y = 'Viability') +
scale_x_continuous(limits = c(-2.2, 2.2)) + scale_y_continuous(limits = c(-2.2, 2.2))
fit_size_plot = ggplot(sim_phens , aes(size , clonal_fitness , group = 1)) + geom_point(alpha = 0.3) + theme_classic(base_size = 20)
fit_size_plot = fit_size_plot + labs(x = 'Spore size', y = 'Clonal fitness') +
geom_point(data = cast_phen_orig, aes(scale(size)   , clonal_fitness  , group = 1), color = 'red', size = 3) +
geom_smooth(color = 'black', method='lm', formula = y ~ poly(x, 2)) +
scale_x_continuous(limits = c(-2.2, 2.2))
png("./figures/clonal_fitness_spore_size.png", heigh = 500, width = 1080)
grid.arrange(suc_size_plot, fit_size_plot, via_size_plot, ncol = 3)
dev.off()
tiff("./figures/clonal_fitness_spore_size.tiff", heigh = 500, width = 1080)
grid.arrange(suc_size_plot, fit_size_plot, via_size_plot, ncol = 3)
dev.off()

write.csv(sim_strains[-1], "scaled_simulated_lineages.csv", row.names = FALSE)
write.csv(sim_phens[-1], "unscaled_simulated_lineages.csv", row.names = FALSE)
write.csv(cast_phen_orig, "mean_observed_unscaled_lineages.csv", row.names = FALSE)
write.csv(cast_phen, "mean_observed_scaled_lineages.csv", row.names = FALSE)

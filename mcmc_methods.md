##G-matrix

To estimate the genetic correlation between traits we first standardized
each trait to a mean of zero and a variance of one, while viability
was first transformed using a logit function and then scaled. The four
traits were used in a multivariate model fitted using MCMCglmm using
a Bayesian modification of the framework described by Fry (2004) to
estimate genetic variances of and genetic correlations among the four
traits measured in the 24 genotypes. In this model, the four traits
are treated as measures of the same underlying trait and the genotype
is used as the unit of repeated measurements. The model estimates the
correlations for the traits at the level of genotypes, which represent
genetic correlations since the traits were measured independently
(and therefore, correlations between different traits measured on
the same genotype must be genetic correlations, not environmental
correlations) (see Fry 2004 for further details). We used weakly
informative independent Gaussian priors for the residual and mixed
effect (genetic) variances. Model convergence was accessed by inspection
of variable traces. G-matrix was estimated as the mean of the posterior
distribution of genetic correlations. This estimate was similar to
results from a maximum likelihood fitted model using SAS or lme4, with
the added advantage that the posterior distribution of covariances
can be used to create confidence intervals for all heritabilities and
genetic correlations, taking into account all sources of uncertainty in
the system and allowing a straight forward test for significant difference
from zero for all correlations (Gelman, 2004). We show the relationship
between traits using simulated lineages, obtained using the posterior
sampled covariance matrices and means to draw values from normal
distributions in the scaled trait space. We plot these values along with
means for each observed lineage to show the agreement between data and
posterior simulations. These simulations help illustrate the uncertainty
for each of the relations between traits, and again their distribution
takes into account every source of underlying uncertainty.

##Fitness

To understand the nature of selection on social traits we first
estimated the means and variances of the unscaled traits using MCMCglmm
and, along with the correlations estimated for the scaled traits, used
them to generate posterior simulations of lineages in the unscaled
space. We constructed a theoretical multiplicative model for total
fitness calculated simply as the product of social success and spore
viability. This model assumes that genotypic fitness is determined by
social success, determined by spore counts, weighted by the viability
of the spores produced. To calculate fitness in each of the simulated
lineages we first returned viability to its original [0, 1] interval via
a inverse logit transformation and scaled success to the same interval,
representing average relative abundance for a given lineage. Fitness,
the product between viability and success, is then restricted between
zero and one. 


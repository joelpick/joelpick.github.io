---
layout: blog
title: Specifying Priors in MCMCglmm
---

Priors can be quite intimidating, and are very important. They are the key difference between frequentist and Bayesian statistics. 

Many people ask how to chose them and specify them, and many people resort to googling and copy and pasting something they find online. However, priors can have a massive impact on your analysis, and so specifying them correctly is important.

# What are Priors 
A prior distribution weights your analysis. You are essentially telling the model how likely you think different values are, and the model weights its likelihood (which comes from the information in your data) by these probabilities to produce a posterior distribution that you then use for inference. 

In the majority of cases in ecology and evolution, we don't have much idea in advance what parameter values are going to be, and so we want to make the priors uninformative - we don't want to influence our analysis. Therefore we pick distributions which are quite flat over a large range of possible values.




# Appropriate Prior Distributions

half cauchy for SD



# Priors in MCMCglmm
a nested list, which contains information for the prior distributions of the fixed effects (B), random effects (G) and residuals (R). 


Simpler models in MCMCglmm often run without prior specification. The default priors are:

## Univariate

### Fixed effects

### Random effects
Have to specify the same number of 

### Residuals 
Have to specify the same number of 

#### Parameter Expanded Priors

## Multivariate
Basically has to match the structure of your random effects

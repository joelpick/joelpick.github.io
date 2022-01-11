---
title: "Intro to linear (mixed models)"
output: 
  html_document:
fig_caption: yes
---



anova, t-test, linear regression, - all the same thing


model equation 

y = beta_0 + beta_1*x_1 + beta_2*x_2 + e

y = response variable
x = predictor variables
beta = parameters - in this case 'fixed effects'
e = residual - N(0,1)

the important bit we really need to understand is e
in someways its simple, its the bit left over that we haven't explained with everything else
however the assumptions we make about it are critical to how we test things in the model.
we call this a parametric test and this is because we make certain assumptions about how e is distributed, and this allows us to test things in the model



 mixed model

 y = beta_0 + beta_1*x_1 + beta_2*x_2 + sigma_1 *z_1 + sigma_2 *z_2 + e

z - latent variables/random effects that are N(0,1)
sigma - sd of random variable

this representation makes it clear that it is a simple extension of the linear model

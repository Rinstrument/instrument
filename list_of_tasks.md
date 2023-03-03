
Tutorial on ADVI algorithm: https://luiarthur.github.io/statorial/varinf/linregpy/.
Overview of the ADVI algorithm: https://luiarthur.github.io/statorial/varinf/introvi/.

Implement method dispatch for inirt model formulas with: 
https://stackoverflow.com/questions/49346939/multiple-dispatch-for-subset-methods-in-r
S4 hidden inside S3 to get the dispatch right

factor loadings and rotations:
https://stats.stackexchange.com/questions/228629/conversion-of-irt-logit-discrimination-parameter-to-factor-loading-metric

chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.lesahoffman.com/PSYC948/948_Lecture6_Binary_Responses.pdf

(a / sqrt(1 + a^2)) or replace 1 with 3.29 depending on the model

Tutorial of IRT in stan:
https://agabrioblog.onrender.com/tutorial/irt-stan/irt-stan/

Add anova() for model comparison: 
https://methodenlehre.github.io/SGSCLM-R-course/cfa-and-sem-with-lavaan.html


M clark has a new set of tutorials that fit the 3pl and 4pl models: 
https://m-clark.github.io/models-by-example/bayesian-irt.html#three-parameter-irt
Woo Hoo!!!!


March 1, 2023: 
 - use the clark article to get the 3pl up and running
 - confirm that the higher-order thing is working correctly
 - is lambda_ind wrong???

March 2, 2023:
 - Didn't get through a lot of these points today, but
  I think I'm putting the final touches on the second-order
  MIRT model with regression..
 - with that in mind: I need to check that regression on
 - thetag works. How should that regression be specified?
 it seems as though its one equation for all thetas potentially
 
 - may need to exagerate effects and then rescale! - sim step
 - finish writing introduction
 - Write model section
 - Simulation study? - choose a model and settings
 - Case study: asti data with covariates
 - Case study from a couple different perspectives with R code
 - Missing data in all places of model now appropriate - need to add
 - hook up the mirt model and confirmatory mirt

 - simulation study: focus on regression coefficient recovery
 - Can show some tables of measures and recovery/plots of random effects?

March 3, 2023:
 - Hook up the regression of second-order model: theta, alpha, delta
 - Fixed theta effects in second-order model
 - Last look over of second-order model
 - DONE???
 - Hook up exploratory MIRT model
 - regressions of the same

...
 April 28, 2023: 
  - Friday morning - submit revisions to journal
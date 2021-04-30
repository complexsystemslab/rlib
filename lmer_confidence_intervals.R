#!/usr/bin/Rscript

library(lme4)
library(lmerTest)

# From ?lmer:

fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

# Note that the summary is more extensive with lmerTest loaded than without:

s1 <- summary(fm1)

# Confidence intervals (works with or without lmerTest loaded):

confint(fm1)

# ... "Computing profile confidence intervals ..."

# MORE EXPLANATION:
# - https://stats.stackexchange.com/questions/117641/how-trustworthy-are-the-confidence-intervals-for-lmer-objects-through-effects-pa


# Via bootstrapping:
# https://rdrr.io/cran/lme4/man/confint.merMod.html

confint.merMod(fm1, method = "boot")

# rlib


## Introduction

Miscellaneous R support functions for data science and machine learning.

Credit Rudolf N. Cardinal

## Installation

Install R and R Studio

and 

```r

source("https://raw.githubusercontent.com/neelsoumya/rlib/master/INSTALL_MANY_MODULES.R")

```

For a restricted installation

```r

install.packages('devtools')
library(devtools)

devtools::install_github('neelsoumya/rlib')

```


## Usage

You can source the latest version from Github. For example, to "source" the miscstat.R file from R, you can do this:

```r

source("https://raw.githubusercontent.com/neelsoumya/rlib/master/miscstat.R")

```

## Example

```r

source("https://raw.githubusercontent.com/neelsoumya/rlib/master/cris_common.R")

```

see

`logistic_regression_roc_curve_withCV.R`

https://github.com/neelsoumya/rlib/blob/master/logistic_regression_roc_curve_withCV.R

cris$visualize_fixed_effects_from_lmer(name of lmer or glmer model)
  
cris$fixed_effects_from_lmer(name of lmer or glmer model)
  
  
## Files
  
  `logistic_regression_roc_curve_withCV.R` logistic regression model with CV cross validation and plotting of log-odds odds ratios and AUPR and AUC curves
  
   https://github.com/neelsoumya/rlib/blob/master/logistic_regression_roc_curve_withCV.R
  
  `survival_analysis_example.rmd` script for survival analysis
  
  `rmarkdown_analysis_template.rmd` template R markdown for reproducible analysis
  
  also see
  
  https://github.com/neelsoumya/teaching_reproducible_science_R/blob/main/rmarkdown.rmd
  
  `heatmap_script.R` script for heatmap
  
  `fn_redorder_dendrogram.R` script for reordering labels in heatmaps
  
  `convert_aupr_to_ggplot.R` convert AUPR plot to ggplot for pretty plotting
  
  `lmer_confidence_intervals.R` confidence intervals for lmer linear mixed effects model
  
  
## Acknowledgements

Rudolf N. Cardinal


## Contact
  
Soumya Banerjee and Rudolf Cardinal  
  

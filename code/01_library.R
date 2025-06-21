# source code for packages to install and load

# install.packages("pacman")

pacman::p_load(here,
               readxl,
               nlme,
               lsmeans,
               performance,
               PerformanceAnalytics,
               ggpubr,
               codyn, # data generation
               psych, # data generation
               # TaxonStand, # data generation
               WorldFlora, # data generation
               DescTools, # data generation
               MASS, # Q1 analysis
               foreign, # Q1 analysis
               Hmisc, # Q1 analysis
               nnet, # Q1 analysis
               gridExtra,
               foreach, # Q1 for each loop creation
               ggExtra,
               grid, # Q1 legend
               gridExtra, # Q1 legend
               patchwork, # Q1 combining plots 
               sf,
               ggalluvial, # Q3 sankey diagram
               tidyverse,
               pscl, # Q1 McFadden R2libra
               broom #table cleaning
               )

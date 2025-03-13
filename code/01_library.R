# source code for packages to install and load

# install.packages("pacman")
pacman::p_load(here,
               nlme,
               lsmeans,
               performance,
               PerformanceAnalytics,
               ggpubr,
               codyn, # data generation
               psych, # data generation
               # TaxonStand, # data generation
               WorldFlora, # data generation
               MASS, # Q1 analysis
               foreign, # Q1 analysis
               Hmisc, # Q1 analysis
               tidyverse) 
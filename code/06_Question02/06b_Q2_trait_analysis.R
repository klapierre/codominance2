
################################################################################
##  06b_Q2_trait_analysis.R: link codominants to trait data
##
##  Authors: Akira Terui
##  Date created: 3/26/2025
################################################################################

# PROCEDURE:
# 1. select sites with codominants with trait data
# 2. calculate species trait distance for selected sites
#   - for codominant pairs (observed)
#   - for all species pairs (null distribution)
# 3. keep how many species/sites are excluded

# setup -------------------------------------------------------------------

source("code/01_library.R")
source("code/02_functions.R")

# read data ---------------------------------------------------------------

## data for codominant species
df_codom <- readRDS("data/Q2ctlGroupsSite.rds") %>% 
  group_by(site_code,
           project_name,
           community_type) %>% 
  mutate(pair_id = row_number()) %>% 
  ungroup() %>% 
  pivot_longer(cols = starts_with("alpha"),
               names_to = "alpha",
               values_to = "species") %>% 
  drop_na(species)

## data for trait data
df_trait <- readRDS("data/allTraits.rds")

# select unique combo of site, project, community-type --------------------

df_flag <- df_codom %>% 
  distinct(species) %>% 
  left_join(df_trait) %>% 
  mutate(flag = if_any(everything(),
                       .fns = is.na))
  
         
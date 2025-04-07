
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
df_codom0 <- readRDS("data/Q2ctlGroupsSite.rds") %>% 
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
## - clonal data are removed because of high uncertainty
df_trait <- readRDS("data/allTraits.rds") %>% 
  dplyr::select(-clonal)

## data for species pool
df_pool0 <- readRDS("data/allSppList.rds") %>% 
  distinct(site_proj_comm,
           genus_species)

# select unique combo of site, project, community-type --------------------

# codominance dataframe ####
## flag species row if any traits are missing
## use "df_codom_cl" for analysis

df_flag <- df_codom0 %>% 
  distinct(species) %>% 
  left_join(df_trait) %>% 
  mutate(flag = if_any(everything(),
                       .fns = is.na)) %>% 
  select(species, flag)
  
df_codom <- df_codom0 %>%
  left_join(df_flag) %>% 
  group_by(site_code,
           project_name,
           community_type,
           pair_id) %>% 
  mutate(flag_site = any(flag)) %>% 
  ungroup()

df_codom_cl <- df_codom %>% 
  filter(!flag_site) %>% 
  select(-c(flag, flag_site))

## check numbers
(n_all <- nrow(df_codom))
(n_cl <- nrow(df_codom_cl))
(n_omit <- df_codom %>% 
  pull(flag_site) %>% 
  sum())

if (n_all != sum(c(n_cl, n_omit))) 
  stop("something wrong")

# species pool dataframe ####
## use "df_pool_cl" for analysis

df_pool <- df_pool0 %>% 
  rename(species = genus_species) %>% 
  left_join(df_trait,
            by = "species") %>% 
  mutate(flag = if_any(everything(),
                       .fns = is.na))
  
df_pool_cl <- df_pool %>% 
  group_by(site_proj_comm) %>% 
  summarize(n = n_distinct(species),
            n_na = sum(flag),
            p_na = n_na / n) %>% 
  right_join(df_pool,
             by = "site_proj_comm") %>% 
  drop_na(LDMC:n_fixation_type) %>% 
  select(-c(flag, clonal)) %>% # NOTE: excluding clonal because of "CHECK" remains there
  mutate(across(.cols = where(is.character),
                .fns = function(x) {
                  if (n_distinct(x) > 2) {
                    return(factor(x))
                  } else {
                    return(as.numeric(factor(x)))
                  }
                }))

# trait distance ----------------------------------------------------------

## pool
usite <- unique(df_pool_cl$site_proj_comm)
cnm <- colnames(df_pool_cl)

df_pool_cl %>% 
  filter(site_proj_comm == usite[1]) %>%
  select(LDMC:n_fixation_type) %>% 
  as.data.frame() %>% 
  FD::gowdis()


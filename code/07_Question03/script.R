
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

# ISSUE:
# Some pairs are duplicated in the outcome, look into codes to fix
# - duplicates appear when there are overlaps in tri-dominant and co-dominant pairs
# - removed duplicates with distinct() function

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/01_library.R")
source("code/02_functions.R")

# read data ---------------------------------------------------------------

## data for codominant species for treatment plots
## - treatment, CO2, disturbance, drought, irr, temp have few replicates - removed
df_codom0 <- readRDS("data/Q3trtGroupsSite.rds") %>% 
  filter(!is.na(alpha2),
         !(trt_type %in% c("CO2", 
                           "disturbance", 
                           "drought", 
                           "irr", 
                           "temp"))) %>% 
  group_by(site_code,
           project_name,
           community_type) %>% 
  mutate(pair_id = row_number()) %>% 
  ungroup() %>% 
  pivot_longer(cols = starts_with("alpha"),
               names_to = "alpha",
               values_to = "species") %>% 
  drop_na(species) %>% 
  mutate(site_proj_comm = paste0(site_code, "_", 
                                 project_name, "_",
                                 community_type))

## data for trait data
## - clonal data are removed because of high uncertainty
## - photosynthetic pathway needs to be fixed - some "uncertain"
df_trait <- readRDS("data/allTraits.rds") %>% 
  dplyr::select(species,
                LDMC, 
                SLA,
                SRL, 
                leaf_N,
                plant_height_vegetative,
                seed_dry_mass,
                growth_form, 
                lifespan,
                photosynthetic_pathway, 
                n_fixation_type) %>% 
  mutate(photosynthetic_pathway = case_when(photosynthetic_pathway == "uncertain" ~ "C3",
                                            photosynthetic_pathway == "possible C4" ~ "C4",
                                            photosynthetic_pathway == "possible CAM" ~ "CAM",
                                            photosynthetic_pathway == "C3-C4 Intermediate" ~ "hybrid",
                                            .default = photosynthetic_pathway))

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
  mutate(flag_pair = any(flag)) %>% 
  ungroup()

df_codom_cl <- df_codom %>% 
  filter(!flag_pair) %>% 
  select(-c(flag, flag_pair))

## check numbers
(n_all <- nrow(df_codom))
(n_cl <- nrow(df_codom_cl))
(n_omit <- df_codom %>% 
    pull(flag_pair) %>% 
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
  select(c(site_proj_comm:species, 
           LDMC:n_fixation_type)) %>% 
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
usite <- unique(df_codom_cl$site_proj_comm)

df_p_trt <- foreach(k = usite,
                    .combine = bind_rows) %do% {
                      
                      ## select one site
                      df_codom_i <- df_codom_cl %>% 
                        filter(site_proj_comm == k) %>% 
                        mutate(pair_id = as.numeric(factor(pair_id)))
                      
                      ## define species pool of the site
                      df_pool_i <- df_pool_cl %>% 
                        filter(site_proj_comm == k) %>% 
                        arrange(species)
                      
                      ## trait distance matrix for the site
                      md <- df_pool_i %>%
                        select(LDMC:n_fixation_type) %>% 
                        as.data.frame() %>% 
                        FD::gowdis() %>% 
                        data.matrix()
                      
                      ## species pool
                      pool <- df_pool_i %>% 
                        pull(species)
                      
                      ## get trait distance of the co-dominants
                      ## use `distinct()` to remove duplicates - duplicates appear when tri-dominants and co-dominants share species
                      df_dist <- get_dist(data = df_codom_i,
                                          pool = pool,
                                          md = md) %>% 
                        mutate(p_na = unique(df_pool_i$p_na),
                               trt_type = unique(df_codom_i$trt_type)) %>% 
                        distinct()
                      
                      return(df_dist)
                    }


# pair with control -------------------------------------------------------

# readRDS("data/traitp_ctr.rds") %>% 
#   left_join(df_p,
#             by = "site_proj_comm") %>% 
#   view()

df_p_trt %>% 
  ggplot(aes(x = p)) +
  geom_density()

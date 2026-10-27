
################################################################################
##  06a_Q2_formatCodominants.R: Get codominating species from codom list.
##
##  Authors: Akira Terui, Jordan Winter, Kimberly Komatsu
##  Date created: January 29, 2025
################################################################################

source("code/01_library.R")
source("code/02_functions.R")


# Read Data --------------------------------------------------------------------

codomSppList <- readRDS("data/allSppList.rds") %>% 
  left_join(readRDS("data/expInfo.rds")) %>% 
  filter(num_group!=4, # remove even communities
         rank<(num_group+1)) # filter to codominating species only
  

# Most common species combinations for each plot through time ------------------

#filter to control plots only
codomControl <- codomSppList %>% 
  filter(trt_type=='control') %>% 
  left_join(read_rds('data/modeSite.rds'))

Q2ctlGroupsSite <- codomControl %>%
  filter(num_group==mode_site) %>% #filter to only plots where the number of codominanting species matches the site pairings
  arrange(exp_unit, genus_species) %>% # arrange species within plots alphabetically
  group_by(exp_unit) %>%
  mutate(alphabetical_rank=paste('alpha', rank(genus_species), sep='')) %>%  # generate an alphabetical rank for each species group
  ungroup() %>% 
  dplyr::select(-rank) %>% 
  pivot_wider(names_from=alphabetical_rank, values_from=genus_species) %>% # generate a list of codominating species groups
  group_by(site_code, project_name, community_type, plot_id, alpha1, alpha2, alpha3) %>% 
  summarise(codom_freq_plot=length(plot_id), .groups="drop") %>% # frequency of each species group within a plot through time
  group_by(site_code, project_name, community_type, plot_id) %>% 
  filter(codom_freq_plot==max(codom_freq_plot)) %>%  # drop any species group that is not most common in a plot through time
  ungroup() %>% 
  group_by(site_code, project_name, community_type, alpha1, alpha2, alpha3) %>% 
  summarise(codom_freq_proj=length(community_type), .groups="drop") %>% # frequency of each species group within an experiment across plots
  group_by(site_code, project_name, community_type) %>% 
  filter(codom_freq_proj==max(codom_freq_proj)) %>%  # drop any species group that is not most common within an experiment across plots
  ungroup()

# saveRDS(Q2ctlGroupsSite, file = "data/Q2ctlGroupsSite.rds")



# All possible species combinations (co-occur within at least one plot at one time point) within each site

allSppList <- readRDS("data/allSppList.rds")

allSppPairs <- allSppList %>%
  distinct(site_proj_comm, plot_id, calendar_year, genus_species) %>% 
  # Generate all species pairs within each plot-year
  group_by(site_proj_comm, plot_id, calendar_year) %>%
  filter(n() >= 2) %>%
  reframe(
    species_1 = combn(sort(genus_species), 2)[1, ],
    species_2 = combn(sort(genus_species), 2)[2, ]
  ) %>%
  ungroup() %>%
  distinct(site_proj_comm, species_1, species_2)

# saveRDS(allSppPairs, file = "data/allSppPairs.rds")



# # Testing which species pairs are missing from the overall dataset ----------
# 
# # Konza subset
# knz <- allSppList %>% 
#   filter(str_detect(site_proj_comm, "KNZ"))
# 
# knzPairs <- allSppPairs %>% 
#   filter(str_detect(site_proj_comm, "KNZ"))
# 
# knzPossible <- allSppList %>% 
#   filter(str_detect(site_proj_comm, "KNZ")) %>% 
#   select(genus_species) %>% 
#   unique() %>% 
#   arrange(genus_species) %>% 
#   pull(genus_species) %>% 
#   combn(2) %>% 
#   t() %>% 
#   as.data.frame() %>% 
#   rename(species_1=V1, species_2=V2)
# 
# knzMissingPairs <- anti_join(knzPossible, knzPairs)
# 
# 
# # All projects
# 
# possible <- allSppList %>% 
#   distinct(site_proj_comm, genus_species) %>% 
#   group_by(site_proj_comm) %>% 
#   filter(n() >= 2) %>% 
#   reframe(species_1 = combn(sort(genus_species), 2)[1, ],
#           species_2 = combn(sort(genus_species), 2)[2, ])
# 
# missingPairs <- anti_join(possible, allSppPairs)
# 
# # possible pairs within sites = 892,678
# # existing pairs in dataset = 530,453 (59.4%)
# # missing pairs from possible = 362,225 (40.6%)


# Prepare data for trait analysis Q2 & Q3 ---------------------------------

df_codom0_ctr <- readRDS("data/Q2ctlGroupsSite.rds") %>% 
  mutate(trt_type = "control")

df_codom0_trt <- readRDS("data/Q3trtGroupsSite.rds")

## data for codominant species for treatment plots
## - trt_type, re-coded to trt_category %in% Resource, Stress, Other, Mult.Trts
df_codom0 <- bind_rows(df_codom0_ctr,
                       df_codom0_trt) %>%  
  mutate(trt_cat = case_when(trt_type == 'control' ~ 'control',
                             trt_type %in% c('CO2','N','P','N*P','mult_nutrient','irr','K') ~ 'resource',
                             trt_type %in% c('drought','temp','herb_removal') ~ 'stress',
                             trt_type == 'multiple_trts' ~ 'mult',
                             .default = 'other')) %>% 
  mutate(site_proj_comm = paste0(site_code, "_", 
                                 project_name, "_",
                                 community_type) %>% 
           str_to_lower(),
         trt_type = str_to_lower(trt_type),
         .before = site_code) %>% 
  dplyr::select(-c(site_code,
                   project_name,
                   community_type,
                   codom_freq_proj)) %>% 
  filter(!is.na(alpha2)) %>%
  group_by(site_proj_comm) %>% 
  mutate(pair_id = row_number()) %>% 
  ungroup() %>% 
  pivot_longer(cols = starts_with("alpha"),
               names_to = "alpha",
               values_to = "species") %>% 
  drop_na(species)

## data for trait data
## - clonal data are removed because of high uncertainty
## - photosynthetic pathway needs to be fixed - some "uncertain"
df_trait <- readRDS("data/allTraits.rds")  %>% 
  janitor::clean_names() %>%  
  dplyr::select(species,
                ldmc, 
                sla,
                srl, 
                leaf_n,
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
df_pool0 <- readRDS("data/allSppPairs.rds") %>% 
  mutate(site_proj_comm = str_to_lower(site_proj_comm)) %>% 
  group_by(site_proj_comm) %>% 
  mutate(pair_id = row_number()) %>%
  ungroup() %>% 
  pivot_longer(
    cols = starts_with("species_"),
    values_to = "species"
  ) %>% 
  select(-name)

# saveRDS(df_trait, "data/data_trait.rds")
# saveRDS(df_codom0, "data/data_codom.rds")
# saveRDS(df_pool0, "data/data_pool.rds")

################################################################################
##  06b_Q2_trait_analysis.R: link codominants to trait data
##
##  Authors: Akira Terui (modified K. Komatsu)
##  Date created: 3/26/2025
################################################################################

# PROCEDURE:
# 1. select sites with codominants with trait data
# 2. calculate species trait distance for selected sites
#   - for codominant pairs (observed)
#   - for all species pairs (null distribution)
# 3. keep how many species/sites are excluded

# ISSUE:
# Some pairs are duplicated in the outcome
# - duplicates appear when there are overlaps in tri-dominant and co-dominant pairs
# - removed duplicates with distinct() function

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/01_library.R")
source("code/02_functions.R")

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34, color='black'),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34, color='black'),
             plot.title = element_text(size=40, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=40))

# read data ---------------------------------------------------------------

## data for codominant species
df_codom0 <- readRDS("data/Q2ctlGroupsSite.rds") %>% 
  filter(!is.na(alpha2)) %>% 
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

# how many plants are identified to species level
completeSpp <- df_pool0 %>% 
  select(genus_species) %>% 
  unique() %>% #start with 7725 unique species identifiers
  separate(col=genus_species, into=c('genus', 'species', 'other1', 'other2',
                                     'other3', 'other4', 'other5', 'other6',
                                     'other7', 'other10'), sep=' ', remove=F) %>%
  filter(!(species %in% c('spp', 'SP.', 'sp', 'species', '-', '.', 
                          '1', '2', '3', 'seedling')), 
         !is.na(species), 
         !(genus %in% c('unknown', 'unk', 'ZZZZ', 'Unknown')), 
         !(genus_species %in% c('forb, unidentified', 'geranium pot', 
                                'rosa spec', 'pycnan vir', 'annual forb', 
                                'fleshy forb', 'forb hibiscus',
                                'cottony forb', 'forb six', 
                                'grass, unidentified', 'un lily', 
                                'unidentified dicot', 'poaceae family', 
                                'cyperaceae family', 'ericaceae family', 
                                'cryptantha blue', 'Stiff leathery', 
                                'Dangle flower', 'Aster opposite',
                                'Long thin', '')),
         is.na(other1)) %>%  #6211 species identified to species level
  select(genus_species) %>% 
  rename(species=genus_species) %>% 
  left_join(df_trait) %>% 
  na.omit() #3817 species have all traits

# test <- completeSpp %>%
#   filter(is.na(LDMC))
# 
# write.csv(test, 'C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\1_first author\\codominance\\data\\species that still need traits_full.csv')

####### NOTE: fix this here, some species probably should have traits, but are missing for some reason. Complete traits for species in the spreadsheet in onedrive, redo trait imputation, and rerun all trait-based code #########


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

## 3816 species have complete trait information
df_pool_cl <- df_pool %>% 
  group_by(site_proj_comm) %>% 
  summarize(n = n_distinct(species),
            n_na = sum(flag),
            p_na = n_na / n) %>% 
  right_join(df_pool,
             by = "site_proj_comm") %>% 
  drop_na(LDMC:n_fixation_type) %>% # remove species with any missing traits
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

df_p <- foreach(k = usite,
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
                    mutate(p_na = unique(df_pool_i$p_na)) %>% 
                    distinct()
                  
                  return(df_dist)
                }

# df_p %>% 
#   ggplot(aes(x = p)) +
#   geom_density(size = 1) +
#   xlab('Relative Deviation\nof Functional Distance') +
#   ylab('Density') +
#   theme_bw() +
#   theme(panel.grid = element_blank())

# ggsave('Fig3_traitDensityCtl.png',
#        width=3,
#        height=3,
#        units='in',
#        dpi=300,
#        bg='white')

# linking p to environmental factors --------------------------------------

## p value and environmental factors in sf
sf_m <- readRDS("data/envData.rds") %>% 
  mutate(site_proj_comm = paste0(site_code, "_", 
                                 project_name, "_",
                                 community_type)) %>% 
  right_join(df_p,
            by = "site_proj_comm") %>% 
  st_as_sf(coords = c("Longitude", "Latitude"),
           crs= 4326)

saveRDS(sf_m, file = "data/traitp_ctr.rds")

# analysis ----------------------------------------------------------------

df_m <- readRDS("data/traitp_ctr.rds") %>% 
  filter(!is.na(anpp),
         n_obs>0)

## model with environmental gradient
glmmTMB::glmmTMB(cbind(n_obs, n_pool - n_obs) ~
                   scale(MAP) +
                   scale(MAT) +
                   scale(anpp) +
                   scale(gamma_rich) +
                   scale(HumanFootprint) +
                   scale(NDeposition) +
                   (1 | site_proj_comm),
                 family = glmmTMB::betabinomial,
                 data = df_m, 
                 weights = 1 - p_na) %>% 
  summary()

df_m2  <- df_m %>% 
  select(database, site_code, project_name, community_type, site_proj_comm,
         MAP:CONTINENT, pair_id, sp1, sp2, dist, ses, p, n_pool, n_obs, p_na) %>% 
  pivot_longer(cols = c(MAP, MAT, anpp, gamma_rich, HumanFootprint, NDeposition),
               values_to = "x",
               names_to = "var") %>% 
  mutate(var = case_when(var == "anpp" ~ "ANPP",
                         var == "gamma_rich" ~ "Gamma Diversity",
                         var == "HumanFootprint" ~ "Human Footprint",
                         var == "NDeposition" ~ "N Deposition",
                         .default = var) %>% 
           factor(levels=c('MAP','MAT','Gamma Diversity','ANPP','Human Footprint','N Deposition')))

envGradientPlot <- ggplot(data=df_m2, aes(y = p,
             x = x,
             color = CONTINENT)) +
  geom_point(alpha = 0.7, size=5) +
  ylab('Relative Deviation of Functional Distance') +
  facet_wrap(facets =~ var, 
             scales = "free",
             strip.position = "bottom")  +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(size=40, vjust=0.5, margin=margin(r=15)),
        axis.title.x = element_blank(),
        legend.position='bottom')

ggsave("FigH_RDFD_envGradient.png", envGradientPlot, width = 30, height = 15, dpi = 400)

################################################################################
##  06b_Q2_trait_analysis.R: link codominants to trait data
##
##  Authors: Akira Terui, Kim Komatsu
##  Date created: 5/30/2025
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
                      
                      ## select one site
                      df_codom_i <- df_codom_cl %>% 
                        filter(site_proj_comm == k) %>% 
                        mutate(pair_id = as.numeric(factor(pair_id)))
                      
                      utrt <- unique(df_codom_i$trt_type)
                      
                      df_dist <- foreach(l = utrt,
                                         .combine = bind_rows) %do% {
                                           ## subset by treatment type
                                           df_codom_il <- df_codom_i %>% 
                                             filter(trt_type == l)
                                           
                                           ## get trait distance of the co-dominants
                                           ## use `distinct()` to remove duplicates - duplicates appear when tri-dominants and co-dominants share species
                                           df_dist_l <- get_dist(data = df_codom_i,
                                                                 pool = pool,
                                                                 md = md) %>% 
                                             mutate(p_na = unique(df_pool_i$p_na),
                                                    trt_type = unique(df_codom_il$trt_type)) %>% 
                                             distinct()  
                                           
                                           return(df_dist_l)
                                         }
                      
                      return(df_dist)
                    }

df_p_trt <- df_p_trt %>% 
  mutate(trt_category = case_when(trt_type %in% c('CO2','N','P','N*P','mult_nutrient','irr','K') ~ 'Resource',
                                  trt_type %in% c('drought','temp','herb_removal') ~ 'Stress',
                                  trt_type == 'multiple_trts' ~ 'Mult. Trts',
                                  .default = 'Other'),
         trt_category = factor(trt_category,
                               levels=c('Stress', 'Resource', 'Other', 'Mult. Trts')))

# comparison, control vs. treatment ---------------------------------------

df_p_ctl <- readRDS("data/traitp_ctr.rds")

df_lm <- df_p_ctl %>% 
  select(site_proj_comm,
         ses,
         p,
         n_pool,
         n_obs) %>% 
  mutate(treatment = "control") %>% 
  bind_rows(
    select(df_p_trt,
           c(site_proj_comm,
             ses,
             p,
             n_pool,
             n_obs,
             treatment = trt_category))
  ) %>% 
  drop_na(treatment) %>% 
  filter(treatment != "CO2") %>% 
  mutate(treatment = fct_relevel(treatment,
                                 "control"))
df_lm %>%
  ggplot(aes(y = treatment,
             x = ses)) +
  ggridges::geom_density_ridges()

lm(ses ~ treatment,
   df_lm) %>% 
  summary()

# figures -----------------------------------------------------------------

df_p_trt$trt_type2 <- factor(df_p_trt$trt_type,
                             levels=c('CO2','drought','temp','herb_removal',
                                      'multiple_trts',
                                      'N','P','N*P','mult_nutrient','irr'),
                             labels=c('Control','Drought','Warming','Herbivore Rem.',
                                      'Mult. Trts',
                                      'N','P','N*P','Mult. Nutrients','Irrigation'))

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=22, vjust=-0.35, margin=margin(t=15)),
             axis.text.x=element_text(size=12),
             axis.title.y=element_text(size=22, angle=90, vjust=0.5, margin=margin(r=15)), 
             axis.text.y=element_text(size=12),
             plot.title = element_text(size=20, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             strip.text=element_text(size=18),
             legend.title=element_text(size=20), legend.text=element_text(size=20))



# figure for control

dens_ctl <- df_p_ctl %>%
  group_map(~ {
    dens <- density(.x$p, n = 500000, bw = "nrd0",
                    from = 0,
                    to = 1)  
    tibble(x = dens$x, y = dens$y)
  }) %>%
  bind_rows()

dens_ctl <- dens_ctl %>%
  mutate(width = lead(x) - x,
         width = ifelse(is.na(width), lag(width), width))  # handle last NA

png("figure_ctl_gradient.png", width = 12, height = 10, units='in', res = 300)
ggplot(df_p_trt) +
  geom_tile(data = dens_ctl,
            aes(x = x, y = y / 2, height = y, width = width, fill = x),
            inherit.aes = FALSE) +
  scale_fill_gradient(low="#FC9F32", high="#1A2766") +
  facet_wrap(~trt_type2, ncol = 5, scales = "free_y") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none') +
  labs(y = "Density",
       x = "Relative Deviation of Functional Distance")
dev.off()


# figure for treatments
dens_df <- df_p_trt %>%
  group_by(trt_type2) %>%
  filter(n() >= 2) %>%  # Remove tiny groups
  filter(!(trt_type2 %in% c('K', 'Other'))) %>%
  group_map(~ {
    dens <- density(.x$p, n = 500000, bw = "nrd0",
                    from = 0,
                    to = 1)  
    tibble(x = dens$x, y = dens$y, trt_type2 = .y$trt_type2)
  }) %>%
  bind_rows() %>%
  arrange(trt_type2, x) %>%
  group_by(trt_type2) %>%
  mutate(
    xend = lead(x),
    yend = lead(y)
  ) %>%
  drop_na()

png("figure_3_traits_bytrt.png", width = 12, height = 10, units='in', res = 300)
ggplot(df_p_trt) +
  geom_density(data = df_p_ctl,
               aes(x = p),
               fill = "lightgrey", color = NA, alpha = 0.75) +
  geom_segment(data = subset(dens_df, !(trt_type2 %in% c('K', 'other', 'NA'))), 
               aes(x = x, y = y, xend = xend, yend = yend, color = x, group = trt_type2), 
               linewidth = 3,
               lineend='round') +
  scale_color_gradient(low="#FC9F32", high="#1A2766") +
  facet_wrap(~trt_type2, ncol=5,
             scales = "free_y") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank()) +
  labs(y = "Density",
       x = "Relative Deviation of Functional Distance") +
  theme(legend.position='none')
dev.off()

### edit in photoshop to overlap control gradient panel with the grey background for all trts ###


# figure for treatment categories (stress vs resource)
dens_df <- df_p_trt %>%
  group_by(trt_category) %>%
  filter(n() >= 2) %>%  # Remove tiny groups
  # filter(!(trt_category %in% c('K', 'Other'))) %>%
  group_map(~ {
    dens <- density(.x$p, n = 500000, bw = "nrd0",
                    from = 0,
                    to = 1)  
    tibble(x = dens$x, y = dens$y, trt_category = .y$trt_category)
  }) %>%
  bind_rows() %>%
  arrange(trt_category, x) %>%
  group_by(trt_category) %>%
  mutate(
    xend = lead(x),
    yend = lead(y)
  ) %>%
  drop_na()

png("figure_3_traits_byCategory.png", width = 10, height = 5, units='in', res = 300)
ggplot(df_p_trt) +
  geom_density(data = df_p_ctl,
               aes(x = p),
               fill = "lightgrey", color = NA, alpha = 0.75) +
  geom_segment(data = subset(dens_df, !(trt_category %in% c('NA'))), 
               aes(x = x, y = y, xend = xend, yend = yend, color = x, group = trt_category), 
               linewidth=3,
               lineend='round') +
  scale_color_gradient(low="#FC9F32", high="#1A2766") +
  facet_wrap(~trt_category, ncol=5) +
  ylim(0,1.4) +
  scale_y_continuous(breaks=seq(0,1.4, by=0.2), limits = c(0, 1.4) ) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank()) +
  labs(y = "Density",
       x = "Relative Deviation of Functional Distance") +
  theme(legend.position='none')
dev.off()

#figure for ctl
dens_ctl <- df_p_ctl %>%
  group_map(~ {
    dens <- density(.x$p, n = 500000, bw = "nrd0",
                    from = 0,
                    to = 1)  
    tibble(x = dens$x, y = dens$y)
  }) %>%
  bind_rows()

dens_ctl <- dens_ctl %>%
  mutate(width = lead(x) - x,
         width = ifelse(is.na(width), lag(width), width))  # handle last NA

png("figure_ctl_gradient_category.png", width = 10, height = 5, units='in', res = 300)
ggplot(df_p_trt) +
  geom_tile(data = dens_ctl,
            aes(x = x, y = y / 2, height = y, width = width, fill = x),
            inherit.aes = FALSE) +
  scale_fill_gradient(low="#FC9F32", high="#1A2766") +
  facet_wrap(~trt_category, ncol = 5) +
  # ylim(0,1.4) +
  scale_y_continuous(breaks=seq(0,1.4, by=0.2), limits = c(0, 1.4)) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none') +
  labs(y = "Density",
       x = "Relative Deviation of Functional Distance")
dev.off()


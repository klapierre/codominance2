################################################################################
##  06b_Q2_trait_analysis.R: link codominants to trait data
##
##  Authors: Akira Terui (modifid by K. Komatsu)
##  Date created: 5/30/2025
################################################################################

# PROCEDURE:
# 1. select sites with codominants with trait data
# 2. calculate species trait distance for selected sites
#   - for codominant pairs (observed)
#   - for all species pairs (null distribution)
# 3. keep how many species/sites are excluded

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/01_library.R")
source("code/02_functions.R")

# read data ---------------------------------------------------------------

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
df_pool0 <- readRDS("data/allSppList.rds") %>% 
  mutate(site_proj_comm = str_to_lower(site_proj_comm)) %>% 
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
  group_by(site_proj_comm,
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
  drop_na(ldmc:n_fixation_type) %>% 
  select(c(site_proj_comm:species, 
           ldmc:n_fixation_type)) %>% 
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

df_ses <- foreach(k = usite,
                  .combine = bind_rows) %do% {
                    
                    ## define species pool of the site
                    df_pool_i <- df_pool_cl %>% 
                      filter(site_proj_comm == k) %>% 
                      arrange(species)
                    
                    ## trait distance matrix for the site
                    md <- df_pool_i %>%
                      select(ldmc:n_fixation_type) %>% 
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
                    
                    df_dist <- get_dist(data = df_codom_i,
                                        pool = pool,
                                        md = md) %>% 
                      left_join(distinct(df_codom_i,
                                         pair_id,
                                         trt_type,
                                         trt_cat),
                                by = "pair_id") %>% 
                      mutate(p_na = unique(df_pool_i$p_na))
                    
                    return(df_dist)
                  } %>% 
  filter(trt_cat != "other") # NOTE: category "other" has only 5 data points

# ## ses-based analysis
# m <- lm(ses ~ trt_cat,
#         data = df_ses, 
#         weights = 1 - p_na)
# summary(m)

## RDFD-based analysis, beta-binomial by treatment category
glmmTMB(cbind(n_obs, n_pool - n_obs) ~ trt_cat,
        data = df_ses,
        family = betabinomial(),
        weights = 1 - p_na) %>% 
  summary()

## RDFD-based analysis, beta-binomial by treatment type
glmmTMB(cbind(n_obs, n_pool - n_obs) ~ trt_type,
        data = filter(df_ses, trt_type!='co2'),
        family = betabinomial(),
        weights = 1 - p_na) %>% 
  summary()

# figures -----------------------------------------------------------------

## set theme
theme_set(theme_bw())
theme_update(plot.title = element_text(size = 20, 
                                       vjust = 2),
             axis.title.x = element_text(size = 18, 
                                         vjust = -0.35, 
                                         margin = margin(t = 15)),
             axis.title.y = element_text(size = 18,
                                         angle = 90,
                                         vjust = 0.5,
                                         margin = margin(r = 15)), 
             axis.text = element_text(size = 9),
             panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(),
             strip.text = element_text(size = 18),
             legend.text = element_text(size = 9))

## set new column for visual
df_ses <- df_ses %>% 
  mutate(trt_label = case_when(trt_cat == "control" ~ "Control",
                               trt_cat == "mult" ~ "Mult. Trt",
                               trt_cat == "resource" ~ "Resource",
                               trt_cat == "stress" ~ "Stress") %>% 
           factor(levels = c("Control",
                             "Stress",
                             "Resource",
                             "Mult. Trt")))

## density gradient
df_dens <- df_ses %>%
  group_by(trt_label) %>% 
  reframe(x = density(p,
                      n = 1000,
                      bw = "nrd0",
                      from = 0,
                      to = 1)$x,
          y = density(p,
                      n = 1000,
                      bw = "nrd0",
                      from = 0,
                      to = 1)$y) %>%
  group_by(trt_label) %>% 
  mutate(xend = ifelse(is.na(lead(x)),
                       x,
                       lead(x)),
         yend = ifelse(is.na(lead(y)),
                       y,
                       lead(y)),
         width = ifelse(is.na(lead(x) - x),
                        abs(lag(x) - x),
                        lead(x) - x)) %>% 
  ungroup()

## draw figure
g_p_main <- df_ses %>% 
  ggplot() +
  geom_tile(data = df_dens %>% 
              filter(trt_label == "Control") %>% 
              dplyr::select(-trt_label),
            aes(x = x,
                y = y / 2,
                height = y,
                width = width),
            fill = "lightgrey",
            alpha = 0.75,
            inherit.aes = FALSE) +
  geom_tile(data = df_dens %>% 
              filter(trt_label == "Control"),
            aes(x = x,
                y = y/2,
                height = y,
                width = width,
                fill = x),
            inherit.aes = FALSE) +
  geom_segment(data = df_dens %>% 
                 filter(trt_label != "Control"), 
               aes(x = x,
                   y = y, 
                   xend = xend,
                   yend = yend,
                   color = x, 
                   group = trt_label),
               linewidth = 2,
               lineend = 'round') +
  facet_wrap(~trt_label,
             ncol = 4) +
  labs(y = "Density",
       x = "Relative Deviation of Functional Distance",
       fill = NULL) +
  scale_fill_gradient(low="#FC9F32", 
                      high="#1A2766", 
                      breaks = c(0, 0.5, 1.0)) +
  scale_color_gradient(low="#FC9F32",  
                       high="#1A2766") +
  guides(color = "none") +
  theme(strip.background = element_blank(),
        legend.position = "none")

main_grob <- ggplotGrob(g_p_main)

legend <- get_legend(g_p_main + theme(legend.position = "bottom",
                                      plot.margin = margin(0, 0, 0, 0),
                                      legend.text = element_text(size = 12)))

# Create left and right title labels as grobs
left_text <- ggplot() + 
  theme_void() +
  annotate("text", 
           x = 0.1, 
           y = 0.5, 
           label = "Similar traits \n(habitat filtering)", 
           angle = 0, 
           size = 5)

right_text <- ggplot() + 
  theme_void() +
  annotate("text", 
           x = 0.5, 
           y = 0.5, 
           label = "Dissimilar traits \n(niche differentiation)", 
           angle = 0,
           size = 5)

# Combine horizontally: left text | colorbar | right text
legend_block <- arrangeGrob(left_text, legend, right_text,
                            ncol=3, widths=c(3,1,3))


## combine and export
g_p <- arrangeGrob(main_grob, legend_block,
                   ncol=1, heights=c(5,1))

ggsave("Fig3_all.png", g_p, width = 9, height = 5, dpi = 400)



### by treatment type ---------------------------------------------------------

df_ses <- df_ses %>% 
  mutate(trt_label = case_when(trt_type == "control" ~ "Control",
                               trt_type == "n" ~ "N",
                               trt_type == "p" ~ "P",
                               trt_type == "k" ~ "K",
                               trt_type == "n*p" ~ "NP",
                               trt_type == "mult_nutrient" ~ "Mult. Nut.",
                               trt_type == "co2" ~ "CO2",
                               trt_type == "irr" ~ "Irrigation",
                               trt_type == "drought" ~ "Drought",
                               trt_type == "temp" ~ "Warming",
                               trt_type == "herb_removal" ~ "Herb. Removal",
                               trt_type == "multiple_trts" ~ "Mult. Trts") %>% 
           factor(levels = c("Control","N","P","K","NP","Mult. Nut.",
                             "CO2","Irrigation","Drought","Warming",
                             "Herb. Removal","Mult. Trts")))
## density gradient
df_dens <- df_ses %>%
  filter(trt_label!='CO2') %>% 
  group_by(trt_label) %>% 
  reframe(x = density(p,
                      n = 1000,
                      bw = "nrd0",
                      from = 0,
                      to = 1)$x,
          y = density(p,
                      n = 1000,
                      bw = "nrd0",
                      from = 0,
                      to = 1)$y) %>%
  group_by(trt_label) %>% 
  mutate(xend = ifelse(is.na(lead(x)),
                       x,
                       lead(x)),
         yend = ifelse(is.na(lead(y)),
                       y,
                       lead(y)),
         width = ifelse(is.na(lead(x) - x),
                        abs(lag(x) - x),
                        lead(x) - x)) %>% 
  ungroup()

## draw figure
g_p_main <- df_ses %>% 
  ggplot() +
  geom_tile(data = df_dens %>% 
              filter(trt_label == "Control") %>% 
              dplyr::select(-trt_label),
            aes(x = x,
                y = y / 2,
                height = y,
                width = width),
            fill = "lightgrey",
            alpha = 0.75,
            inherit.aes = FALSE) +
  geom_tile(data = df_dens %>% 
              filter(trt_label == "Control"),
            aes(x = x,
                y = y/2,
                height = y,
                width = width,
                fill = x),
            inherit.aes = FALSE) +
  geom_segment(data = df_dens %>% 
                 filter(trt_label != "Control"), 
               aes(x = x,
                   y = y, 
                   xend = xend,
                   yend = yend,
                   color = x, 
                   group = trt_label),
               linewidth = 2,
               lineend = 'round') +
  facet_wrap(~trt_label,
             ncol = 4) +
  labs(y = "Density",
       x = "Relative Deviation of Functional Distance",
       fill = NULL) +
  scale_fill_gradient(low="#FC9F32", 
                      high="#1A2766", 
                      breaks = c(0, 0.5, 1.0)) +
  scale_color_gradient(low="#FC9F32",  
                       high="#1A2766") +
  guides(color = "none") +
  theme(strip.background = element_blank(),
        legend.position = "none")

main_grob <- ggplotGrob(g_p_main)

legend <- get_legend(g_p_main + theme(legend.position = "bottom",
                                      plot.margin = margin(0, 0, 0, 0),
                                      legend.text = element_text(size = 12)))

# Create left and right title labels as grobs
left_text <- ggplot() + 
  theme_void() +
  annotate("text", 
           x = 0.1, 
           y = 0.5, 
           label = "Similar traits \n(habitat filtering)", 
           angle = 0, 
           size = 5)

right_text <- ggplot() + 
  theme_void() +
  annotate("text", 
           x = 0.5, 
           y = 0.5, 
           label = "Dissimilar traits \n(niche differentiation)", 
           angle = 0,
           size = 5)

# Combine horizontally: left text | colorbar | right text
legend_block <- arrangeGrob(left_text, legend, right_text,
                            ncol=3, widths=c(3,1,3))


## combine and export
g_p <- arrangeGrob(main_grob, legend_block,
                   ncol=1, heights=c(5,1))

ggsave("FigG_traits_byTrt.png", g_p, width = 9, height = 9, dpi = 400)

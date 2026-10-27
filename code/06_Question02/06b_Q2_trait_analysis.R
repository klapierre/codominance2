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

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/01_library.R")
source("code/02_functions.R")

# read data ---------------------------------------------------------------

## data for trait data
df_trait <- readRDS("data/data_trait.rds")

## data for codominant species
df_codom0 <- readRDS("data/data_codom.rds") %>% 
  filter(trt_type == "control")

## data for species pool
df_pool0 <- readRDS("data/data_pool.rds")  

# select unique combo of site, project, community-type --------------------

# codominance dataframe ####
## flag species row if any traits are missing
## use "df_codom_cl" for analysis
df_codom <- df_codom0 %>% 
  distinct(species) %>% 
  left_join(df_trait) %>% 
  mutate(
    flag = if_any(everything(), .fns = is.na)
  ) %>% 
  select(species, flag, ldmc:n_fixation_type) %>% 
  right_join(df_codom0) %>% 
  arrange(site_proj_comm, pair_id) %>% 
  group_by(
    site_proj_comm,
    pair_id
  ) %>% 
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
  left_join(df_trait,
            by = "species") %>% 
  mutate(flag = if_any(everything(), .fns = is.na)) %>% 
  group_by(
    site_proj_comm,
    pair_id
  ) %>% 
  mutate(
    flag_pair = any(flag)
  ) %>% 
  group_by(
    site_proj_comm
  ) %>% 
  mutate(
    n_na = sum(flag_pair) / 2,
    p_na = n_na / n_distinct(pair_id)
  ) %>% 
  ungroup()

df_pool_cl <- df_pool %>% 
  filter(!flag_pair) %>% 
  select(
    site_proj_comm:species, 
    p_na,
    ldmc:n_fixation_type
  )

# trait distance ----------------------------------------------------------

## pool
usite <- sort(unique(df_codom_cl$site_proj_comm))

## internal function for trait distance calculation
fdfun <- function(data) {
  data %>% 
    select(ldmc:n_fixation_type) %>% 
    mutate(
      across(
        .cols = where(is.character),
        .fns = function(x) {
          if (n_distinct(x) > 2) {
            return(factor(x))
          } else {
            return(as.numeric(factor(x)))
          }
        }
      )
    ) %>% 
    as.data.frame() %>% 
    FD::gowdis()
}

## calculate RDFD
df_p <- foreach(k = usite,
                .combine = bind_rows) %do% {
                  
                  print(which(k == usite))
                  
                  ## define species pool of the site
                  df_pool_k <- df_pool_cl %>% 
                    filter(site_proj_comm == k) %>% 
                    arrange(species)
                  
                  ## get a vector of trait distance for the pool
                  idxp <- sort(unique(df_pool_k$pair_id))
                  
                  d_pool <- 
                    sapply(
                      idxp,
                      FUN = function(x) {
                        df_pool_k %>% 
                          filter(pair_id == x) %>% 
                          fdfun()
                      }
                    )
                  
                  ## select codominants for the site
                  df_codom_k <- df_codom_cl %>% 
                    filter(site_proj_comm == k)
                  
                  ## a vector of trait distance for codominants
                  idx <- sort(unique(df_codom_k$pair_id))
                  
                  df_d_codom <- 
                    lapply(
                      X = idx,
                      FUN = function(x) {
                        df_pair <- df_codom_k %>% 
                          filter(pair_id == x)
                        
                        mfd <- fdfun(df_pair) %>% 
                          as.matrix()
                        
                        get_dist(df_pair, mfd) %>% 
                          mutate(pair_id = x, .before = sp1)
                      }
                    ) %>% 
                    bind_rows()
                  
                  ## RDFD = n_obs / n_pool
                  n_obs <- sapply(
                    X = df_d_codom[["dist"]],
                    FUN = function(y) {
                      sum(y > d_pool)
                    }
                  )
                  
                  df_d_codom %>% 
                    mutate(
                      site_proj_comm = k,
                      n_obs = n_obs,
                      n_pool = length(d_pool),
                      rdfd = n_obs / n_pool,
                      p_na = unique(df_pool_k$p_na)
                    ) %>% 
                    relocate(site_proj_comm)
                }

saveRDS(df_p, file = "data/traitp_ctr.rds")

# analysis ----------------------------------------------------------------

## RDFD value and environmental factors
df_rich <- readRDS('data/modeSiteRichness.rds') %>% 
  mutate(
    site_proj_comm = paste(site_code, 
                           project_name, 
                           community_type, 
                           sep='_') %>% 
      str_to_lower()
  ) %>% 
  select(
    database, 
    site_proj_comm, 
    site_rich_mode
  )

df_m <- readRDS("data/envData.rds") %>% 
  janitor::clean_names() %>% 
  mutate(
    site_proj_comm = paste0(site_code, "_", 
                            project_name, "_",
                            community_type) %>% 
      str_to_lower(),
    .before = site_code
  ) %>% 
  dplyr::select(-c(site_code,
                   project_name,
                   community_type)) %>% 
  left_join(df_rich) %>% 
  right_join(
    df_p,
    by = "site_proj_comm"
  )

## model with environmental gradient
glmmTMB::glmmTMB(cbind(n_obs, n_pool - n_obs) ~
                   scale(map) +
                   scale(mat) +
                   # scale(anpp) +
                   scale(gamma_rich) +
                   scale(site_rich_mode) +
                   scale(human_footprint) +
                   scale(n_deposition) +
                   (1 | site_proj_comm),
                 family = glmmTMB::betabinomial,
                 data = df_m, 
                 weights = 1 - p_na) %>% 
  summary()

# figure ------------------------------------------------------------------

theme_set(theme_bw())
theme_update(
  axis.title.x = 
    element_text(
      size = 40, 
      vjust = -0.35, 
      margin = margin(t = 15)),
  axis.text.x = 
    element_text(
      size = 34, 
      color = 'black'
    ),
  axis.title.y = 
    element_text(
      size = 40, 
      angle = 90, 
      vjust = 0.5, 
      margin = margin(r = 15)), 
  axis.text.y = 
    element_text(
      size = 34, 
      color='black'
    ),
  plot.title = 
    element_text(
      size = 40, 
      vjust=2
    ),
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  legend.title = element_blank(), 
  legend.text = element_text(size = 40)
)

df_m2  <- df_m %>% 
  dplyr::select(
    database, 
    site_proj_comm,
    map:site_rich_mode, 
    pair_id, 
    sp1, 
    sp2, 
    dist, 
    # ses, 
    rdfd,
    n_pool,
    n_obs,
    p_na
  ) %>% 
  pivot_longer(
    cols = c(map,
             mat, 
             # anpp, 
             gamma_rich, 
             site_rich_mode,
             human_footprint, 
             n_deposition),
    values_to = "x",
    names_to = "var"
  ) %>% 
  mutate(
    var = case_when(var == "map" ~ "MAP",
                    var == "mat" ~ "MAT",
                    # var == "anpp" ~ "ANPP",
                    var == "gamma_rich" ~ "Gamma Diversity",
                    var == "site_rich_mode" ~ "Alpha Diversity",
                    var == "human_footprint" ~ "Human Footprint",
                    var == "n_deposition" ~ "N Deposition",
                    .default = var) %>% 
      factor(levels = c('MAP','MAT','Gamma Diversity','Alpha Diversity','Human Footprint','N Deposition'))
  )

envGradientPlot <- ggplot(data=df_m2, aes(y = rdfd,
                                          x = x)) +
  geom_point(alpha = 0.7, size=5) +
  # geom_smooth(
  #   data = subset(df_m2, var == "MAP"),
  #   aes(x = x, y = rdfd, group = 1),
  #   method = "lm",
  #   se = F,
  #   color = "black",
  #   linewidth = 1.5,
  #   inherit.aes = FALSE
  # ) +
  ylab('Relative Deviation of Functional Distance') +
  facet_wrap(facets =~ var, 
             scales = "free",
             strip.position = "bottom")  +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(size=40, vjust=0.5, margin=margin(r=15)),
        axis.title.x = element_blank(),
        legend.position='bottom')

# ggsave("FigS7_RDFD_envGradient.png", envGradientPlot, width = 30, height = 15, dpi = 400)



# set aside ---------------------------------------------------------------

# # how many plants are identified to species level
# completeSpp <- df_pool0 %>%
#   select(genus_species) %>%
#   unique() %>% #start with 7725 unique species identifiers
#   separate(col=genus_species, into=c('genus', 'species', 'other1', 'other2',
#                                      'other3', 'other4', 'other5', 'other6',
#                                      'other7', 'other10'), sep=' ', remove=F) %>%
#   filter(!(species %in% c('spp', 'SP.', 'sp', 'species', '-', '.',
#                           '1', '2', '3', 'seedling')),
#          !is.na(species),
#          !(genus %in% c('unknown', 'unk', 'ZZZZ', 'Unknown')),
#          !(genus_species %in% c('forb, unidentified', 'geranium pot',
#                                 'rosa spec', 'pycnan vir', 'annual forb',
#                                 'fleshy forb', 'forb hibiscus',
#                                 'cottony forb', 'forb six',
#                                 'grass, unidentified', 'un lily',
#                                 'unidentified dicot', 'poaceae family',
#                                 'cyperaceae family', 'ericaceae family',
#                                 'cryptantha blue', 'Stiff leathery',
#                                 'Dangle flower', 'Aster opposite',
#                                 'Long thin', '')),
#          is.na(other1)) %>%  #6211 species identified to species level
#   select(genus_species) %>%
#   rename(species=genus_species) %>%
#   left_join(df_trait) %>%
#   na.omit() #3817 species have all traits
# 
# write.csv(test, 'C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\1_first author\\codominance\\data\\species that still need traits_full.csv')



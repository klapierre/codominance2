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

## data for trait data
df_trait <- readRDS("data/data_trait.rds")

## data for codominant species
df_codom0 <- readRDS("data/data_codom.rds")

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

df_p <- foreach(k = usite,
                .combine = bind_rows) %do% {
                  
                  print(which(k == usite))
                  
                  ## define species pool of the site
                  df_pool_k <- df_pool_cl %>% 
                    filter(site_proj_comm == k) %>% 
                    arrange(pair_id)
                  
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
                        
                        df_trt <- df_pair %>% 
                          group_by(pair_id) %>% 
                          summarize(
                            trt_cat = unique(trt_cat), 
                            .groups = "drop"
                          )
                        
                        mfd <- fdfun(df_pair) %>% 
                          as.matrix()
                        
                        get_dist(df_pair, mfd) %>% 
                          mutate(
                            pair_id = x, 
                            .before = sp1
                          ) %>% 
                          left_join(df_trt)
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
                      p_na = unique(df_pool_k$p_na)
                    ) %>% 
                    relocate(site_proj_comm)
                } %>% 
  filter(trt_cat != "other") # NOTE: category "other" has only 5 data points

saveRDS(df_p, file = "data/traitp_trt.rds")

## RDFD-based analysis, beta-binomial by treatment category
glmmTMB(cbind(n_obs, n_pool - n_obs) ~ trt_cat,
        data = df_p,
        family = betabinomial(),
        weights = 1 - p_na) %>% 
  summary()

## RDFD-based analysis, beta-binomial by treatment type
glmmTMB(cbind(n_obs, n_pool - n_obs) ~ trt_type,
        data = filter(df_p, trt_type !='co2'),
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
df_p <- df_p %>% 
  mutate(p = n_obs / n_pool) %>% 
  mutate(trt_label = case_when(trt_cat == "control" ~ "Control",
                               trt_cat == "mult" ~ "Mult. Trt",
                               trt_cat == "resource" ~ "Resource",
                               trt_cat == "stress" ~ "Stress") %>% 
           factor(levels = c("Control",
                             "Stress",
                             "Resource",
                             "Mult. Trt")))

## density gradient
df_dens <- df_p %>%
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
g_p_main <- df_p %>% 
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

df_p <- df_p %>% 
  mutate(trt_label = case_when(trt_type == "control" ~ "(A) Unmanipulated",
                               trt_type == "n" ~ "(B) N",
                               trt_type == "p" ~ "(C) P",
                               trt_type == "k" ~ "(D) K",
                               trt_type == "n*p" ~ "(E) NP",
                               trt_type == "mult_nutrient" ~ "(F) Mult. Nutrients",
                               trt_type == "co2" ~ "CO2",
                               trt_type == "irr" ~ "(G) Irrigation",
                               trt_type == "drought" ~ "(H) Drought",
                               trt_type == "temp" ~ "(I) Warming",
                               trt_type == "herb_removal" ~ "(J) Herb. Removal",
                               trt_type == "multiple_trts" ~ "(K) Mult. Treatments") %>% 
           factor(levels = c("(A) Unmanipulated","(B) N","(C) P","(D) K","(E) NP","(F) Mult. Nutrients",
                             "CO2","(G) Irrigation","(H) Drought","(I) Warming",
                             "(J) Herb. Removal","(K) Mult. Treatments")))
## density gradient
df_dens <- df_p %>%
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
g_p_main <- df_p %>% 
  ggplot() +
  geom_tile(data = df_dens %>% 
              filter(trt_label == "(A) Unmanipulated") %>% 
              dplyr::select(-trt_label),
            aes(x = x,
                y = y / 2,
                height = y,
                width = width),
            fill = "lightgrey",
            alpha = 0.75,
            inherit.aes = FALSE) +
  geom_tile(data = df_dens %>% 
              filter(trt_label == "(A) Unmanipulated"),
            aes(x = x,
                y = y/2,
                height = y,
                width = width,
                fill = x),
            inherit.aes = FALSE) +
  geom_segment(data = df_dens %>% 
                 filter(trt_label != "(A) Unmanipulated"), 
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
           label = "Similar\ntraits", 
           angle = 0, 
           size = 6)

right_text <- ggplot() + 
  theme_void() +
  annotate("text", 
           x = 0.5, 
           y = 0.5, 
           label = "Dissimilar\ntraits", 
           angle = 0,
           size = 6)

# Combine horizontally: left text | colorbar | right text
legend_block <- arrangeGrob(left_text, legend, right_text,
                            ncol=3, widths=c(3,1,3))


## combine and export
g_p <- arrangeGrob(main_grob, legend_block,
                   ncol=1, heights=c(5,1))

ggsave("Fig3_traits_byTrt.png", g_p, width = 10, height = 9, dpi = 400)

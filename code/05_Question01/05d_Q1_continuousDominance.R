################################################################################
##  05d_Q1_continuousDominance.R: Calculate dominance as integers in control plots and 
##  compare to environmental data.
##
##  Authors: Kimberly Komatsu
##  Date created: Sept 15, 2026
################################################################################

source("code/01_library.R")
source("code/02_functions.R")

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))



# Read data ---------------------------------------------------------------

numCodomPlotYear <- readRDS("data/numCodomPlotYear.rds") %>% 
  separate(exp_unit, into=c('site_code', 'project_name', 'community_type', 
                            'plot_id', 'treatment', 'calendar_year'), 
           sep='::', remove=F) %>% 
  left_join(readRDS("data/expInfo.rds")) %>% 
  left_join(readRDS("data/envData.rds")) %>% 
  filter(trt_type=='control') %>% # control plots only
  transform(z_anpp = as.numeric(scale(anpp)),
            z_MAP = as.numeric(scale(MAP)),
            z_MAT = as.numeric(scale(MAT)),
            z_gammarich = as.numeric(scale(gamma_rich)),
            z_Ndep = as.numeric(scale(NDeposition)),
            z_humanfootprint = as.numeric(scale(HumanFootprint)))




# Calculate mode across years for all plots ----------------------------------------------------------

# create a dataframe of codominant group for plots with a single timepoint (because can't calculate mode of singleton)
singletonCodomPlotYear <- numCodomPlotYear %>% 
  group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>% 
  mutate(length=length(plot_id)) %>% 
  ungroup() %>% 
  filter(length==1) %>% 
  rename(plot_codom=num_group) %>% 
  dplyr::select(database, site_code, project_name, community_type, plot_id, trt_type, treatment, plot_codom)

# calculate mode across years for all plots
modePlotTrue <- numCodomPlotYear %>%  
  group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>% 
  reframe(plot_codom = DescTools::Mode(num_group)) %>% # mode function must be capital here 
  ungroup() %>% 
  filter(!is.na(plot_codom)) %>%
  group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>% 
  summarise(plot_codom = round(mean(plot_codom), digits=0), .groups='drop') %>% # calculate mean for ties
  rbind(singletonCodomPlotYear)

# for plots with singleton ties for modes, calculate mean and round to nearest integer
multipleMode <- numCodomPlotYear %>% 
  select(database, site_code, project_name, community_type, plot_id, trt_type, treatment, num_group, calendar_year) %>% 
  unique() %>% 
  full_join(modePlotTrue) %>% 
  filter(is.na(plot_codom)) %>%
  group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>% 
  summarise(plot_codom = round(mean(num_group), digits=0), .groups='drop')

# bind dataframes for averaged ties and true modes at plot level
modePlot <- rbind(modePlotTrue, multipleMode)


# Calculate mode across all control plots for each experiment ------------------

# create a dataframe of codominant group for experiments with a single plot (because can't calculate mode of singleton)
singletonCodomPlot <- modePlot %>% 
  filter(trt_type == 'control') %>% 
  group_by(database, site_code, project_name, community_type) %>% 
  mutate(length = length(community_type)) %>% 
  ungroup() %>% 
  filter(length == 1) %>% 
  rename(mode_site = plot_codom) %>% 
  dplyr::select(database, site_code, project_name, community_type, mode_site) 

# calculate mode across plots for each experiment, dropping those with ties
modeSiteTrue <- modePlot %>%
  filter(trt_type == 'control') %>%
  group_by(database, site_code, project_name, community_type) %>% # mode generated from these
  reframe(mode_site = DescTools::Mode(plot_codom)) %>%  
  ungroup() %>% 
  group_by(database, site_code, project_name, community_type) %>% 
  summarise(mode_site = round(mean(mode_site), digits=0), .groups='drop') %>% # calculate mean for ties
  filter(!is.na(mode_site)) %>% 
  rbind(singletonCodomPlot)

# for plots with ties for modes, calculate mean and round to nearest integer
multipleModeProj <- modePlot %>% 
  filter(trt_type == 'control') %>% 
  select(database, site_code, project_name, community_type, plot_id, plot_codom) %>% 
  full_join(modeSiteTrue) %>% 
  filter(is.na(mode_site)) %>%
  group_by(database, site_code, project_name, community_type) %>% 
  summarise(mode_site = round(mean(plot_codom), digits=0), .groups='drop')

# bind dataframes for average ties and true modes at site level
modeSite <- rbind(modeSiteTrue, multipleModeProj) %>% 
  left_join(readRDS("data/envData.rds")) %>% 
  mutate(dominance_group=case_when(mode_site==1~'a',
                                   mode_site==2~'b',
                                   mode_site==3~'b',
                                   TRUE~'c')) %>% 
  transform(z_anpp = as.numeric(scale(anpp)),
            z_MAP = as.numeric(scale(MAP)),
            z_MAT = as.numeric(scale(MAT)),
            z_gammarich = as.numeric(scale(gamma_rich)),
            z_Ndep = as.numeric(scale(NDeposition)),
            z_humanfootprint = as.numeric(scale(HumanFootprint)))
  



# Histogram --------------------------------------------------------------

# individual plots
ggplot(data=numCodomPlotYear, aes(x=num_codominants_continuous)) +
  geom_histogram() +
  xlab('Number of Dominant Species') + ylab('Count')

# site modes
ggplot(data=modeSite, aes(x=mode_site, fill=dominance_group)) +
  geom_histogram(binwidth=1, color='black') +
  xlab('Number of Dominant Species') + ylab('Count') +
  scale_x_continuous(breaks=seq(0,20, 1)) +
  scale_fill_manual(values=c("#007BA7", "#A63922", "#D8B573")) +
  theme(legend.position='none')
# ggsave("FigS3_histogramInteger.png", width = 8, height = 8, dpi = 400)

# Mixed-Effects Model ----------------------------------------------------

# individual plots
model <- glmmTMB(num_codominants_continuous ~ z_anpp*(z_MAP + z_MAT + z_gammarich) + z_Ndep + z_humanfootprint +
                   (1 | site_code / project_name / plot_id),
                 family = truncated_nbinom2,
                 data = numCodomPlotYear)
summary(model)

# site modes
model2 <- glmmTMB(mode_site ~ z_anpp*(z_MAP + z_MAT + z_gammarich) + z_Ndep + z_humanfootprint +
                   (1 | site_code),
                 family = truncated_nbinom2,
                 data = modeSite)
summary(model2)


# Figures -----------------------------------------------------------------

# ANPP * MAP figure

anpp_support <- transformCodom %>%
  filter(!is.na(z_anpp), !is.na(z_MAP)) %>%
  mutate(anpp_band = ggplot2::cut_number(z_anpp, n = 10)) %>%
  group_by(anpp_band) %>%
  summarise(z_anpp = median(z_anpp),
            z_MAP_low = quantile(z_MAP, 0.025),
            z_MAP_high = quantile(z_MAP, 0.975),
            .groups = "drop") %>%
  slice(c(1,5,10)) # low, middle, and high ANPP bands of the 20 bands

pred_MAP <- ggpredict(model,
                      terms = c("z_MAP [all]", paste0("z_anpp [", paste(anpp_support$z_anpp, collapse = ", "), "]"))) %>%
  as.data.frame() %>%
  mutate(z_anpp = as.numeric(as.character(group)))

pred_MAP_clipped <- pred_MAP %>%
  mutate(line_id = vapply(z_anpp,
                          function(x) which.min(abs(x - anpp_support$z_anpp)),
                          integer(1))) %>%
  left_join(anpp_support %>% transmute(line_id = row_number(), z_MAP_low, z_MAP_high), by = "line_id") %>%
  filter(x >= z_MAP_low, x <= z_MAP_high)

MAP_mean <- mean(transformCodom$MAP, na.rm = TRUE)
MAP_sd <- sd(transformCodom$MAP, na.rm = TRUE)
ANPP_mean <- mean(transformCodom$anpp, na.rm = TRUE)
ANPP_sd <- sd(transformCodom$anpp, na.rm = TRUE)

pred_MAP_raw <- pred_MAP_clipped %>%
  mutate(MAP_raw = x * MAP_sd + MAP_mean,
         ANPP_raw = z_anpp * ANPP_sd + ANPP_mean,
         ANPP_label = paste0(round(ANPP_raw, 0), " ANPP"))

ggplot(pred_MAP_raw, aes(x = MAP_raw, y = predicted, color = ANPP_label, fill = ANPP_label, group = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  labs(x = "MAP", y = "Predicted number of dominant species", color = "ANPP", fill = "ANPP")


# MAT figure

pred_MAT <- ggpredict(model, terms = "z_MAT [all]")

ggplot(pred_MAT, aes(x = x, y = predicted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  labs(x = "MAT",
       y = "Predicted number of dominant species")


# gamma div figure

pred_gammaDiv <- ggpredict(model, terms = "z_gammarich [all]")

ggplot(pred_gammaDiv, aes(x = x, y = predicted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  labs(x = "Gamma Diversity",
       y = "Predicted number of dominant species")
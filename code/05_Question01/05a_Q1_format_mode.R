
################################################################################
##  05a_Q1_site_modes.R: Calculate mode of codominance numbers in control plots and 
##  compare to environmental data.
##
##  Authors: Kimberly Komatsu, Jordan Winter
##  Date created: 
################################################################################

source("code/01_library.R")
source("code/02_functions.R")


# Read data ---------------------------------------------------------------

numCodomPlotYear <- readRDS("data/numCodomPlotYear.rds") %>% 
  separate(exp_unit, into=c('site_code', 'project_name', 'community_type', 
                            'plot_id', 'treatment', 'calendar_year'), 
           sep='::', remove=F) %>% 
  left_join(readRDS("data/expInfo.rds"))

 
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
  mutate(lumpMode = ifelse(mode_site == 3, 2, mode_site)) 

# how many sites?
siteCount <- modeSite %>%
  mutate(lat = round(Latitude, 1),
         long = round(Longitude, 1)) %>% 
  select(lat, long) %>% 
  unique()

# saveRDS(modePlot, file = "data/modePlot.rds")
# saveRDS(modeSite, file = "data/modeSite.rds")
# write.csv(modeSite, 'C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\1_first author\\codominance\\data\\modeSite_20251119.csv', row.names=F)

################################################################################
##  07a_Q3a_trt_modes.R: Calculate mode of codominance numbers in treatment plots and 
##  compare to control plots and treatment data.
##
##  Authors: K. Komatsu
##  Date created: May 2, 2025
################################################################################

source("code/01_library.R")
source("code/02_functions.R")


# Read data ---------------------------------------------------------------

#mode codominance across control plots and years at all sites
modesite <- readRDS("data/modeSite.rds") 


#codominance number across all plots and years
numCodomPlotYear <- readRDS("data/numCodomPlotYear.rds") %>% 
  separate(exp_unit, into=c('site_code', 'project_name', 'community_type', 'plot_id', 'treatment', 'calendar_year'), sep='::', remove=F) %>% 
  left_join(readRDS("data/expInfo.rds")) %>% 
  #create trt_type2 variable that groups some of the most common treatments together
  mutate(disturb=ifelse(trt_type %in% c("mow_clip","burn","burn*graze","disturbance","burn*mow_clip"), 1, 0), #unify codes across datasets
         # tCO2=ifelse(trt_type %in% c("CO2"), 1, 0),
         drought=ifelse(trt_type %in% c("drought"), 1, 0),
         # therb_removal=ifelse(trt_type %in% c("herb_removal"), 1, 0),
         irg=ifelse(trt_type %in% c("irr"), 1, 0),
         # ttemp=ifelse(trt_type %in% c("temp"), 1, 0),
         # tn=ifelse(trt_type %in% c("N"), 1, 0),
         # tp=ifelse(trt_type %in% c("P"), 1, 0),
         multtrts=ifelse(trt_type %in% c("CO2*temp", "burn*graze","burn*mow_clip","drought*CO2*temp","drought*mow_clip",
                                         "drought*temp*mow_clip","herb_removal*mow_clip","irr*CO2","irr*CO2*temp","irr*mow_clip",
                                         "irr*herb_removal","irr*temp*mow_clip","N*CO2*temp","N*irr*CO2","N*irr*mow_clip",
                                         "N*P*burn*graze", "mult_nutrient*irr","N*irr*CO2*temp", "N*CO2","N*mow_clip","N*burn",
                                         "N*burn*graze","N*disturbance","P*burn*graze","P*burn*mow_clip","N*drought",
                                         "N*herb_removal","P*herb_removal","N*irr","N*irr*temp","N*temp","mult_nutrient*temp",
                                         "N*P*temp","mult_nutrient*mow_clip","N*burn*mow_clip","N*P*burn","N*P*mow_clip",
                                         "P*burn","P*mow_clip","mult_nutrient*herb_removal","mult_nutrient*herb_removal*mow_clip",
                                         "temp*mow_clip","drought*temp","irr*temp","mult_nutrient","N*P"),1,0)) %>%
  mutate(trt_type2=ifelse(disturb==1, 'disturbance', ifelse(multtrts==1, 'multiple trts', trt_type)))
 
#subset to final experiment year for each experiment and only treatment plots
numCodomPlot <- numCodomPlotYear %>% 
  group_by(site_code, project_name, community_type, plot_id) %>% 
  mutate(max_year=max(as.numeric(calendar_year)),
         experiment_length=max_year-as.numeric(calendar_year)+1) %>% 
  ungroup() %>% 
  group_by(site_code, project_name, community_type, plot_id) %>% 
  filter(calendar_year==max(calendar_year)) %>% 
  ungroup() %>% 
  filter(trt_type2!='control')

# nutnetCtlOnly <- numCodomPlot %>% 
#   filter(database=='nutnet') %>% 
#   select(site_code, project_name, community_type, treatment) %>% 
#   unique() %>% 
#   group_by(site_code, project_name, community_type) %>% 
#   summarise(trt_num=length(treatment), .groups='drop') %>% 
#   filter(trt_num==1)

 
# # Calculate mode across years for all plots ----------------------------------------------------------
#  
# # create a dataframe of codominant group for plots with a single timepoint (because can't calculate mode of singleton)
# singletonCodomPlotYear <- numCodomPlotYear %>% 
#   group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>% 
#   mutate(length=length(plot_id)) %>% 
#   ungroup() %>% 
#   filter(length==1) %>% 
#   rename(plot_codom=num_group) %>% 
#   dplyr::select(database, site_code, project_name, community_type, plot_id, trt_type, treatment, plot_codom) 
#  
# # calculate mode across years for all plots
# modePlotTrue <- numCodomPlotYear %>%  
#   group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>% 
#   reframe(plot_codom = DescTools::Mode(num_group)) %>% # mode function must be capital here 
#   ungroup() %>% 
#   filter(!is.na(plot_codom)) %>%
#   group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>% 
#   summarise(plot_codom=round(mean(plot_codom), digits=0), .groups='drop') %>% # calculate mean for ties
#   rbind(singletonCodomPlotYear)
# 
# # for plots with singleton ties for modes, calculate mean and round to nearest integer
# multipleMode <- numCodomPlotYear %>% 
#   select(database, site_code, project_name, community_type, plot_id, trt_type, treatment, num_group, calendar_year) %>% 
#   unique() %>% 
#   full_join(modePlotTrue) %>% 
#   filter(is.na(plot_codom)) %>%
#   group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>% 
#   summarise(plot_codom=round(mean(num_group), digits=0), .groups='drop')
#  
# # bind dataframes for averaged ties and true modes at plot level
# modePlot <- rbind(modePlotTrue, multipleMode)


# Calculate mode across all treatment plots for each experiment ------------------

# create a dataframe of codominant group for experiments with a single plot (because can't calculate mode of singleton)
singletonCodomPlot <- numCodomPlot %>% 
  group_by(database, site_code, project_name, community_type, trt_type2) %>% 
  mutate(length=length(community_type)) %>% 
  ungroup() %>% 
  filter(length==1) %>% 
  rename(mode_trt=num_codominants) %>% 
  dplyr::select(database, site_code, project_name, community_type, trt_type2, mode_trt) 

# calculate mode across plots for each experiment, dropping those with ties
modeTrtTrue <- numCodomPlot %>%
   group_by(database, site_code, project_name, community_type, trt_type2) %>% # mode generated from these
   reframe(mode_trt = DescTools::Mode(num_codominants)) %>%  
   ungroup() %>% 
   group_by(database, site_code, project_name, community_type, trt_type2) %>% 
   summarise(mode_trt=round(mean(mode_trt), digits=0), .groups='drop') %>% # calculate mean for ties
   filter(!is.na(mode_trt)) %>% 
   rbind(singletonCodomPlot)

# for plots with ties for modes, calculate mean and round to nearest integer
multipleModeProj <- numCodomPlot %>% 
  select(database, site_code, project_name, community_type, plot_id, trt_type2, num_codominants) %>% 
  full_join(modeTrtTrue) %>% 
  filter(is.na(mode_trt)) %>%
  group_by(database, site_code, project_name, community_type, trt_type2) %>% 
  summarise(mode_trt=round(mean(num_codominants), digits=0), .groups='drop')

# bind dataframes for average ties and true modes at site level
modeTrt <- rbind(modeTrtTrue, multipleModeProj) %>% 
  left_join(readRDS("data/envData.rds")) %>% 
  mutate(lumpMode = ifelse(mode_trt == 3, 2, mode_trt)) 



# Histogram- count of codom level per variable ----------------------------

df_hist <- modeTrt %>% 
  select(mode_trt, lumpMode, MAP, MAT, gamma_rich, anpp, HumanDisturbance, N_Deposition) %>% 
  pivot_longer(cols = c("MAP", "MAT", "gamma_rich", "anpp", "HumanDisturbance", "N_Deposition"), 
               names_to = "variable", values_to = "value") %>% 
  mutate(lumpMode = ifelse(lumpMode == 1, "Monodominated", 
                    ifelse(lumpMode == 2, "Codominated", "Even")),
         mode.levels = factor(lumpMode, levels = c("Monodominated", "Codominated", "Even")))


# Visualizing distribution of codominance across sites ----------------------------------------------------------

mapData <- modeTrt %>% 
  mutate(codom_category = ifelse(mode_site == 1, "Monodominated", 
                          ifelse(mode_site == 2, "Codominated", 
                          ifelse(mode_site == 3, "Tridominated", 
                                 "Even")))) %>% 
  filter(!is.na(N_Deposition))

mapData$codom_category <- factor(mapData$codom_category, levels = c("Monodominated", "Codominated", "Tridominated", "Even"))

autumnalPalette <- c("#02385A", "#A63922", "#D8B573", 'grey')

mapData %>% 
  ggplot(aes(x="", y=mode_site, fill = codom_category)) +
  geom_bar(stat = "identity", width=1) +
  coord_polar("y",start=0)

ggplot(mapData, aes(x = codom_category, fill = codom_category)) +
  geom_histogram(stat = "count") +
  stat_count(binwidth = 1, 
             geom = 'text', 
             color = 'white', 
             aes(label = after_stat(count)),
             position = position_stack(vjust = 0.5))+
  xlab("Site Dominance")+
  ylab("Count")+
  scale_fill_manual(values = autumnalPalette) +
  theme(legend.position = "none")

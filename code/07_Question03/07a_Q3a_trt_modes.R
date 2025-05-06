################################################################################
##  07a_Q3a_trt_modes.R: Calculate mode of codominance numbers in treatment plots and 
##  compare to control plots and treatment data.
##
##  Authors: K. Komatsu
##  Date created: May 2, 2025
################################################################################

source("code/01_library.R")
source("code/02_functions.R")

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_text(size=20), legend.text=element_text(size=20))


# Read data ---------------------------------------------------------------

#mode codominance across control plots and years at all sites
modeSite <- readRDS("data/modeSite.rds") %>% 
  mutate(lump_mode_site=ifelse(mode_site %in% c(2,3), 2, mode_site))


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
                                         "temp*mow_clip","drought*temp","irr*temp","C*stone","irr*plant_mani", 
                                         'irr*plant_mani*herb_removal',"mult_nutrient*drought","mult_nutrient*fungicide",
                                         "mult_nutrient*plant_mani","mult_nutrient*plant_mani*herb_removal","N*fungicide",
                                         "N*P*burn*mow_clip","N*P*plant_mani","N*plant_mani","N*plant_mani*disturbance",
                                         "N*plant_mani*mow_clip","N*stone","N*temp*fungicide","P*plant_mani","plant_mani*disturbance",
                                         "plant_mani*herb_removal","plant_mani*mow_clip","precip_vari*temp","temp*fungicide"),1,0)) %>%
  mutate(trt_type3=ifelse(disturb==1, 'disturbance', ifelse(multtrts==1, 'multiple_trts', trt_type))) %>%  
  mutate(trt_type2=ifelse(trt_type3 %in% c('N','P','K','N*P','mult_nutrient','herb_removal','disturbance','CO2',
                                           'irr','drought','temp','multiple_trts', 'control'), trt_type3, 'other')) %>%
  filter(site_code!='ufrec.us', #filter out this site, which has no control plots
         !grepl("plant_mani", trt_type)) #filter out any treatment that directly manipulates plant species 
         
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
  rename(mode_trt=num_group) %>% 
  dplyr::select(database, site_code, project_name, community_type, trt_type2, mode_trt) 

# calculate mode across plots for each experiment, dropping those with ties
modeTrtTrue <- numCodomPlot %>%
   group_by(database, site_code, project_name, community_type, trt_type2) %>% # mode generated from these
   reframe(mode_trt = DescTools::Mode(num_group)) %>%  
   ungroup() %>% 
   group_by(database, site_code, project_name, community_type, trt_type2) %>% 
   summarise(mode_trt=round(mean(mode_trt), digits=0), .groups='drop') %>% # calculate mean for ties
   filter(!is.na(mode_trt)) %>% 
   rbind(singletonCodomPlot)

# for plots with ties for modes, calculate mean and round to nearest integer
multipleModeProj <- numCodomPlot %>% 
  select(database, site_code, project_name, community_type, plot_id, trt_type2, num_group) %>% 
  full_join(modeTrtTrue) %>% 
  filter(is.na(mode_trt)) %>%
  group_by(database, site_code, project_name, community_type, trt_type2) %>% 
  summarise(mode_trt=round(mean(num_group), digits=0), .groups='drop')

# bind dataframes for average ties and true modes at site level
modeTrt <- rbind(modeTrtTrue, multipleModeProj) %>% 
  left_join(readRDS("data/envData.rds")) %>% 
  mutate(lump_mode_trt = ifelse(mode_trt == 3, 2, mode_trt)) %>% 
  full_join(modeSite) %>% #join site codominance data (ctl plot data)
  mutate(lump_mode_trt_cat = ifelse(lump_mode_trt==1, 'Monodominated',
                             ifelse(lump_mode_trt==2, 'Codominated',
                                    'Even')),
         lump_mode_site_cat = ifelse(lump_mode_site==1, 'Monodominated',
                             ifelse(lump_mode_site==2, 'Codominated',
                                    'Even')))

modeTrt$lump_mode_trt_cat <- factor(modeTrt$lump_mode_trt_cat, levels = c("Monodominated", "Codominated", "Even"))
modeTrt$lump_mode_site_cat <- factor(modeTrt$lump_mode_site_cat, levels = c("Monodominated", "Codominated", "Even"))
modeTrt$trt_type2 <- factor(modeTrt$trt_type2, levels = c('N','P','K','N*P','mult_nutrient',
                                                          'herb_removal','disturbance',
                                                          'CO2','irr','drought','temp','other',
                                                          'multiple_trts'))

# log linear model --------------------------------------------------------------------------------
summaryTableTrt <- xtabs(~ trt_type2 + lump_mode_site_cat + lump_mode_trt_cat, data = modeTrt)

m1 <- loglm(~trt_type2+lump_mode_site_cat+lump_mode_trt_cat, data=summaryTableTrt) #independence
m2 <- loglm(~lump_mode_trt_cat*(lump_mode_site_cat+trt_type2), data=summaryTableTrt) #conditional independence
m3 <- loglm(~trt_type2+lump_mode_site_cat*lump_mode_trt_cat, data=summaryTableTrt) #joint independence (trt_type2)
m4 <- loglm(~lump_mode_site_cat+trt_type2*lump_mode_trt_cat, data=summaryTableTrt) #joint independence (lump_mode_site_cat)

anova(m1,m2,m3,m4)
#model 4 (conditional independence - trt codom number depends on both trt type and site codom)


# effect of site codominance
modeSiteCodom <- xtabs(~ lump_mode_site_cat + lump_mode_trt_cat, data = modeTrt)

print(chisq <- chisq.test(modeSiteCodom))
# X-squared = 197.43, df = 4, p-value < 2.2e-16

mosaicplot(modeSiteCodom, shade = TRUE, las=2,
           main = "modeSiteCodom")


# effect of treatment
modeTrtCodom <- xtabs(~ lump_mode_trt_cat + trt_type2, data = modeTrt)

print(chisq <- chisq.test(modeTrtCodom))
# X-squared = 53.607, df = 24, p-value = 0.0004808

mosaicplot(modeTrtCodom, shade = TRUE, las=2,
           main = "modeTrtCodom")


# Histogram- count of codom level per treatment and site codom number ----------------------------

autumnalPalette <- c("#02385A", "#A63922", "#D8B573", 'grey')

#overall
ggplot(modeTrt, aes(x = lump_mode_site_cat, fill = lump_mode_trt_cat)) +
  geom_bar(stat = "count", position='stack') +
  stat_count(geom = 'text', 
             color = 'white', 
             aes(label = after_stat(count)),
             position = position_stack(vjust = 0.5)) +
  xlab("Site Dominance") +
  ylab("Count") +
  scale_fill_manual(values = autumnalPalette, name='Treatment Dominance')

# ggsave(file='Fig4a_trtHistograms_overall.png', width=10, height=10, units='in', dpi=300, bg='white')


#stacked by treatment (don't run factor of trts above to get them all)
ggplot(modeTrt, aes(x = lump_mode_site_cat, fill = lump_mode_trt_cat)) +
  geom_bar(stat = "count", position='stack') +
  stat_count(geom = 'text', 
             color = 'white', 
             aes(label = after_stat(count)),
             position = position_stack(vjust = 0.5)) +
  xlab("Site Dominance") +
  ylab("Count") +
  scale_fill_manual(values = autumnalPalette, name='Treatment Dominance') +
  facet_grid(rows='trt_type2', scales='free_y')

# ggsave(file='Fig4b_trtHistograms_byTrt_all.png', width=10, height=30, units='in', dpi=300, bg='white')

#stacked, only subset with higher replication across sites
ggplot(subset(modeTrt, !(trt_type2 %in% c('C','fungicide','light','lime','plant_mani','precip_vari','stone') & !(is.na(trt_type2)))), 
       aes(x = lump_mode_site_cat, fill = lump_mode_trt_cat)) +
  geom_bar(stat = "count", position='stack') +
  stat_count(geom = 'text', 
             color = 'white', 
             aes(label = after_stat(count)),
             position = position_stack(vjust = 0.5)) +
  xlab("Site Dominance") +
  ylab("Count") +
  scale_fill_manual(values = autumnalPalette, name='Treatment Dominance') +
  facet_grid(rows='trt_type2', scales='free_y')

# ggsave(file='Fig4c_trtHistograms_byTrt_subset.png', width=10, height=30, units='in', dpi=300, bg='white')


#faceted, only subset with higher replication across sites
ggplot(modeTrt, aes(x = lump_mode_trt_cat, fill = lump_mode_trt_cat)) +
  geom_bar(stat = "count", position='identity') +
  # stat_count(geom = 'text', 
  #            color = 'white', 
  #            aes(label = after_stat(count)),
  #            position = position_stack(vjust = 0.5)) +
  xlab("Treatment Dominance") +
  ylab("Count") +
  scale_fill_manual(values = autumnalPalette, name='Treatment Dominance') +
  theme(legend.position='none') +
  facet_grid(rows=vars(trt_type2), cols=vars(lump_mode_site_cat), scales='free_y') 

# ggsave(file='Fig4d_trtHistograms_byTrt_facet.png', width=30, height=30, units='in', dpi=300, bg='white')
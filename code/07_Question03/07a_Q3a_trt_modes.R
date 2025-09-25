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
             legend.title=element_text(size=16), legend.text=element_text(size=14))


# Read data ---------------------------------------------------------------

#mode codominance across control plots and years at all sites
modeSite <- readRDS("data/modeSite.rds") %>% 
  mutate(lump_mode_site=ifelse(mode_site %in% c(2,3), 2, mode_site))


#codominance number across all plots and years
numCodomPlotYear <- readRDS("data/numCodomPlotYear.rds") %>% 
  separate(exp_unit, into=c('site_code', 'project_name', 'community_type', 'plot_id', 'treatment', 'calendar_year'), sep='::', remove=F) %>% 
  left_join(readRDS("data/expInfo.rds")) %>% 
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
  filter(trt_type!='control')

# nutnetCtlOnly <- numCodomPlot %>% 
#   filter(database=='nutnet') %>% 
#   select(site_code, project_name, community_type, treatment) %>% 
#   unique() %>% 
#   group_by(site_code, project_name, community_type) %>% 
#   summarise(trt_num=length(treatment), .groups='drop') %>% 
#   filter(trt_num==1)
## 60 NutNet sites are control plots only, therefore not included in analyses below
 
# Calculate mode across treatments for all plots ----------------------------------------------------------

# create a dataframe of codominant group for plots with a single treatment (because can't calculate mode of singleton)
singletonCodomPlotYear <- numCodomPlot %>%
  group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>%
  mutate(length=length(plot_id)) %>%
  ungroup() %>%
  filter(length==1) %>%
  rename(plot_codom=num_group) %>%
  dplyr::select(database, site_code, project_name, community_type, plot_id, trt_type, treatment, plot_codom)

# calculate mode across years for all plots
modePlotTrue <- numCodomPlot %>%
  group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>%
  reframe(plot_codom = DescTools::Mode(num_group)) %>% # mode function must be capital here
  ungroup() %>%
  filter(!is.na(plot_codom)) %>%
  group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>%
  summarise(plot_codom=round(mean(plot_codom), digits=0), .groups='drop') %>% # calculate mean for ties
  rbind(singletonCodomPlotYear)

# for plots with singleton ties for modes, calculate mean and round to nearest integer
multipleMode <- numCodomPlot %>%
  select(database, site_code, project_name, community_type, plot_id, trt_type, treatment, num_group, calendar_year) %>%
  unique() %>%
  full_join(modePlotTrue) %>%
  filter(is.na(plot_codom)) %>%
  group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>%
  summarise(plot_codom=round(mean(num_group), digits=0), .groups='drop')

# bind dataframes for averaged ties and true modes at plot level
modePlot <- rbind(modePlotTrue, multipleMode)


# Calculate mode across all trt types for each experiment ------------------

# create a dataframe of codominant group for experiments with a single plot (because can't calculate mode of singleton)
singletonCodomPlot <- modePlot %>% 
  group_by(database, site_code, project_name, community_type, trt_type) %>% 
  mutate(length=length(community_type)) %>% 
  ungroup() %>% 
  filter(length==1) %>% 
  rename(mode_trt=plot_codom) %>% 
  dplyr::select(database, site_code, project_name, community_type, trt_type, mode_trt) 

# calculate mode across plots for each experiment, dropping those with ties
modeTrtTrue <- modePlot %>%
   group_by(database, site_code, project_name, community_type, trt_type) %>% # mode generated from these
   reframe(mode_trt = DescTools::Mode(plot_codom)) %>%  
   ungroup() %>% 
   group_by(database, site_code, project_name, community_type, trt_type) %>% 
   summarise(mode_trt=round(mean(mode_trt), digits=0), .groups='drop') %>% # calculate mean for ties
   filter(!is.na(mode_trt)) %>% 
   rbind(singletonCodomPlot)

# for plots with ties for modes, calculate mean and round to nearest integer
multipleModeProj <- modePlot %>% 
  select(database, site_code, project_name, community_type, plot_id, trt_type, plot_codom) %>% 
  full_join(modeTrtTrue) %>% 
  filter(is.na(mode_trt)) %>%
  group_by(database, site_code, project_name, community_type, trt_type) %>% 
  summarise(mode_trt=round(mean(plot_codom), digits=0), .groups='drop')

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
  # %>% mutate(drop=ifelse(database=='nutnet' & trt_type=='herb_removal', 1, 
  #             ifelse(database=='nutnet' & trt_type=='multiple_trts', 1, 0))) %>% 
  # filter(drop==0) #drops herbivore removal trts from all NutNet sites (which often have low grazing pressure compared to GEx)

modeTrt$lump_mode_trt_cat <- factor(modeTrt$lump_mode_trt_cat, levels = c("Monodominated", "Codominated", "Even"))
modeTrt$lump_mode_site_cat <- factor(modeTrt$lump_mode_site_cat, levels = c("Monodominated", "Codominated", "Even"))
modeTrt$trt_type <- factor(modeTrt$trt_type, levels = c('N','P','K','N*P','mult_nutrient',
                                                          'herb_removal','disturbance',
                                                          'CO2','irr','drought','temp','other',
                                                          'multiple_trts'))

# saveRDS(modeTrt, file = "data/modeTrt.rds")


# trt vs ctl chisquared model --------------------------------------------------------------------------------

modeSite2 <- modeSite %>% 
  mutate(lump_mode_cat=ifelse(lump_mode_site==1, 'Monodominated',
                       ifelse(lump_mode_site==2, 'Codominated', 'Even'))) %>%
  select(database, site_code, project_name, community_type, lump_mode_cat) %>%  
  mutate(trt_type='control')

modeAll <- modeTrt %>% 
  select(database, site_code, project_name, community_type, trt_type, lump_mode_trt_cat) %>% 
  rename(lump_mode_cat=lump_mode_trt_cat) %>% 
  rbind(modeSite2) %>% 
  mutate(trt_ctl=ifelse(trt_type=='control', 'Control', 'Treatment'))

summaryTableAll <- xtabs(~ trt_ctl + lump_mode_cat, data = modeAll)

print(chisq <- chisq.test(summaryTableAll))
# X-squared = 10.195, df = 2, p-value = 0.006112

mosaicplot(summaryTableAll, shade = TRUE, las=2,
           main = "summaryTableAll")

summaryTableAllPercent <- as.data.frame(summaryTableAll) %>% 
  group_by(trt_ctl) %>%
  mutate(percent=Freq/sum(Freq)) %>% 
  ungroup()

ggplot(summaryTableAllPercent, aes(x=trt_ctl, y=lump_mode_cat)) +
  geom_tile(aes(fill=100*percent), color='white') +
  geom_text(aes(label=Freq), size=6, color='grey40') +
  # geom_text(aes(label=round(100*percent, digits=0))) +
  scale_fill_gradient(low='#F8FBF8', high='#031B88') +
  scale_y_discrete(limits=rev) +
  xlab('') + ylab('') + labs(fill='Column\nPercentage')

# ggsave(file='Fig4b_heatMapTrt.png', width=7, height=3, units='in', dpi=300, bg='white')
  

# log linear model --------------------------------------------------------------------------------
summaryTableTrt <- xtabs(~ trt_type + lump_mode_site_cat + lump_mode_trt_cat, data = modeTrt)

m1 <- loglm(~trt_type+lump_mode_site_cat+lump_mode_trt_cat, data=summaryTableTrt) #independence
m2 <- loglm(~lump_mode_trt_cat*(lump_mode_site_cat+trt_type), data=summaryTableTrt) #conditional independence
m3 <- loglm(~trt_type+lump_mode_site_cat*lump_mode_trt_cat, data=summaryTableTrt) #joint independence (trt_type)
m4 <- loglm(~lump_mode_site_cat+trt_type*lump_mode_trt_cat, data=summaryTableTrt) #joint independence (lump_mode_site_cat)
m5 <- loglm(~lump_mode_trt_cat*lump_mode_site_cat, data=summaryTableTrt) #only lump_mode_site_cat
m6 <- loglm(~lump_mode_trt_cat*trt_type, data=summaryTableTrt) #only trt_type

anova(m1,m2,m3,m4,m5,m6)
#model 2 (conditional independence - trt codom number depends on both trt type and site codom)


# effect of site codominance
modeSiteCodom <- xtabs(~ lump_mode_site_cat + lump_mode_trt_cat, data = modeTrt)

print(chisq <- chisq.test(modeSiteCodom))
# X-squared = 197.14, df = 4, p-value < 2.2e-16

mosaicplot(modeSiteCodom, shade = TRUE, las=2,
           main = "modeSiteCodom")

summaryModeSiteCodom <- as.data.frame(modeSiteCodom) %>% 
  group_by(lump_mode_site_cat ) %>%
  mutate(percent=Freq/sum(Freq)) %>% 
  ungroup()

ggplot(summaryModeSiteCodom, aes(x=lump_mode_site_cat , y=lump_mode_trt_cat)) +
  geom_tile(aes(fill=percent)) +
  geom_text(aes(label=Freq), size=6, color='grey40') +
  # geom_text(aes(label=round(100*percent, digits=0))) +
  scale_fill_gradient(low='#F8FBF8', high='#031B88') +
  scale_y_discrete(limits=rev) +
  xlab('Control') + ylab('Treatment') + labs(fill='Column\nPercentage')

# ggsave(file='Fig5a_heatMapSiteTrt.png', width=10, height=5, units='in', dpi=300, bg='white')


# effect of treatment
modeTrtCodom <- xtabs(~ lump_mode_trt_cat + trt_type, data = modeTrt)

print(chisq <- chisq.test(modeTrtCodom))
# X-squared = 57.796, df = 24, p-value = 0.0001299

mosaicplot(modeTrtCodom, shade = TRUE, las=2,
           main = "modeTrtCodom")

summaryModeTrtCodom <- as.data.frame(modeTrtCodom) %>% 
  group_by(trt_type) %>%
  mutate(percent=Freq/sum(Freq)) %>% 
  ungroup() %>% 
  mutate(trt_type_nice=ifelse(trt_type=='mult_nutrient', 'Mult. Nutrients', 
                       ifelse(trt_type=='herb_removal', 'Herbivore Rem.',
                       ifelse(trt_type=='disturbance', 'Disturbance',
                       ifelse(trt_type=='irr', 'Irrigation',
                       ifelse(trt_type=='drought', 'Drought',
                       ifelse(trt_type=='temp', 'Warming',
                       ifelse(trt_type=='other', 'Other',
                       ifelse(trt_type=='multiple_trts', 'Mult. Trts', as.character(trt_type))))))))))

summaryModeTrtCodom$trt_type_nice <- factor(summaryModeTrtCodom$trt_type_nice, 
                                            levels = c('N','P','K','N*P','Mult. Nutrients',
                                                       'Herbivore Rem.','Disturbance',
                                                       'CO2','Irrigation','Drought','Warming','Other',
                                                       'Mult. Trts'))

ggplot(summaryModeTrtCodom, aes(x=trt_type_nice , y=lump_mode_trt_cat)) +
  geom_tile(aes(fill=percent)) +
  geom_text(aes(label=Freq), size=6, color='grey40') +
  # geom_text(aes(label=round(100*percent, digits=0))) +
  scale_fill_gradient(low='#F8FBF8', high='#031B88') +
  scale_y_discrete(limits=rev) +
  xlab('') + ylab('Treatment') + labs(fill='Column\nPercentage') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# ggsave(file='Fig5b_heatMapTrtTrt.png', width=10, height=4, units='in', dpi=300, bg='white')


# Histogram- count of codom level per treatment and site codom number ----------------------------

autumnalPalette <- c("#007BA7", "#A63922", "#D8B573", 'grey')

#overall by trt codom category
ggplot(subset(modeTrt, !is.na(lump_mode_trt_cat) & !is.na(lump_mode_site_cat)), aes(x = lump_mode_trt_cat, fill = lump_mode_site_cat)) +
  geom_bar(stat = "count", position='stack') +
  stat_count(geom = 'text',
             color = 'white',
             aes(label = after_stat(count)),
             position = position_stack(vjust = 0.5)) +
  xlab("Treatment Dominance") +
  ylab("Count") +
  scale_fill_manual(values = autumnalPalette, name='Control Dominance') +
  theme(legend.position=c(0.78,0.85))

# ggsave(file='Fig4a_trtHistograms_overall.png', width=6, height=6, units='in', dpi=300, bg='white')


# Sankey diagram with trt

longModeTrt <- modeTrt %>%
  select(lump_mode_trt_cat, lump_mode_site_cat) %>%
  na.omit() %>%
  count(lump_mode_trt_cat, lump_mode_site_cat) %>%
  group_by(lump_mode_trt_cat, lump_mode_site_cat) %>%
  mutate(group = cur_group_id()) %>%
  ungroup() %>%
  pivot_longer(cols = c(lump_mode_trt_cat, lump_mode_site_cat),
               names_to = "ctl_trt",
               values_to = "category") %>%
  mutate(ctl_trt = factor(ctl_trt,
                          levels = c("lump_mode_site_cat", "lump_mode_trt_cat"),
                          labels = c("Ambient", "Treatment")))

ggplot(longModeTrt, aes(x = ctl_trt, stratum = category, alluvium = group, y = n, fill = category)) +
  geom_flow(stat = "alluvium", lode.guidance = "forward", alpha = 0.7) +
  geom_stratum(width = 1/3) +
  scale_fill_manual(values=autumnalPalette) +
  labs(x=NULL, y='Count', fill=NULL) +
  coord_cartesian(xlim=c(1.3, 1.7))

# ggsave(file='Fig4_trtSankey_overall.png', width=8, height=6, units='in', dpi=300, bg='white')



# Group into treatments categories for Sankey

longModeTrtCategory <- modeTrt  %>% 
  mutate(trt_category=ifelse(trt_type %in% c('CO2','N','P','N*P','mult_nutrient','irr','K'), 'Resource',
                      ifelse(trt_type %in% c('drought','temp','herb_removal'), 'Stress',
                      ifelse(trt_type=='multiple_trts', 'Mult. Trts', 'Other')))) %>%
  select(lump_mode_trt_cat, lump_mode_site_cat, trt_category) %>%
  na.omit() %>%
  count(lump_mode_trt_cat, lump_mode_site_cat, trt_category) %>%
  group_by(lump_mode_trt_cat, lump_mode_site_cat, trt_category) %>%
  mutate(group = cur_group_id()) %>%
  ungroup() %>%
  pivot_longer(cols = c(lump_mode_trt_cat, trt_category, lump_mode_site_cat),
               names_to = "ctl_trt",
               values_to = "category") %>%
  mutate(ctl_trt = factor(ctl_trt,
                          levels = c("lump_mode_site_cat", "trt_category", "lump_mode_trt_cat"),
                          labels = c("Ambient", "Treatment Category", "Treatment")))

# ggplot(longModeTrtCategory, aes(x = ctl_trt, stratum = category, alluvium = group, y = n, fill = category)) +
#   geom_flow(stat = "alluvium", lode.guidance = "forward", alpha = 0.7) +
#   geom_stratum(width = 1/3) 
#   # scale_fill_manual(values=autumnalPalette) +
#   # labs(x=NULL, y='Count', fill=NULL) +
#   # coord_cartesian(xlim=c(1.3, 1.7))

# ggsave(file='FigY_trtSankey_overall_category.png', width=8, height=6, units='in', dpi=300, bg='white')



# Pie charts by treatment category
  
stressWide <- modeTrt  %>% 
  mutate(trt_category=ifelse(trt_type %in% c('CO2','N','P','N*P','mult_nutrient','irr','K'), 'Resource',
                      ifelse(trt_type %in% c('drought','temp','herb_removal'), 'Stress',
                      ifelse(trt_type=='multiple_trts', 'Mult. Trts', 'Other')))) %>%
  filter(trt_category=='Stress') %>% 
  select(lump_mode_trt_cat) %>%
  na.omit() %>%
  count(lump_mode_trt_cat) %>% 
  mutate(proportion = round((n/sum(n)), digits=3)) %>% 
  arrange(proportion) %>%
  mutate(labels=scales::percent(proportion))

stressFig <- ggplot(stressWide, aes(x="", y=proportion, fill=lump_mode_trt_cat)) +
  geom_col() +
  coord_polar(theta="y") +
  scale_fill_manual(values=autumnalPalette)  +
  ggtitle('Stress') +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')


resourceWide <- modeTrt  %>% 
  mutate(trt_category=ifelse(trt_type %in% c('CO2','N','P','N*P','mult_nutrient','irr','K'), 'Resource',
                             ifelse(trt_type %in% c('drought','temp','herb_removal'), 'Stress',
                                    ifelse(trt_type=='multiple_trts', 'Mult. Trts', 'Other')))) %>%
  filter(trt_category=='Resource') %>% 
  select(lump_mode_trt_cat) %>%
  na.omit() %>%
  count(lump_mode_trt_cat) %>% 
  mutate(proportion = round((n/sum(n)), digits=3)) %>% 
  arrange(proportion) %>%
  mutate(labels=scales::percent(proportion))

resourceFig <- ggplot(resourceWide, aes(x="", y=proportion, fill=lump_mode_trt_cat)) +
  geom_col() +
  coord_polar(theta="y") +
  scale_fill_manual(values=autumnalPalette)  +
  ggtitle('Resource Addition') +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')


otherWide <- modeTrt  %>% 
  mutate(trt_category=ifelse(trt_type %in% c('CO2','N','P','N*P','mult_nutrient','irr','K'), 'Resource',
                             ifelse(trt_type %in% c('drought','temp','herb_removal'), 'Stress',
                                    ifelse(trt_type=='multiple_trts', 'Mult. Trts', 'Other')))) %>%
  filter(trt_category=='Other') %>% 
  select(lump_mode_trt_cat) %>%
  na.omit() %>%
  count(lump_mode_trt_cat) %>% 
  mutate(proportion = round((n/sum(n)), digits=3)) %>% 
  arrange(proportion) %>%
  mutate(labels=scales::percent(proportion))

otherFig <- ggplot(otherWide, aes(x="", y=proportion, fill=lump_mode_trt_cat)) +
  geom_col() +
  coord_polar(theta="y") +
  scale_fill_manual(values=autumnalPalette)  +
  ggtitle('Other') +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')


multWide <- modeTrt  %>% 
  mutate(trt_category=ifelse(trt_type %in% c('CO2','N','P','N*P','mult_nutrient','irr','K'), 'Resource',
                             ifelse(trt_type %in% c('drought','temp','herb_removal'), 'Stress',
                                    ifelse(trt_type=='multiple_trts', 'Mult. Trts', 'Other')))) %>%
  filter(trt_category=='Mult. Trts') %>% 
  select(lump_mode_trt_cat) %>%
  na.omit() %>%
  count(lump_mode_trt_cat) %>% 
  mutate(proportion = round((n/sum(n)), digits=3)) %>% 
  arrange(proportion) %>%
  mutate(labels=scales::percent(proportion))

multFig <- ggplot(multWide, aes(x="", y=proportion, fill=lump_mode_trt_cat)) +
  geom_col() +
  coord_polar(theta="y") +
  scale_fill_manual(values=autumnalPalette)  +
  ggtitle('Mult. Trts') +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')


#grouped figure
ggarrange(stressFig, resourceFig, otherFig, multFig,
          ncol = 2, nrow = 2)

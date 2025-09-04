################################################################################
##  08a_plotSizeRelationship.R: Is codominance related to plot size at the site level?
##
##  Author: Kimberly Komatsu
##  Date created: September 2, 2025
################################################################################

source('code\\01_library.R')
source('code\\02_functions.R')

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34, color='black'),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34, color='black'),
             plot.title = element_text(size=40, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))



### calculate richness at plot level

plotSize <- readRDS("data/modeSite.rds") %>% 
  separate(exp_unit, into=c('site_code', 'project_name', 'community_type', 
                            'plot_id', 'treatment', 'calendar_year'), 
           sep='::', remove=F) %>% 
  left_join(readRDS("data/expInfo.rds"))



### calculate mode of modes for richness

# Calculate mode across years for all plots ----------------------------------------------------------

# create a dataframe of codominant group for plots with a single timepoint (because can't calculate mode of singleton)
singletonCodomPlotYear <- plotRichness %>% 
  group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>% 
  mutate(length=length(plot_id)) %>% 
  ungroup() %>% 
  filter(length==1) %>% 
  rename(plot_rich=richness) %>% 
  dplyr::select(database, site_code, project_name, community_type, plot_id, trt_type, treatment, plot_rich) 

# calculate mode across years for all plots
modePlotTrue <- plotRichness %>%  
  group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>% 
  reframe(plot_rich = DescTools::Mode(richness)) %>% # mode function must be capital here 
  ungroup() %>% 
  filter(!is.na(plot_rich)) %>%
  group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>% 
  summarise(plot_rich=round(mean(plot_rich), digits=0), .groups='drop') %>% # calculate mean for ties
  rbind(singletonCodomPlotYear)

# for plots with singleton ties for modes, calculate mean and round to nearest integer
multipleMode <- plotRichness %>% 
  select(database, site_code, project_name, community_type, plot_id, trt_type, treatment, richness, calendar_year) %>% 
  unique() %>% 
  full_join(modePlotTrue) %>% 
  filter(is.na(plot_rich)) %>%
  group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>% 
  summarise(plot_rich=round(mean(richness), digits=0), .groups='drop')

# bind dataframes for averaged ties and true modes at plot level
modePlot <- rbind(modePlotTrue, multipleMode) %>% 
  left_join(readRDS("data/modePlot.rds")) %>% 
  mutate(group=ifelse(plot_codom %in% c(2,3), 'codominated', 
               ifelse(plot_codom==1, 'monodominated', 'even')))


# Calculate mode across all control plots for each experiment ------------------

# create a dataframe of codominant group for experiments with a single plot (because can't calculate mode of singleton)
singletonCodomPlot <- modePlot %>% 
  select(-plot_codom) %>% 
  filter(trt_type=='control') %>% 
  group_by(database, site_code, project_name, community_type) %>% 
  mutate(length=length(community_type)) %>% 
  ungroup() %>% 
  filter(length==1) %>% 
  rename(site_rich=plot_rich) %>% 
  dplyr::select(database, site_code, project_name, community_type, site_rich) 

# calculate mode across plots for each experiment, dropping those with ties
modeSiteTrue <- modePlot %>%
  select(-plot_codom) %>% 
  filter(trt_type=='control') %>%
  group_by(database, site_code, project_name, community_type) %>% # mode generated from these
  reframe(site_rich = DescTools::Mode(plot_rich)) %>%  
  ungroup() %>% 
  group_by(database, site_code, project_name, community_type) %>% 
  summarise(site_rich=round(mean(site_rich), digits=0), .groups='drop') %>% # calculate mean for ties
  filter(!is.na(site_rich)) %>% 
  rbind(singletonCodomPlot)

# for plots with ties for modes, calculate mean and round to nearest integer
multipleModeProj <- modePlot %>% 
  select(-plot_codom) %>% 
  filter(trt_type=='control') %>% 
  select(database, site_code, project_name, community_type, plot_id, plot_rich) %>% 
  full_join(modeSiteTrue) %>% 
  filter(is.na(site_rich)) %>%
  group_by(database, site_code, project_name, community_type) %>% 
  summarise(site_rich=round(mean(plot_rich), digits=0), .groups='drop')

# bind dataframes for average ties and true modes at site level
modeSite <- rbind(modeSiteTrue, multipleModeProj) %>% 
  left_join(readRDS("data/envData.rds")) %>% 
  left_join(readRDS("data/modeSite.rds")) %>% 
  mutate(group=ifelse(mode_site %in% c(2,3), 'codominated', 
                      ifelse(mode_site==1, 'monodominated', 'even')))



# ordinal logistic regressions ------------------

allModel <- polr(group ~ richness, data = plotRichness, Hess=T) 
nullAllModel <- polr(group ~ 1, data = plotRichness, Hess=T)
summary(allModel)
lrt <- anova(allModel, nullAllModel)

plotModel <- polr(group ~ plot_rich, data = modePlot, Hess=T) 
nullPlotModel <- polr(group ~ 1, data = modePlot, Hess=T)
summary(plotModel)
lrt <- anova(plotModel, nullPlotModel)

siteModel <- polr(group ~ site_rich, data = modeSite, Hess=T) 
nullSiteModel <- polr(group ~ 1, data = modeSite, Hess=T)
summary(siteModel)
lrt <- anova(siteModel, nullSiteModel)



# anovas (backwards testing) ------------------
summary(allBackwardModel <- aov(richness ~ group, data=plotRichness))
TukeyHSD(allBackwardModel)

summary(plotBackwardModel <- aov(plot_rich ~ group, data=modePlot))
TukeyHSD(plotBackwardModel)

summary(siteBackwardModel <- aov(site_rich ~ group, data=modeSite))
TukeyHSD(siteBackwardModel)


# figures ------------------

plotRichness$group <- factor(plotRichness$group, levels=c('even', 'codominated', 'monodominated'))
modePlot$group <- factor(modePlot$group, levels=c('even', 'codominated', 'monodominated'))
modeSite$group <- factor(modeSite$group, levels=c('even', 'codominated', 'monodominated'))

allFig <- ggplot(data=plotRichness, aes(x=group, y=richness)) +
  geom_boxplot() +
  xlab('') + ylab('Richness (all plots and years)') +
  coord_flip()

plotModeFig <- ggplot(data=modePlot, aes(x=group, y=plot_rich)) +
  geom_boxplot() +
  xlab('') + ylab('Richness (plot mode)') +
  coord_flip()

siteModeFig <- ggplot(data=modeSite, aes(x=group, y=site_rich)) +
  geom_boxplot() +
  xlab('') + ylab('Richness (site mode)') +
  coord_flip()


grid.arrange(allFig, plotModeFig, siteModeFig, ncol=3, nrow=1)
#export at 4000x500



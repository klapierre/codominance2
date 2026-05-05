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



### get plot size data for each database

corre <- read.csv('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\1_first author\\codominance\\data\\CoRRE\\corre_plot_size.csv') %>% 
  rename(plot_size=plot_size_m2) %>% 
  select(site_code, project_name, community_type, plot_size)

gex <- read.csv('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\1_first author\\codominance\\data\\GEx\\GEx-metadata-with-other-env-layers-v2.csv') %>% 
  rename(plot_size=PlotSize,
         site_code=site) %>% 
  mutate(project_name=0, community_type=0) %>% 
  select(site_code, project_name, community_type, plot_size)

plotSize <- rbind(corre, gex) %>% 
  full_join(readRDS("data/modeSite.rds")) %>% 
  mutate(plot_size=ifelse(database=='nutnet', 1, plot_size),
         group=ifelse(lumpMode==1, 'monodominated',
               ifelse(lumpMode==2, 'codominated',
               'even'))) %>% 
  filter(!is.na(plot_size))



# ordinal logistic regressions ------------------

siteModel <- polr(as.factor(group) ~ plot_size, data = plotSize, Hess=T) 
nullSiteModel <- polr(as.factor(group) ~ 1, data = plotSize, Hess=T)
ctable <- coef(summary(siteModel))
pvals <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
cbind(ctable, "p value" = pvals)
exp(coef(siteModel))
exp(confint(siteModel))
anova(siteModel, nullSiteModel)


# anovas (backwards testing) ------------------

summary(siteBackwardModel <- aov(plot_size ~ as.factor(group), data = plotSize))
TukeyHSD(siteBackwardModel)


# figures ------------------

plotSize$group <- factor(plotSize$group, levels=c('even', 'codominated', 'monodominated'))

ggplot(data=barGraphStats(data=plotSize, variable="plot_size", byFactorNames=c("group")), aes(x=group, y=mean)) +
  geom_jitter(data=plotSize, aes(x=group, y=plot_size), size=4, color='grey',
              width = 0.3, height = 0 ) +
  geom_point(size=9) +
  # geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0, size=6) +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se), width=0.2, size=3) +
  xlab('') + ylab(expression("Plot Size (m"^2*")")) +
  coord_flip()
#export at 1500x500

ggplot(data=plotSize, aes(x=group, y=plot_size)) +
  geom_boxplot() +
  xlab('') + ylab(expression("Plot Size (m"^2*")")) +
  coord_flip()
#export at 1500x500

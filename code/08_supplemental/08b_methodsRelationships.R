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
             legend.title=element_blank(), legend.text=element_text(size=30))

### get data collection method for each database (does not include GEx) --------

method <- read.csv('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\1_first author\\codominance\\data\\CoRRE\\CoRRE_RelativeCoverMarch2024.csv') %>% 
  select(site_code, project_name, community_type, data_type) %>% 
  unique() %>% 
  mutate(database='corre')


### get plot size data for each database ---------------------------------------

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


### combine datasets -----------------------------------------------------------

sizeMethod <- plotSize %>% 
  full_join(method) %>% 
  mutate(data_type2=case_match(data_type,
                               "% cover" ~ "cover",
                               "ANPP" ~ "biomass",
                               "Biomass" ~ "biomass",
                               "Cover" ~ "cover",
                               "abundance" ~ "density",
                               "frequency" ~ "density",
                               "line" ~ "hits",
                               "percent cover" ~ "cover",
                               "percent pin hits" ~ "hits",
                               "pin" ~ "hits",
                               "pinhits" ~ "hits",
                               "pinpoint" ~ "hits",
                               "point_intercept" ~ "hits",
                               "proportion subplots present" ~ "other",
                               "Relative Cover" ~ "cover",
                               "relative density" ~ "density",
                               "0" ~ "other",
                               .default = data_type)) %>% 
  mutate(data_type2=ifelse(database=='nutnet', 'cover', data_type2)) %>% 
  filter(!is.na(database), !is.na(group))



# ordinal logistic regressions ------------------

# plot size
sizeModel <- polr(as.factor(group) ~ plot_size, data = sizeMethod, Hess=T) 
nullSizeModel <- polr(as.factor(group) ~ 1, data = sizeMethod, Hess=T)
ctable <- coef(summary(sizeModel))
pvals <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
cbind(ctable, "p value" = pvals)
exp(coef(sizeModel))
exp(confint(sizeModel))
anova(sizeModel, nullSizeModel)

# sampling method
methodModel <- polr(as.factor(group) ~ data_type2, data = filter(sizeMethod, !is.na(data_type2)), Hess=T) 
nullMethodModel <- polr(as.factor(group) ~ 1, data = filter(sizeMethod, !is.na(data_type2)), Hess=T)
ctable <- coef(summary(methodModel))
pvals <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
cbind(ctable, "p value" = pvals)
exp(coef(methodModel))
exp(confint(methodModel))
anova(methodModel, nullMethodModel)


# figures ------------------

sizeMethod$group <- factor(sizeMethod$group, levels=c('even', 'codominated', 'monodominated'))

# plot size figure

ggplot(data=barGraphStats(data=sizeMethod, variable="plot_size", byFactorNames=c("group")), aes(x=group, y=mean)) +
  geom_jitter(data=sizeMethod, aes(x=group, y=plot_size), size=4, color='grey',
              width = 0.3, height = 0 ) +
  geom_point(size=9) +
  # geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0, size=6) +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se), width=0.2, size=3) +
  xlab('') + ylab(expression("Plot Size (m"^2*")")) +
  coord_flip()
#export at 1500x500


# method figure

methodPlotData <- sizeMethod %>%
  mutate(dominance=factor(group, levels = c("monodominated", "codominated", "even"))) %>%
  count(data_type2, dominance) %>%
  group_by(data_type2) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

method_labels <- methodPlotData %>%
  group_by(data_type2) %>%
  summarise(n_total=sum(n), .groups="drop") %>%
  mutate(label=paste0(str_to_sentence(data_type2), "\n(", n_total, ")"))

ggplot(filter(methodPlotData, !is.na(data_type2)), aes(x=data_type2, y=proportion, fill=dominance)) +
  geom_col(width=0.75, color="white") +
  scale_fill_manual(values=c("monodominated"="#007BA7", "codominated"="#A63922", "even"="#D8B573")) +
  scale_x_discrete(labels=method_labels$label) +
  labs(x="Sampling Method", y="Proportion")
#export at 1500x1000
  


### get experiment duration data for each database ---------------------------------------

expDuration <- readRDS('data/codomSppList.rds') %>% 
  group_by(database, site_proj_comm) %>% 
  filter(treatment_year==max(treatment_year)) %>% 
  ungroup() %>% 
  select(-exp_unit, -calendar_year, -plot_id, -group, -num_group, -genus_species, -rank) %>% 
  unique() %>% 
  rename(exp_length=treatment_year)

gex <- read.csv('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\1_first author\\codominance\\data\\GEx\\GEx_cleaned_11June2020.csv') %>% 
  rename(site_code=site) %>% 
  select(site_code, exage) %>% 
  unique()

expDuration2 <- full_join(expDuration, gex) %>% 
  mutate(exp_length=ifelse(is.na(exp_length), exage, exp_length)) %>% 
  filter(!is.na(exp_length), exp_length>0) %>% 
  select(-exage) %>% 
  full_join(readRDS('data/modeTrt.rds')) %>% 
  filter(!is.na(lump_mode_trt_cat))

# ordinal logistic regressions ------------------

siteModel <- polr(as.factor(lump_mode_trt_cat) ~ exp_length, data = expDuration2, Hess=T) 
nullSiteModel <- polr(as.factor(lump_mode_trt_cat) ~ 1, data = expDuration2, Hess=T)
ctable <- coef(summary(siteModel))
pvals <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
cbind(ctable, "p value" = pvals)
exp(coef(siteModel))
exp(confint(siteModel))
anova(siteModel, nullSiteModel)


# anovas (backwards testing) ------------------

summary(siteBackwardModel <- aov(exp_length ~ as.factor(lump_mode_trt_cat), data = expDuration2))
TukeyHSD(siteBackwardModel)


# figures ------------------

expDuration2$lump_mode_trt_cat <- factor(expDuration2$lump_mode_trt_cat, levels=c('Even', 'Codominated', 'Monodominated'))

ggplot(data=expDuration2, aes(x=lump_mode_trt_cat, y=mean)) +
  geom_jitter(data=expDuration2, aes(x=lump_mode_trt_cat, y=exp_length), size=4, color='grey',
              width = 0.3, height = 0 ) +
  geom_point(size=9) +
  # geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0, size=6) +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se), width=0.2, size=3) +
  xlab('') + ylab('Experiment Length (years)') +
  coord_flip()
#export at 1500x500

ggplot(data=expDuration2, aes(x=lump_mode_trt_cat, y=exp_length)) +
  geom_boxplot() +
  xlab('') + ylab(expression("Experiment Length (years)")) +
  coord_flip()
#export at 1500x500

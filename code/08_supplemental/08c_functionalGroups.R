################################################################################
##  08c_functionalGroups.R: Are codominant species more likely to be in the same
##  functional group?
##
##  Authors: Kimberly Komatsu
##  Date created: 4/29/2026
################################################################################


# setup -------------------------------------------------------------------

rm(list = ls())
source("code/01_library.R")
source("code/02_functions.R")

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34, color='black'),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34, color='black'),
             plot.title = element_text(size=40, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=40))

# read data ---------------------------------------------------------------

# Q2ctlGroupsSite <- readRDS("data/Q2ctlGroupsSite.rds") %>% mutate(trt_type='control')
# Q3trtGroupsSite <- readRDS("data/Q3trtGroupsSite.rds")

functionalGroup <- readRDS("data/allTraits.rds") %>% 
  select(species, growth_form)

families <- read.csv('https://pasta.lternet.edu/package/data/eml/edi/1533/4/5ebbc389897a6a65dd0865094a8d0ffd') %>% 
  select(family, species) %>% unique()

RDFD <- readRDS("data/traitp_trt.rds")


# by functional group ---------------------------------------------------------------

#generate codom long list and associate functional groups
codomGroupsSite <- RDFD %>% 
  rownames_to_column() %>% 
  pivot_longer(cols=c(sp1, sp2), names_to='codom', values_to='species') %>%
  left_join(functionalGroup) %>% 
  select(-c(LDMC:stem_support), -c(photosynthetic_pathway:n_fixation_type)) %>% 
  mutate(growth_form=ifelse(growth_form %in% c('vine', 'woody'), 'other', growth_form)) %>% 
  select(-species, -codom) %>%
  group_by(rowname, site_proj_comm, trt_type, p) %>% 
  arrange(growth_form) %>%
  mutate(id=row_number()) %>% 
  ungroup() %>% 
  pivot_wider(names_from=id, names_prefix='sp', values_from=growth_form) %>% 
  mutate(groups=paste(sp1, sp2, sep='_'),
         ctl_trt=ifelse(trt_type=='control', 'Unmanipulated', 'Treatment'))

## RDFD beta-binomial by functional group pairings
m <- glmmTMB(cbind(n_obs, n_pool - n_obs) ~ ctl_trt*groups,
        data = codomGroupsSite,
        family = betabinomial(),
        weights = 1 - p_na)

Anova(m)
emmeans(m, ~ctl_trt*groups)

#by ctl_trt
codomGroupsPercent <- codomGroupsSite %>%
  group_by(ctl_trt, groups) %>% 
  summarize(count=length(groups), .groups='drop') %>%
  group_by(ctl_trt) %>% 
  mutate(sum=sum(count),
         percent=100*(count/sum)) %>% 
  ungroup()

ggplot(data=codomGroupsSite, aes(x=groups, y=p, fill=fct_rev(ctl_trt))) +
  geom_boxplot() +
  scale_fill_manual(values=c('white', '#666666')) +
  scale_x_discrete(limits=c('graminoid_graminoid', 'forb_graminoid', 'graminoid_other',
                            'forb_forb', 'forb_other', 'other_other'),
                   labels=c('G-G\n(49%,38%)', 'G-F\n(24%,35%)', 'G-O\n(18%,12%)', 
                            'F-F\n(6%,9%)', 'F-O\n(1%,5%)', 'O-O\n(3%,1%)')) +
  xlab('') + ylab('RDFD')

ggsave("FigG_RDFD_functGroup.png", width = 30, height = 10, dpi = 400)

# #without ctl_trt (no effect of ctl_trt in model)
# codomGroupsPercent <- codomGroupsSite %>%
#   group_by(groups) %>% 
#   summarize(count=length(groups), .groups='drop') %>%
#   mutate(sum=sum(count),
#          percent=100*(count/sum)) %>% 
#   ungroup()
# 
# ggplot(data=codomGroupsSite, aes(x=groups, y=p)) +
#   geom_boxplot() +
#   scale_x_discrete(limits=c('graminoid_graminoid', 'forb_graminoid', 'graminoid_other',
#                             'forb_forb', 'forb_other', 'other_other'),
#                    labels=c('G-G\n(41%)', 'G-F\n(32%)', 'G-O\n(14%)', 
#                             'F-F\n(8%)', 'F-O\n(4%)', 'O-O\n(2%)')) +
#   xlab('') + ylab('RDFD')


# by plant family ---------------------------------------------------------------

#generate codom long list and associate families
codomFamiliesSite <- RDFD %>% 
  rownames_to_column() %>% 
  pivot_longer(cols=c(sp1, sp2), names_to='codom', values_to='species') %>% 
  left_join(families) %>% 
  select(-species, -codom) %>%
  group_by(rowname, site_proj_comm, trt_type, dist) %>% 
  mutate(id=row_number()) %>% 
  ungroup() %>% 
  pivot_wider(names_from=id, names_prefix='sp', values_from=family) %>% 
  mutate(groups=ifelse(sp1==sp2, 'Same Family', 'Diff. Family'),
         ctl_trt=ifelse(trt_type=='control', 'Unmanipulated', 'Treatment'))

## RDFD beta-binomial by family pairings
n <- glmmTMB(cbind(n_obs, n_pool - n_obs) ~ ctl_trt*groups,
             data = codomFamiliesSite,
             family = betabinomial(),
             weights = 1 - p_na)

Anova(n)
emmeans(n, ~ctl_trt*groups)

#by ctl_trt
codomFamiliesPercent <- codomFamiliesSite %>%
  group_by(ctl_trt, groups) %>% 
  summarize(count=length(groups), .groups='drop') %>%
  group_by(ctl_trt) %>% 
  mutate(sum=sum(count),
         percent=100*(count/sum)) %>% 
  ungroup()

ggplot(data=codomFamiliesSite, aes(x=groups, y=p, fill=fct_rev(ctl_trt))) +
  geom_boxplot() +
  scale_fill_manual(values=c('white', '#666666')) +
  scale_x_discrete(limits=c('Same Family', 'Diff. Family'),
                   labels=c('Same Family\n(47%, 38%)', 'Diff. Family\n(53%, 62%)')) +
  xlab('') + ylab('RDFD') +
  theme(legend.position=c(0.15,0.92), legend.text=element_text(size=35))

ggsave("FigG_RDFD_family.png", width = 15, height = 10, dpi = 400)

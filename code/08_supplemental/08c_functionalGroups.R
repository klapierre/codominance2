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

Q2ctlGroupsSite <- readRDS("data/Q2ctlGroupsSite.rds") %>% mutate(trt_type='control')
Q3trtGroupsSite <- readRDS("data/Q3trtGroupsSite.rds")

functionalGroup <- readRDS("data/allTraits.rds") %>% 
  select(species, growth_form)

#generate codom long list and associate functional groups
codomGroupsSite <- rbind(Q2ctlGroupsSite, Q3trtGroupsSite) %>% 
  filter(!is.na(alpha2)) %>% 
  rownames_to_column() %>% 
  pivot_longer(cols=c(alpha1, alpha2, alpha3), names_to='codom', values_to='species') %>% 
  filter(!is.na(species)) %>% 
  group_by(rowname, site_code, project_name, community_type, trt_type) %>% 
  mutate(num_codom=length(species)) %>%
  ungroup() %>% 
  left_join(functionalGroup) %>% 
  mutate(growth_form=ifelse(growth_form %in% c('vine', 'woody', 'fern', 'lycophyte'), 'other', growth_form)) %>% 
  select(-species, -codom) %>%
  group_by(rowname, site_code, project_name, community_type, trt_type) %>% 
  arrange(growth_form) %>%
  mutate(id=row_number()) %>% 
  ungroup() %>% 
  pivot_wider(names_from=id, names_prefix='alpha', values_from=growth_form) %>% 
  mutate(groups=paste(alpha1, alpha2, alpha3, sep='_'),
         ctl_trt=ifelse(trt_type=='control', 'ctl', 'trt')) %>%
  group_by(num_codom, groups) %>% 
  summarize(count=length(groups), .groups='drop') %>% 
  mutate(drop=ifelse(groups %in% c('NA_NA_NA','forb_NA_NA','graminoid_NA_NA','other_NA_NA'), 1,
              ifelse(num_codom==3 & groups %in% c('forb_forb_NA','forb_graminoid_NA','forb_other_NA','graminoid_graminoid_NA',
                                              'graminoid_other_NA','other_other_NA'), 1, 0))) %>% 
  filter(drop!=1) %>% select(-drop) %>% 
  mutate(groups=recode(groups,
                       'forb_forb_NA' = 'forb_forb',
                       'forb_graminoid_NA' = 'forb_graminoid',
                       'forb_other_NA' = 'forb_other',
                       'graminoid_graminoid_NA' = 'graminoid_graminoid',
                       'graminoid_other_NA' = 'graminoid_other',
                       'other_other_NA' = 'other_other')) %>% 
  group_by(num_codom) %>% 
  mutate(sum=sum(count),
         percent=100*(count/sum)) %>% 
  ungroup()

ggplot(data=codomGroupsSite, aes(x=as.factor(num_codom), y=percent, fill=groups)) +
  geom_bar(stat='identity')




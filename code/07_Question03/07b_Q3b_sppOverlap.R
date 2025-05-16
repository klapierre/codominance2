################################################################################
##  07b_Q3b_sppOverlap.R: Determine the overlap of dominant species across treatment and control. 
##
##  Authors: K. Komatsu
##  Date created: May 6, 2025
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

allSppList         <- readRDS("data/allSppList.rds") 
modeSite           <- readRDS("data/modeSite.rds") 
modeTrt            <- readRDS("data/modeTrt.rds") 
expInfo            <- readRDS("data/expInfo.rds") 
Q2ctlGroupsSite    <- readRDS("data/Q2ctlGroupsSite.rds") %>% 
                        select(-codom_freq_proj) %>% 
                        rename(ctlm_alpha1=alpha1,
                               ctlm_alpha2=alpha2,
                               ctlm_alpha3=alpha3) 


# Most common species combinations for each treatment in final year of experiment ------------------

Q3trtGroupsSite <- allSppList %>% 
  left_join(expInfo) %>% 
  filter(trt_type!='control') %>% # filter out control plots
  filter(site_code!='ufrec.us', # filter out this site, which has no control plots
         !grepl("plant_mani", trt_type)) %>% # filter out any treatment that directly manipulates plant species  
  group_by(database, site_code, project_name, community_type) %>% 
  mutate(max_year=max(as.numeric(calendar_year)),
         experiment_length=max_year-as.numeric(calendar_year)+1) %>% # get experiment length
  ungroup() %>% 
  group_by(database, site_code, project_name, community_type) %>% 
  filter(calendar_year==max(calendar_year)) %>% # filter to final experiment year
  ungroup() %>%
  group_by(exp_unit) %>%
  mutate(group2=round(mean(num_group), digits=0)) %>% # average codom and round to nearest integer to fix ties
  ungroup() %>% 
  select(-num_group, -group) %>% 
  rename(num_group=group2) %>% 
  unique() %>% 
  filter(num_group!=4, # remove even communities
         rank<(num_group+1)) %>% # filter to codominating species only
  arrange(exp_unit, genus_species) %>% # arrange species within plots alphabetically
  group_by(exp_unit) %>%
  mutate(alphabetical_rank=paste('alpha', rank(genus_species), sep='')) %>%  # generate an alphabetical rank for each species group
  ungroup() %>% 
  dplyr::select(-rank) %>% 
  pivot_wider(names_from=alphabetical_rank, values_from=genus_species) %>% # generate a list of codominating species groups
  group_by(site_code, project_name, community_type, trt_type, plot_id, alpha1, alpha2, alpha3) %>% 
  summarise(codom_freq_plot=length(plot_id), .groups="drop") %>% # frequency of each species group within a plot through time
  group_by(site_code, project_name, community_type, trt_type, plot_id) %>% 
  filter(codom_freq_plot==max(codom_freq_plot)) %>%  # drop any species group that is not most common in a plot through time
  ungroup() %>% 
  group_by(site_code, project_name, community_type, trt_type, alpha1, alpha2, alpha3) %>% 
  summarise(codom_freq_proj=length(community_type), .groups="drop") %>% # frequency of each species group within an experimental treatment across plots
  group_by(site_code, project_name, community_type, trt_type) %>% 
  filter(codom_freq_proj==max(codom_freq_proj)) %>%  # drop any species group that is not most common within an experimental treatment across plots
  ungroup()

# saveRDS(Q3trtGroupsSite, file = "data/Q3trtGroupsSite.rds")



# Overlap between trt species and ctl species ----------------------------------------------------

allGroupsSite <- Q3trtGroupsSite %>% 
  select(-codom_freq_proj) %>% 
  rename(trtm_alpha1=alpha1,
         trtm_alpha2=alpha2,
         trtm_alpha3=alpha3) %>% 
  full_join(Q2ctlGroupsSite, relationship = 'many-to-many') %>% 
  rowwise() %>% 
  mutate(any_match = length(intersect(na.omit(c_across(starts_with("trtm_"))),
                                      na.omit(c_across(starts_with("ctlm_"))))) > 0) %>% 
  mutate(num_match = length(intersect(na.omit(c_across(starts_with("trtm_"))),
                                      na.omit(c_across(starts_with("ctlm_")))))) %>% 
  mutate(num_antimatch = length(setdiff(na.omit(c_across(starts_with("trtm_"))),
                                        na.omit(c_across(starts_with("ctlm_")))))) %>% 
  mutate(trt_na = sum(is.na(c_across(starts_with("trtm_")))),
         ctl_na = sum(is.na(c_across(starts_with("ctlm_"))))) %>% 
  mutate(trt_codom = sum(!is.na(c_across(starts_with("trtm_")))),
         ctl_codom = sum(!is.na(c_across(starts_with("ctlm_"))))) %>% 
  ungroup() %>% 
  mutate(match2=ifelse(num_match==0, 'none',
                ifelse(ctl_codom==1 & trt_codom==1 & num_match==1, 'full',
                ifelse(ctl_codom==2 & trt_codom==2 & num_match==2, 'full',
                ifelse(ctl_codom==3 & trt_codom==3 & num_match==3, 'full',     
                      'partial'))))) %>% 
  mutate(match=ifelse(is.na(ctlm_alpha1), 'NA',
               ifelse(is.na(trtm_alpha1), 'NA', 
                      match2))) %>% 
  select(site_code, project_name, community_type, trt_type, trtm_alpha1, trtm_alpha2, trtm_alpha3,
         ctlm_alpha1, ctlm_alpha2, ctlm_alpha3, trt_codom, ctl_codom, match) %>% 
  filter(trt_codom!=0, ctl_codom!=0)


# saveRDS(allGroupsSite, file = "data/allGroupsSite.rds")



# Effect of trt type on overlap between trt species and ctl species ----------------------------------------------------

overlapTable <- xtabs(~ match + trt_type, data = allGroupsSite)

print(chisq <- chisq.test(overlapTable))
# X-squared = 95.146, df = 24, p-value = 1.995e-10

mosaicplot(overlapTable, shade = TRUE, las=2,
           main = "overlapTable")

summaryOverlapTable <- as.data.frame(overlapTable) %>% 
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
                       ifelse(trt_type=='multiple_trts', 'Mult. Trts', as.character(trt_type)))))))))) %>% 
  mutate(match_nice=str_to_sentence(match))

summaryOverlapTable$trt_type_nice <- factor(summaryOverlapTable$trt_type_nice, 
                                            levels = c('N','P','K','N*P','Mult. Nutrients',
                                                       'Herbivore Rem.','Disturbance',
                                                       'CO2','Irrigation','Drought','Warming','Other',
                                                       'Mult. Trts'))

summaryOverlapTable$match_nice <- factor(summaryOverlapTable$match_nice, levels = c('Full', 'Partial', 'None'))

ggplot(summaryOverlapTable, aes(x=trt_type_nice , y=match_nice)) +
  geom_tile(aes(fill=percent)) +
  geom_text(aes(label=Freq), size=6, color='grey40') +
  # geom_text(aes(label=round(100*percent, digits=0))) +
  scale_fill_gradient(low='#F8FBF8', high='#031B88') +
  scale_y_discrete(limits=rev) +
  xlab('') + ylab('Species Overlap') + labs(fill='Column\nPercentage') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# ggsave(file='Fig6_heatMapOverlapTrt.png', width=10, height=4, units='in', dpi=300, bg='white')

overallOverlap <- summaryOverlapTable %>% 
  group_by(match_nice) %>% 
  summarize(count=sum(Freq), .groups='drop')

ggplot(overallOverlap, aes(x="", y=count, fill=match_nice)) +
  geom_col() +
  coord_polar(theta="y") +
  scale_fill_manual(values=c('#B0AAD1', '#F7695F', '#6EA1C9'))  +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        plot.title = element_text(vjust = 0.5),
        legend.position = 'none')  

# ggsave(file='Fig6a_pieOverlapTrt.png', width=4, height=4, units='in', dpi=300, bg='white')
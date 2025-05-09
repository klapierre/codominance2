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
               ifelse(ctl_codom==1 & num_match==1, 'full',
               ifelse(ctl_codom %in% c(2,3) & trt_codom==1 & num_match==1, 'full',
               ifelse(ctl_codom %in% c(2,3) & trt_codom==2 & num_match==2, 'full',
               ifelse(ctl_codom==3 & trt_codom==3 & num_match==3, 'full',     
                      'partial')))))) %>% 
  mutate(match=ifelse(is.na(ctlm_alpha1), 'NA',
               ifelse(is.na(trtm_alpha1), 'NA', 
                      match2))) %>% 
  select(site_code, project_name, community_type, trt_type, trtm_alpha1, trtm_alpha2, trtm_alpha3, ctlm_alpha1, ctlm_alpha2, ctlm_alpha3, trt_codom, ctl_codom, match)


# saveRDS(allGroupsSite, file = "data/allGroupsSite.rds")



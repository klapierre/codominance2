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
                        rename(ctl_alpha1=alpha1,
                               ctl_alpha2=alpha2,
                               ctl_alpha3=alpha3) 


# Most common species combinations for each treatment in final year of experiment ------------------

Q3trtGroupsSite <- allSppList %>% 
  left_join(expInfo) %>% 
  filter(trt_type!='control') %>% #filter out control plots
  #create trt_type variable that groups some of the most common treatments together
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
  mutate(trt_type=ifelse(trt_type3 %in% c('N','P','K','N*P','mult_nutrient','herb_removal','disturbance','CO2',
                                           'irr','drought','temp','multiple_trts', 'control'), trt_type3, 'other')) %>%
  filter(site_code!='ufrec.us', #filter out this site, which has no control plots
         !grepl("plant_mani", trt_type)) %>% #filter out any treatment that directly manipulates plant species  
  group_by(database, site_code, project_name, community_type) %>% 
  mutate(max_year=max(as.numeric(calendar_year)),
         experiment_length=max_year-as.numeric(calendar_year)+1) %>% 
  ungroup() %>% 
  group_by(database, site_code, project_name, community_type) %>% 
  filter(calendar_year==max(calendar_year)) %>% 
  ungroup() %>%
  group_by(exp_unit) %>%
  mutate(group2=round(mean(num_group), digits=0)) %>% #average codom and round to nearest integer to fix ties
  ungroup() %>% 
  select(-num_group, -group) %>% 
  rename(num_group=group2) %>% 
  unique() %>% 
  filter(num_group!=4, #remove even communities
         rank<(num_group+1)) %>% #filter to codominating species only
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
  rename(trt_alpha1=alpha1,
         trt_alpha2=alpha2,
         trt_alpha3=alpha3) %>% 
  full_join(Q2ctlGroupsSite, relationship = 'many-to-many') %>% 
  rowwise() %>% 
  mutate(any_match = length(intersect(
    na.omit(c_across(starts_with("trt_"))),
    na.omit(c_across(starts_with("ctl_"))))) > 0) %>%
  ungroup()


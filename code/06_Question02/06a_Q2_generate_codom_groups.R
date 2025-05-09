
################################################################################
##  06a_Q2_formatCodominants.R: Get codominating species from codom list.
##
##  Authors: Akira Terui, Jordan Winter, K. Komatsu
##  Date created: January 29, 2025
################################################################################

source("code/01_library.R")
source("code/02_functions.R")


# Read Data --------------------------------------------------------------------

codomSppList <- readRDS("data/codomSppList.rds") %>% 
  left_join(readRDS("data/expInfo.rds"))


# Most common species combinations for each plot through time ------------------

#filter to control plots only
codomControl <- codomSppList %>% 
  filter(trt_type=='control')

Q2ctlGroupsSite <- codomControl %>%  
  # group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>% 
  # reframe(plot_codom = DescTools::Mode(num_group)) %>% # mode function must be capital here 
  # ungroup() %>% 
  # left_join(codomControl) %>% 
  # filter(num_group==plot_codom) %>% # remove years where number of codominants in a plot was not the mode
  arrange(exp_unit, genus_species) %>% # arrange species within plots alphabetically
  group_by(exp_unit) %>%
  mutate(alphabetical_rank=paste('alpha', rank(genus_species), sep='')) %>%  # generate an alphabetical rank for each species group
  ungroup() %>% 
  dplyr::select(-rank) %>% 
  pivot_wider(names_from=alphabetical_rank, values_from=genus_species) %>% # generate a list of codominating species groups
  group_by(site_code, project_name, community_type, plot_id, alpha1, alpha2, alpha3) %>% 
  summarise(codom_freq_plot=length(plot_id), .groups="drop") %>% # frequency of each species group within a plot through time
  group_by(site_code, project_name, community_type, plot_id) %>% 
  filter(codom_freq_plot==max(codom_freq_plot)) %>%  # drop any species group that is not most common in a plot through time
  ungroup() %>% 
  group_by(site_code, project_name, community_type, alpha1, alpha2, alpha3) %>% 
  summarise(codom_freq_proj=length(community_type), .groups="drop") %>% # frequency of each species group within an experiment across plots
  group_by(site_code, project_name, community_type) %>% 
  filter(codom_freq_proj==max(codom_freq_proj)) %>%  # drop any species group that is not most common within an experiment across plots
  ungroup()

# saveRDS(Q2ctlGroupsSite, file = "data/Q2ctlGroupsSite.rds")

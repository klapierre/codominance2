
################################################################################
##  06a_Q2_formatCodominants.R: Get codominating species from codom list.
##
##  Authors: Akira Terui, Jordan Winter, Kimberly Komatsu
##  Date created: January 29, 2025
################################################################################

source("code/01_library.R")
source("code/02_functions.R")


# Read Data --------------------------------------------------------------------

codomSppList <- readRDS("data/allSppList.rds") %>% 
  left_join(readRDS("data/expInfo.rds")) %>% 
  filter(num_group!=4, # remove even communities
         rank<(num_group+1)) # filter to codominating species only
  

# Most common species combinations for each plot through time ------------------

#filter to control plots only
codomControl <- codomSppList %>% 
  filter(trt_type=='control') %>% 
  left_join(read_rds('data/modeSite.rds'))

Q2ctlGroupsSite <- codomControl %>%
  filter(num_group==mode_site) %>% #filter to only plots where the number of codominanting species matches the site pairings
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



# All possible species combinations (co-occur within at least one plot at one time point) within each site

allSppList <- readRDS("data/allSppList.rds")

allSppPairs <- allSppList %>%
  distinct(site_proj_comm, plot_id, calendar_year, genus_species) %>% 
  # Generate all species pairs within each plot-year
  group_by(site_proj_comm, plot_id, calendar_year) %>%
  filter(n() >= 2) %>%
  reframe(
    species_1 = combn(sort(genus_species), 2)[1, ],
    species_2 = combn(sort(genus_species), 2)[2, ]
  ) %>%
  ungroup() %>%
  distinct(site_proj_comm, species_1, species_2)

# saveRDS(allSppPairs, file = "data/allSppPairs.rds")
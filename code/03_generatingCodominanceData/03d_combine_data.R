################################################################################
##  03d_combine_data.R: Combine data from CoRRE, GEx, and NutNet and clean.
##
##  Author: Ashley LaRoque (modified K. Komatsu)
##  Date created: 
################################################################################

source("code/01_library.R")
source("code/02_functions.R")


# list of all spp and ranks: details on line 158
# expInfo <- readRDS("data/expInfo.rds")

# categorical groups of codoms: details on line 232
# numCodomPlotYear <- readRDS("data/numCodomPlotYear.rds")

# list of all spp and ranks: details on line 249
# allSppList <- readRDS("data/allSppList.rds")

# list of all spp and ranks: details on line 259
# codomSppList <- readRDS("data/codomSppList.rds")



#kim's laptop
setwd('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\first author\\2024_codominance')


# Read data ---------------------------------------------------------------

corre <- read.csv('data/CoRRE/corre_codominants_list_20250312.csv') %>%
  dplyr::select(-block, -genus_species, -relcov, -rank) %>%
  unique() %>%
  left_join(read_csv('data/corre/corre_richEven_20240208.csv')) %>%
  left_join(read_csv('data/corre/corre_plot_size.csv')) %>%
  left_join(read_csv('data/corre/corre_siteBiotic_2021.csv')) %>%
  left_join(read_csv('data/corre/corre_siteLocationClimate_2021.csv')) %>%
  left_join(read_csv('data/corre/corre_ExperimentInfo_2021.csv')) %>% 
  group_by(site_code, project_name, community_type, calendar_year, treatment) %>%
  mutate(plot_number=length(plot_id),
         database='corre') %>%
  ungroup() %>%
  group_by(site_code, project_name, community_type) %>%
  mutate(experiment_length = max(treatment_year)) %>%
  ungroup() %>%
  dplyr::select(database, exp_unit, site_code, project_name, community_type, plot_id, 
                calendar_year, treatment_year, experiment_length, treatment, trt_type,
                plot_size_m2, plot_number, plot_permenant, MAP, MAT, rrich, anpp, Cmax,
                num_codominants, richness, Evar) %>%
  rename(gamma_rich=rrich) %>% 
  filter(project_name != "NutNet") #remove NutNet data from CoRRE to avoid duplicates with NutNet database


unique(corre$trt_type)


#gex

gex <- read_csv('data/GEx/gex_codominants_list_20250312.csv') %>%
  dplyr::select(-genus_species, -relcov, -rank) %>%
  unique() %>% 
  left_join(read_csv('data/gex/gex_richEven_20240213.csv'))%>%
  left_join(read_csv('data/gex/gex-metadata-with-other-env-layers-v2.csv')) %>%
  unique() %>% 
  mutate(database='gex', 
         project_name='NA', 
         community_type='NA',
         trt_type=ifelse(trt=='G', 'control', 'herb_removal'), 
         plot_permenant='NA', 
         MAT=bio1/10)%>%
  rename(site_code = site,
         plot_id = block, 
         calendar_year = year, 
         treatment_year = exage,
         plot_id = block, 
         treatment = trt, 
         plot_size_m2 = PlotSize, 
         MAP = bio12,
         gamma_rich = sprich,
         anpp = ANPP) %>%
  group_by(site_code, plot_id, treatment) %>%
  mutate(experiment_length = max(treatment_year)) %>%
  ungroup() %>%
  group_by(site_code, project_name, community_type, calendar_year, treatment) %>%
  mutate(plot_number = length(plot_id)) %>%
  ungroup() %>%
  dplyr::select(database, exp_unit, site_code, project_name, community_type, plot_id, 
                calendar_year, treatment_year, treatment, trt_type, plot_size_m2, 
                plot_number, plot_permenant, MAP, MAT, gamma_rich, anpp, experiment_length,
                Cmax, num_codominants, richness, Evar)

unique(gex$trt_type)



#NutNet
nutnetANPP <- read_csv('data/NutNet/comb-by-plot-clim-soil-diversity_2023-11-07.csv') %>%
  filter(trt=='Control') %>%
  group_by(site_code, year) %>%
  summarise(anpp=mean((vascular_live_mass+nonvascular_live_mass), na.rm=T)) %>%
  ungroup() %>%
  filter(!is.nan(anpp), anpp>0) %>% #drops the data from some sites in years they didn't collect, and two years at one site that reported 0 growth
  group_by(site_code) %>%
  summarise(anpp=mean(anpp, na.rm=T)) %>%
  ungroup()

nutnetSiteInfo <- read_csv('data/NutNet/comb-by-plot-clim-soil-diversity_2023-11-07.csv') %>%
  group_by(site_code) %>%
  mutate(experiment_length=max(year_trt)) %>%
  ungroup() %>%
  group_by(site_code, year) %>%
  mutate(plot_number=length(plot)) %>%
  ungroup() %>%
  rename(MAP=MAP_v2, MAT=MAT_v2, gamma_rich=site_richness) %>%
  left_join(nutnetANPP) %>%
  dplyr::select(site_code, trt, year, MAP, MAT, gamma_rich, anpp, experiment_length, plot_number) %>%
  unique()

nutnet <- read_csv('data/NutNet/NutNet_codominants_list_plot_20250312.csv') %>%
  dplyr::select(exp_unit, site_code, plot, year, year_trt, trt, Cmax, num_codominants) %>%
  unique() %>%
  left_join(read_csv('data/NutNet/nutnet_plot_richEven_20240213.csv')) %>% 
  left_join(nutnetSiteInfo) %>%
  mutate(database='nutnet', project_name='NA', community_type='NA', plot_size_m2=1, plot_permenant='y',
         trt_type=ifelse(year_trt<1, 'control',
                  ifelse(trt=='Control', 'control',
                  ifelse(trt=='Fence', 'herb_removal', 
                  ifelse(trt=='NPK+Fence', 'mult_nutrient*herb_removal',
                  ifelse(trt=='N', 'N',
                  ifelse(trt=='P', 'P', 
                  ifelse(trt=='K', 'K', 
                  ifelse(trt=='NP', 'N*P', 
                         'mult_nutrient'))))))))) %>%
  rename(plot_id=plot, calendar_year=year, treatment_year=year_trt, treatment=trt) %>%
  dplyr::select(database, exp_unit, site_code, project_name, community_type, plot_id,
                calendar_year, treatment_year, treatment, trt_type, plot_size_m2, plot_number,
                plot_permenant, MAP, MAT, gamma_rich, anpp, experiment_length, Cmax, 
                num_codominants, richness, Evar)

unique(nutnet$trt_type)



# -----combine datasets-----

individualExperiments <- rbind(corre, gex, nutnet)

unique(individualExperiments$trt_type) # identify treatments

expInfo <- individualExperiments %>%
  mutate(trt_type2=ifelse(project_name=='IRG' & treatment=='i', 'irr',
                   ifelse(project_name=='IRG' & treatment=='c', 'control',
                          trt_type))) %>% 
  dplyr::select(exp_unit, database, site_code, project_name, community_type, plot_id, treatment, trt_type2,
                plot_size_m2, plot_number, plot_permenant, MAP, MAT, gamma_rich, anpp) %>%
  unique() %>% 
  rename(trt_type=trt_type2)

saveRDS(expInfo, file = "data/expInfo.rds") # saving derived data for analyses


#-----abundance cutoffs of codominance-----

correAbund <- read_csv('data/CoRRE/corre_codominantsRankAll_20250312.csv') %>% 
  dplyr::select(exp_unit, site_code, project_name, community_type, plot_id, treatment, calendar_year, genus_species, relcov, rank) %>% 
  filter(project_name != "NutNet") %>% #remove NutNet data from CoRRE to avoid duplicates with NutNet database
  mutate(database='corre')

gexAbund <- read_csv('data/GEx/gex_codominantsRankAll_20250312.csv') %>% 
  rename(site_code=site, treatment=trt, calendar_year=year) %>% 
  mutate(plot_id=paste(block, treatment, sep='_'),
         project_name=0, community_type=0) %>% 
  dplyr::select(exp_unit, site_code, project_name, plot_id, community_type, treatment, calendar_year, genus_species, relcov, rank) %>% 
  mutate(database='gex')

nutnetAbund <- read_csv('data/NutNet/NutNet_codominantsRankAll_20250312.csv') %>% 
  rename(calendar_year=year, plot_id=plot, treatment=trt) %>% 
  mutate(project_name=0, community_type=0) %>% 
  dplyr::select(exp_unit, site_code, project_name, community_type, plot_id, treatment, calendar_year, genus_species, relcov, rank) %>% 
  mutate(database='nutnet')

allAbund <- rbind(correAbund, gexAbund, nutnetAbund) %>% 
  dplyr::select(-site_code, -project_name, -community_type, -plot_id, -treatment, -calendar_year) %>% 
  full_join(individualExperiments) %>% 
  unique()

replicatesAll <- allAbund %>% dplyr::select(exp_unit) %>% unique() #66,783 individual data points (plot*year combinations)
replicatesSpatial <- allAbund %>% dplyr::select(site_code, project_name, community_type, plot_id) %>% unique() #11,127 individual plots
replicatesExperiment <- allAbund %>% dplyr::select(site_code, project_name, community_type) %>% unique() #542 experiments


#cutoff at 20% abundance for mean of all codominant species (only applied to plots with fewer than 4 codominant species, otherwise they will be set to "even" communities below); for monodominated plots, the single monodominant species cover must be greater than 20%
toFilter <- allAbund %>% 
  filter(num_codominants<4) %>% 
  group_by(exp_unit) %>% 
  filter(rank<(num_codominants+1)) %>% 
  summarise(mean_cover=mean(relcov)) %>% 
  ungroup() %>% 
  filter(mean_cover<20) %>% #2187 (3.3%) of these have mean cover of codominant species less than 20% and are being set to "even"
  mutate(category='even')

filterComplete <- toFilter %>%  
  full_join(allAbund) %>% 
  mutate(site_proj_comm=paste(site_code, project_name, community_type, sep='_')) %>% 
  mutate(num_codominants_fix=ifelse(is.na(category), num_codominants, 4)) %>% #make anything that fits our two filters into an even community (4+ codom), grouped below
  dplyr::select(-num_codominants) %>% 
  rename(num_codominants=num_codominants_fix)

filterMeanPlotLevel <- filterComplete %>% 
  select(exp_unit, num_codominants) %>% 
  unique()


# Create categorical codom groups -----------------------------------------

# group number of codominants into 4 categories 
df_grouped <- filterComplete %>% 
  mutate(group = case_when(num_codominants == 1 ~ "monodominated",
                           num_codominants == 2 ~ "codominated",
                           num_codominants == 3 ~ "tridominated",
                           num_codominants >= 4 ~ "even"),
         num_group = case_when(num_codominants == 1 ~ 1,
                               num_codominants == 2 ~ 2,
                               num_codominants == 3 ~ 3,
                               num_codominants >= 4 ~ 4))

df_plotLevel <- df_grouped %>% 
  select(exp_unit, num_codominants, num_group, group) %>% 
  unique()

saveRDS(df_plotLevel, file = "data/numCodomPlotYear.rds") # saving derived data for analyses

summary(as.factor(df_plotLevel$group)) #37,396 plots monodominated, 8996 plots even

# visualize groups
ggplot(df_plotLevel,
       aes(x=group)) +
  geom_histogram(stat='count', aes(x = factor(group, level = c('monodominated', 'codominated', 'tridominated', 'even')))) +
  theme_minimal()


# Generate species lists -----------------------------------------

allSppList <- df_grouped %>% 
  dplyr::select(database, site_proj_comm, exp_unit, site_code, project_name, community_type, calendar_year, treatment_year, 
                group, num_group, genus_species, rank)

saveRDS(allSppList, file = "data/allSppList.rds") # saving derived data for analyses

codomSppList <- df_grouped %>% 
  group_by(exp_unit) %>% 
  filter(rank<num_group+1) %>% # only keep species that were top dominants
  filter(group %in% c('codominated', 'tridominated')) %>% # only keep plots where the are codominanting spp (not even communities or monodominated communties)
  ungroup() %>% 
  dplyr::select(database, site_proj_comm, exp_unit, site_code, project_name, community_type, calendar_year, treatment_year, 
                group, num_group, genus_species, rank)

saveRDS(codomSppList, file = "data/codomSppList.rds") # saving derived data for analyses
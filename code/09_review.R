################################################################################
##  09_review.R: Sample dataset for code review.
##
##  Authors: Kimberly Komatsu
##  Date created: 6/16/2026
################################################################################


# setup -------------------------------------------------------------------

rm(list = ls())
source("code/01_library.R")
source("code/02_functions.R")


# develop sample dataset --------------------------------------------------

expList <- read.csv('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\1_first author\\codominance\\data/corre/corre_ExperimentInfo_March2024.csv') %>% 
  filter(public=='1') %>% 
  select(site_code, project_name, community_type) %>% 
  mutate(database='corre') %>% 
  unique()

# environmental data for each project
envData <- readRDS("data/envData.rds") %>% 
  semi_join(expList)

# experiment information, including treatments and environmental characteristics
expInfo <- readRDS("data/expInfo.rds") %>% 
  semi_join(expList)

#species names
sppNames <- read.csv('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\1_first author\\codominance\\data\\CoRRE\\corre2trykey_2021.csv') %>%
  select(genus_species, species_matched) %>%
  unique()

# experiment relative cover data
corre <- read.csv('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\1_first author\\codominance\\data\\CoRRE\\CoRRE_RawAbundanceMarch2024.csv') %>%
  # select(-X) %>%
  left_join(sppNames) %>% 
  semi_join(expList)


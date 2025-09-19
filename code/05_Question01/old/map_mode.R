### LOAD PACKAGES ###
library(tidyverse)
library(codyn)


#### LOAD NETWORK DATA ####
## CoRRE ####
CoRRERichProd <- read.csv("data/CoRRE/CoRRE_siteBiotic_2021.csv") %>%
  # keep community types w/in projects and sites separate so there are distinct ANPP and richness values for each, but the same coordinates, MAP, and MAT #
  mutate(site_proj_comm = paste(site_code, project_name, community_type, sep="_")) %>%
  rename("GDiv" = "rrich", "ANPP"="anpp")
#'rrich' is gamma div#

CoRRECoordClim <- read.csv("data/CoRRE/CoRRE_siteLocationClimate_2021.csv") %>%
  # remove Sil and SORBAS sites #
  filter(site_code!="Sil" & site_code!="SORBAS")

CoRREfull <- merge(CoRRECoordClim,CoRRERichProd, by=c("site_code"), all=T) %>%
  select(site_code, site_proj_comm, Latitude, Longitude, MAP, MAT, GDiv, ANPP) %>%
  # remove NutNet sites from CoRRE database #
  filter(site_proj_comm!="CDR_NutNet_0" & 
           site_proj_comm!="cbgb.us_NutNet_0" &
           site_proj_comm!="shps.us_NutNet_0" &
           site_proj_comm!="sier.us_NutNet_0" &
           site_proj_comm!="temple.us_NutNet_0" &
           site_proj_comm!="veluwe.nl_NutNet_0" &
           site_proj_comm!="yarra.au_NutNet_0")

#################

## GEx ####
GExfull <- read.csv("data/GEx/GEx-metadata-with-other-env-layers-v2.csv") %>%
  rename("site_code"="site", "Latitude"="Final.Lat", "Longitude"="Final.Long", "N_deposition"="N.deposition1993",
         "MAP"="precip", "GDiv"="sprich") %>%
  # 'bio1' is (MAT*100) from WorldClim in C #
  mutate(MAT = bio1/10) %>%
  select(site_code, Latitude, Longitude, MAP, MAT, GDiv, ANPP, N_deposition) 
# keep 'NA' for Kruger (Mananga, Shitbotawna, and Satara) sites as MAP not good predictor #
# 'precip' is MAP #
# 'sprich' is gamma div #

################

## NutNet ####
NutNetCoordClim <- read.csv("data/NutNet/comb-by-plot-clim-soil_2023-11-07.csv") %>%
  filter(trt=="Control") %>%
  mutate(N_Dep = as.numeric(N_Dep)) %>%
  # change site_code for 'yarra.au' since the same site_code exists in CoRRE so need to make sure they are not the same experiment #
  mutate(site_code = ifelse(site_code=="yarra.au", "yarra.au_NutNetdf", site_code)) %>%
  # change "NULL" habitat to NA # 
  mutate(habitat = ifelse(habitat=="NULL", NA, habitat)) %>%
  group_by(site_code, habitat)%>%
  # get average for site_codes to simplify dataframe - doesn't change the values since all the same anyways, just simplifies #
  summarise(Latitude=mean(latitude), Longitude=mean(longitude), GDiv=mean(site_richness),MAP=mean(MAP_v2), MAT=mean(MAT_v2), N_deposition=mean(N_Dep))
# 'site_richness' is gamma div #

NutNetprod <- read.csv("data/NutNet/full-biomass_2023-11-07.csv") %>% 
  filter(trt=="Control") %>%
  select(-year_trt, -trt, -live, -subplot, -site_name, -block) %>%
  # change site_code for 'yarra.au' since the same site_code exists in CoRRE so need to make sure they are not the same experiment (do the same for others after confirm) #
  mutate(site_code = ifelse(site_code=="yarra.au", "yarra.au_NutNetdf", site_code)) %>%
  group_by(year, site_code, plot) %>%
  spread(category,mass) %>%
  rename("FORB_Phlox_diffusa" = "FORB + PHLOX DIFFUSA") %>%
  group_by(year, site_code, plot) %>%
  # Ingrid said "LIVE" is all live biomass that doesn't fit into one of the other categories #
  mutate(TOTALPROD = sum(GRAMINOID, WOODY, FORB, LEGUME, PTERIDOPHYTE, VASCULAR, LIVE, ANNUAL, PERENNIAL, BRYOPHYTE, CACTUS, FORB_Phlox_diffusa, LICHEN,na.rm=T)) %>%
  select(year, site_code, plot, TOTAL, TOTALPROD) %>%
  gather(category, mass, 4:5) %>%
  group_by(year, site_code, plot) %>%
  summarise(anpp = ifelse(category=="TOTAL"|category=="TOTALPROD", mass, NA)) %>%
  # remove negative ANPP value from twostep.us before averaging #
  filter(anpp > 0) %>%
  group_by(site_code) %>%
  summarise(ANPP = mean(anpp))


NutNetfULL <- merge.data.frame(NutNetCoordClim, NutNetprod, by=c("site_code"), all=T)

# COMBINE DATA FROM ALL NETWORKS TOFETHER #
FinalAllNetworks <- bind_rows(CoRREfull, GExfull, NutNetfULL)

# write dataframe to csv to then input to ArcGIS Pro and join with HDI raster and NDep raster data from Oak Ridge and export table to .xlsx #
saveRDS(FinalAllNetworks, file = "data/forArcGIS_CodominanceAllNetworkData.rds")

#################################################################################################################

# read in data for mode #
df_mode_q1 <- readRDS("data_formatted/df_mode_q1.rds")

df_mode_ashley <- df_mode_q1 %>% 
  mutate(site_proj_comm = ifelse(project_name == "0",
                                 site_code,
                                 paste(site_code, project_name, community_type, sep= "_")))

# read in dataframe generated from ArcGIS that contains all site info, NDep, and HDI #
Codominance_AllSiteData <- read_csv("data/Codominance_AllSiteData.csv")
Codominance_AllSiteData <- Codominance_AllSiteData %>% 
  mutate(site_proj_comm = ifelse(is.na(site_proj_comm), site_code, site_proj_comm))

df_combined <- df_mode_ashley %>% 
  left_join(Codominance_AllSiteData, by = c("site_proj_comm", "site_code")) %>% 
  filter(site_proj_comm!="CDR_NutNet_0" & 
           site_proj_comm!="cbgb.us_NutNet_0" &
           site_proj_comm!="shps.us_NutNet_0" &
           site_proj_comm!="sier.us_NutNet_0" &
           site_proj_comm!="temple.us_NutNet_0" &
           site_proj_comm!="veluwe.nl_NutNet_0" &
           site_proj_comm!="yarra.au_NutNet_0")


write.csv(df_combined, "data_formatted/df_combined.csv", row.names=FALSE)
 # or for use in R 
saveRDS(df_combined, file = "data_formatted/df_combined.rds")

############################################################################################################






# 
#ggplot(df_combined, aes(x = as.factor(mode_yr), y = N_Deposition)) +
#geom_boxplot() +
# coord_flip()
# 
# 
# 
# corre_df <- read.csv("corre_codominantsRankAll_202402091.csv")
# gex_codominantsRankAll_202402091 <- read_csv("gex_codominantsRankAll_202402091.csv")
# NutNet_codominantsRankAll_20240213 <- read_csv("NutNet_codominantsRankAll_20240213.csv")
# 
# 
# 
# 
# # Brought in site level biotic and abiotic data from CoBuddies Group
# 
# 
# controlonlyCoRRe <- subset(corre_df, treatment %in% c("C",
#                                                       "c",
#                                                       "N0F0",
#                                                       "u_u_c",
#                                                       "Camb_Namb",
#                                                       "Control",
#                                                       "control",
#                                                       "CNA",
#                                                       "0",
#                                                       "CONTROL",
#                                                       "mixed_CO",
#                                                       "t1",
#                                                       "9",
#                                                       "9_f_u_n",
#                                                       "N0S0H0",
#                                                       "N0P0S0",
#                                                       "con",
#                                                       "CA_N1",
#                                                       "1_0_CO",
#                                                       "A",
#                                                       "uuuu",
#                                                       "_000",
#                                                       "CURRENT",
#                                                       "amb",
#                                                       "Open_Ungrazed",
#                                                       "1NF",
#                                                       "N0B0",
#                                                       "ref_rich16",
#                                                       "gcc",
#                                                       "N0",
#                                                       "CT",
#                                                       "AcAt",
#                                                       "N0M0",
#                                                       "x",
#                                                       "N0P0",
#                                                       "0_CONTROL",
#                                                       "CK_W0",
#                                                       "ambient",
#                                                       "ct",
#                                                       "CK",
#                                                       "N1P0",
#                                                       "ambient_control",
#                                                       "Cont",
#                                                       "XXX",
#                                                       "P3N0",
#                                                       "T0F0",
#                                                       "OO",
#                                                       "UnwarmedControl",
#                                                       "Reference",
#                                                       "CC",
#                                                       "U CC",
#                                                       "0_CONTROL_1",
#                                                       "0N0P")) %>% 
#   mutate(site_proj_comm = paste(site_code, project_name, community_type, sep="_")) %>% 
#   # remove Sil and SORBAS sites #
#   filter(site_code!="Sil" & site_code!="SORBAS") %>% 
#   # remove NutNet sites from CoRRE database #
#   filter(site_proj_comm!="CDR_NutNet_0" & 
#            site_proj_comm!="cbgb.us_NutNet_0" &
#            site_proj_comm!="shps.us_NutNet_0" &
#            site_proj_comm!="sier.us_NutNet_0" &
#            site_proj_comm!="temple.us_NutNet_0" &
#            site_proj_comm!="veluwe.nl_NutNet_0" &
#            site_proj_comm!="yarra.au_NutNet_0")
# 
# controlonlyCoRRe$block <- as.character(controlonlyCoRRe$block)
# 
# 
# controlonlyGex <- subset(gex_codominantsRankAll_202402091, trt %in% "U") %>% 
#   group_by(site, year) %>% 
#   #summarise(AVG_Site_Mode = mean(num_codominants)) %>% 
#   mutate(site_proj_comm = site) %>% 
#   ungroup()#averaging across blocks within a site within a year
# 
# 
# controlnlyNut <- subset(NutNet_codominantsRankAll_20240213, trt %in% "Control") %>% 
#   mutate(site_proj_comm = site_code)
# controlnlyNut$block <- as.character(controlnlyNut$block)
# 
# MegaBind <- bind_rows(controlnlyNut,controlonlyCoRRe,controlonlyGex)  %>% 
#   #mutate(site_proj_comm = ifelse(site_proj_comm == NA, site_code, site_proj_comm)) 
#   mutate(site_proj_comm = ifelse(is.na(site_proj_comm), site, site_proj_comm))  
# 
# 
# ModeCodom <- MegaBind %>%
#   group_by(site_proj_comm) %>% 
#   summarise(SiteDomMode = Mode(num_codominants)) %>% 
#   rename(site_proj_comm = "site_proj_comm")
# 
# 
# 
# MapData <- left_join(Codominance_AllSiteData, ModeCodom, by = "site_proj_comm")
# 
# write.csv(MapData,"MapData.csv")
# 
# NewData <- MapData %>% 
#   mutate(Doms = ifelse(SiteDomMode == 1, "Monodominated", ifelse(SiteDomMode == 2, "Codominated", "Even"))) %>% 
#   na.omit(MapData$SiteDomMode)
# NewData$Doms <- factor(NewData$Doms, levels = c("Monodominated", "Codominated", "Even"))
# 
# autumnalPalette <- c("#02385A", "#A63922", "#D8B573")
# 
# NewData %>% 
#   ggplot(aes(x="", y=SiteDomMode, fill = Doms)) +
#   geom_bar(stat = "identity", width=1) +
#   coord_polar("y",start=0)
# 
# ggplot(NewData, aes(x = Doms, fill = Doms)) +
#   geom_histogram(stat = "count") +
#   stat_count(binwidth = 1, 
#              geom = 'text', 
#              color = 'white', 
#              aes(label = after_stat(count)),
#              position = position_stack(vjust = 0.5))+
#   xlab("Site Dominance")+
#   ylab("Count")+
#   scale_fill_manual(values = autumnalPalette) +
#   theme(legend.position = "none")

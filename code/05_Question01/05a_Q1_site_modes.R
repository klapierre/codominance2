
################################################################################
##  05_Q1_site_modes.R: Calculate mode of codominance numbers in control plots and 
##  compare to environmental data.
##
##  Authors: Ashley LaRoque, Elise Grabda (modified K. Komatsu)
##  Date created: 
################################################################################

source("code/01_library.R")
source("code/02_functions.R")


# # includes categorical groups 'format_data'
#  df_grouped<- readRDS("data_formatted/df_grouped.rds")
# 
# # mode of codoms from 'format_data'
#  df_mode_q1 <- readRDS("data_formatted/df_mode_q1.rds")
#  
# # mode combined with map data from 'map_mode'
#  df_combined <- readRDS('data_formatted/df_combined.rds')
 
 
 
 # Read data ---------------------------------------------------------------
 
 numCodomPlotYear <- readRDS("data/numCodomPlotYear.rds") %>% 
   separate(exp_unit, into=c('site_code', 'project_name', 'community_type', 'plot_id', 'treatment', 'calendar_year'), sep='::', remove=F) %>% 
   left_join(readRDS("data/envData.rds")) %>% 
   left_join(readRDS("data/expInfo.rds"))
 
 
 # Calculate mode ----------------------------------------------------------
 
 # calculate mode across years for all plots
 modePlot <- numCodomPlotYear %>%  
   group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>% 
   reframe(plot_codom = Mode(num_group)) %>% # mode function must be capital here 
   ungroup()  
 
 
 # calculate mode across all control plots for each experiment
 modeSite <- modePlot %>%
   filter(trt_type=='control') %>%  
   group_by(site_code, project_name, community_type) %>% # mode generated from these
   reframe(mode_site = Mode(plot_codom)) %>%  
   ungroup()
 
 # saveRDS(modeSite, file = "data/modeSite.rds")


# boxplots for each predictor #
ggplot(df_combined, aes(x=as.factor(mode_yr), y=abs(Latitude)))+
  geom_boxplot()+
  coord_flip()

ggplot(doms, aes(x=as.factor(mode_yr), y=MAP))+
  geom_boxplot()+
  coord_flip()

ggplot(doms, aes(x=as.factor(mode_yr), y=MAT))+
  geom_boxplot()+
  coord_flip()

ggplot(doms, aes(x=as.factor(mode_yr), y=GDiv))+
  geom_boxplot()+
  coord_flip()

ggplot(doms, aes(x=as.factor(mode_yr), y=ANPP))+
  geom_boxplot()+
  coord_flip()

ggplot(doms, aes(x=as.factor(mode_yr), y=HumanDisturbance))+
  geom_boxplot()+
  coord_flip()

ggplot(doms, aes(x=as.factor(mode_yr), y=N_Deposition))+
  geom_boxplot()+
  coord_flip()

factors <-doms[,c(10,11,12,13,14,15)]
chart.Correlation(factors, method = "spearman")

summary()

# Analyses #
library(nnet, MASS)
a <- multinom(factor(df_combined$mode_yr,
                     levels = c(4,3,2,1)) ~ MAP+GDiv+HumanDisturbance+N_Deposition+MAP*MAT*ANPP+GDiv*N_Deposition*HumanDisturbance ,
              data=df_combined)
stepAIC(a, direction = "backward")

coef <- summary(a)$coefficients
coef

stderr <- summary(a)$standard.errors
stderr

z <- coef/stderr
p_values <- 2 * (1 - pnorm(abs(z)))
p_values

exp(coef(a))

head(round(fitted(a),2))

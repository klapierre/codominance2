
################################################################################
##  05_Q1_site_modes.R: Calculate mode of codominance numbers in control plots and 
##  compare to environmental data.
##
##  Authors: Ashley LaRoque, Elise Grabda (modified K. Komatsu)
##  Date created: 
################################################################################

source("code/01_library.R")
source("code/02_functions.R")


 # Read data ---------------------------------------------------------------
 
 numCodomPlotYear <- readRDS("data/numCodomPlotYear.rds") %>% 
   separate(exp_unit, into=c('site_code', 'project_name', 'community_type', 'plot_id', 'treatment', 'calendar_year'), sep='::', remove=F) %>% 
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
   ungroup() %>% 
   left_join(readRDS("data/envData.rds"))
 
 # saveRDS(modeSite, file = "data/modeSite.rds")


 # Visualizing each predictor with boxplots ----------------------------------------------------------
ggplot(modeSite, aes(x=as.factor(mode_site), y=abs(Latitude)))+
  geom_boxplot()+
  coord_flip()

ggplot(modeSite, aes(x=as.factor(mode_site), y=MAP))+
  geom_boxplot()+
  coord_flip()

ggplot(modeSite, aes(x=as.factor(mode_site), y=MAT))+
  geom_boxplot()+
  coord_flip()

ggplot(modeSite, aes(x=as.factor(mode_site), y=gamma_rich))+
  geom_boxplot()+
  coord_flip()

ggplot(modeSite, aes(x=as.factor(mode_site), y=anpp))+
  geom_boxplot()+
  coord_flip()

ggplot(modeSite, aes(x=as.factor(mode_site), y=HumanDisturbance))+
  geom_boxplot()+
  coord_flip()

ggplot(modeSite, aes(x=as.factor(mode_site), y=N_Deposition))+
  geom_boxplot()+
  coord_flip()

factors <-modeSite[,c(6:13)]
chart.Correlation(factors, method = "spearman")


# Multinomial Analysis ----------------------------------------------------------

a <- multinom(factor(modeSite$mode_site,
                     levels = c(4,3,2,1)) ~ MAP+gamma_rich+HumanDisturbance+N_Deposition+MAP*MAT*anpp+gamma_rich*N_Deposition*HumanDisturbance ,
              data=modeSite)
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



# Visualizing distribution of codominance across sites ----------------------------------------------------------

mapData <- modeSite %>% 
  mutate(codom_category = ifelse(mode_site == 1, "Monodominated", 
                          ifelse(mode_site == 2, "Codominated", 
                          ifelse(mode_site == 3, "Tridominated", 
                                 "Even")))) %>% 
  filter(!is.na(N_Deposition))

mapData$codom_category <- factor(mapData$codom_category, levels = c("Monodominated", "Codominated", "Tridominated", "Even"))

autumnalPalette <- c("#02385A", "#A63922", "#D8B573", 'grey')

mapData %>% 
  ggplot(aes(x="", y=mode_site, fill = codom_category)) +
  geom_bar(stat = "identity", width=1) +
  coord_polar("y",start=0)

ggplot(mapData, aes(x = codom_category, fill = codom_category)) +
  geom_histogram(stat = "count") +
  stat_count(binwidth = 1, 
             geom = 'text', 
             color = 'white', 
             aes(label = after_stat(count)),
             position = position_stack(vjust = 0.5))+
  xlab("Site Dominance")+
  ylab("Count")+
  scale_fill_manual(values = autumnalPalette) +
  theme(legend.position = "none")

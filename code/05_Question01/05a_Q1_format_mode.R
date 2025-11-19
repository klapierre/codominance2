
################################################################################
##  05a_Q1_site_modes.R: Calculate mode of codominance numbers in control plots and 
##  compare to environmental data.
##
##  Authors: Kimberly Komatsu, Jordan Winter
##  Date created: 
################################################################################

source("code/01_library.R")
source("code/02_functions.R")


# Read data ---------------------------------------------------------------
<<<<<<< HEAD
 
numCodomPlotYear <- readRDS("data/numCodomPlotYear.rds") %>% 
  separate(exp_unit, into=c('site_code', 'project_name', 'community_type', 
                            'plot_id', 'treatment', 'calendar_year'), 
           sep='::', remove=F) %>% 
  left_join(readRDS("data/expInfo.rds"))

 
# Calculate mode across years for all plots ----------------------------------------------------------
 
# create a dataframe of codominant group for plots with a single timepoint (because can't calculate mode of singleton)
singletonCodomPlotYear <- numCodomPlotYear %>% 
   group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>% 
   mutate(length=length(plot_id)) %>% 
   ungroup() %>% 
   filter(length==1) %>% 
   rename(plot_codom=num_group) %>% 
   dplyr::select(database, site_code, project_name, community_type, plot_id, trt_type, treatment, plot_codom) 
 
# calculate mode across years for all plots
modePlotTrue <- numCodomPlotYear %>%  
   group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>% 
   reframe(plot_codom = DescTools::Mode(num_group)) %>% # mode function must be capital here 
   ungroup() %>% 
   filter(!is.na(plot_codom)) %>%
   group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>% 
   summarise(plot_codom = round(mean(plot_codom), digits=0), .groups='drop') %>% # calculate mean for ties
   rbind(singletonCodomPlotYear)

# for plots with singleton ties for modes, calculate mean and round to nearest integer
multipleMode <- numCodomPlotYear %>% 
   select(database, site_code, project_name, community_type, plot_id, trt_type, treatment, num_group, calendar_year) %>% 
   unique() %>% 
   full_join(modePlotTrue) %>% 
   filter(is.na(plot_codom)) %>%
   group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>% 
   summarise(plot_codom = round(mean(num_group), digits=0), .groups='drop')
 
# bind dataframes for averaged ties and true modes at plot level
modePlot <- rbind(modePlotTrue, multipleMode)


# Calculate mode across all control plots for each experiment ------------------

# create a dataframe of codominant group for experiments with a single plot (because can't calculate mode of singleton)
singletonCodomPlot <- modePlot %>% 
  filter(trt_type == 'control') %>% 
  group_by(database, site_code, project_name, community_type) %>% 
  mutate(length = length(community_type)) %>% 
  ungroup() %>% 
  filter(length == 1) %>% 
  rename(mode_site = plot_codom) %>% 
  dplyr::select(database, site_code, project_name, community_type, mode_site) 

# calculate mode across plots for each experiment, dropping those with ties
modeSiteTrue <- modePlot %>%
   filter(trt_type == 'control') %>%
   group_by(database, site_code, project_name, community_type) %>% # mode generated from these
   reframe(mode_site = DescTools::Mode(plot_codom)) %>%  
   ungroup() %>% 
   group_by(database, site_code, project_name, community_type) %>% 
   summarise(mode_site = round(mean(mode_site), digits=0), .groups='drop') %>% # calculate mean for ties
   filter(!is.na(mode_site)) %>% 
   rbind(singletonCodomPlot)

# for plots with ties for modes, calculate mean and round to nearest integer
multipleModeProj <- modePlot %>% 
  filter(trt_type == 'control') %>% 
  select(database, site_code, project_name, community_type, plot_id, plot_codom) %>% 
  full_join(modeSiteTrue) %>% 
  filter(is.na(mode_site)) %>%
  group_by(database, site_code, project_name, community_type) %>% 
  summarise(mode_site = round(mean(plot_codom), digits=0), .groups='drop')

# bind dataframes for average ties and true modes at site level
modeSite <- rbind(modeSiteTrue, multipleModeProj) %>% 
  left_join(readRDS("data/envData.rds")) %>% 
  mutate(lumpMode = ifelse(mode_site == 3, 2, mode_site)) 

# how many sites?
siteCount <- modeSite %>%
  mutate(lat = round(Latitude, 1),
         long = round(Longitude, 1)) %>% 
  select(lat, long) %>% 
  unique()

# saveRDS(modePlot, file = "data/modePlot.rds")
# saveRDS(modeSite, file = "data/modeSite.rds")
# write.csv(modeSite, 'C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\1_first author\\codominance\\data\\modeSite_20251119.csv', row.names=F)

=======


# Multinomial Analysis ----------------------------------------------------------

#multinom function comes from nnet
#mutinomial model using monodominance as the reference state
multinom.baseline1 <- multinom(factor(df_iap$lumpMode,
                                      levels = c(1,2,4)) ~ Aridity + 
                                 cv_Precip +
                                 ANPP + 
                                 GDiv * 
                                 NDep * 
                                 HumanFootprint ,
                               data = df_iap)

#checking p-value manually
(coef <- summary(multinom.baseline1)$coefficients)

(stderr <- summary(multinom.baseline1)$standard.errors)

z <- coef/stderr
(p_values <- 2 * (1 - pnorm(abs(z))))



#PR2() from 'pscl' package provides many R2 estimates
pR2(multinom.baseline1) #McFaddenR2 = 0.337

#Putting into more convenient table format, P-values and OR's match model output 
results <- tidy(multinom.baseline1, conf.int = TRUE, exponentiate = TRUE) #tidy() is a 'broom' function

results_table <- results %>%
  mutate(OR_CI = sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high),
         p_value = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))) %>%
  select(y.level, term, OR_CI,p_value) %>%
  rename( "Outcome Category" = y.level,
          "Predictor" = term,
          "Odds Ratio (95% CI) vs. Monodominated" = OR_CI,
          "P-value" = p_value)

# Format: for figure of multinomial model predictions----------------------------------------------------------

# Generate mean (later used in model predicted data)
df_om <- df_iap %>% 
  na.omit() %>% 
  mutate(Aridity_mean = mean(Aridity),
         GDiv_mean = mean(GDiv),
         HumanFootprint_mean = mean(HumanFootprint),
         NDep_mean = mean(NDep),
         ANPP_mean = mean(ANPP),
         cv_Precip_mean = mean(cv_Precip)) 

# Variables of interest
var <- c("Aridity", "GDiv", "HumanFootprint", "NDep", "ANPP", "cv_Precip") 

# Generate sequence of values looped for each variable 
df_seq <- foreach(v = var, .combine = bind_cols) %do% {
  
  df_v <- df_om %>% 
    dplyr::select(v)
  
  out <- as_tibble(seq(from = min(df_v), 
                       to = max(df_v), 
                       length.out = 100))
  
}

# Set column names 
colnames(df_seq) <- var

# Add identifier column to link dataframes 
df_seq <- df_seq %>% 
  add_column(seq = seq(from = 1, to = 100, length.out = 100))

# Combine dataframes- variable means and sequence 
df_r <- df_om %>% 
  dplyr::select(ends_with("_mean")) %>% 
  slice(1:100) %>% 
  add_column(seq = seq(from = 1, to = 100, length.out = 100)) %>% 
  full_join(df_seq, by = "seq") 

# Mean of variables of interest
var2 <- c("Aridity_mean", "GDiv_mean", "HumanFootprint_mean", "NDep_mean", "ANPP_mean", "cv_Precip_mean") 

# Predict data using model, sequence, and mean
df_predicted <- foreach(v = var, 
                        v2 = var2,
                        .combine = bind_cols) %do% {
                          # Select variable sequence
                          df_v <- df_r %>% 
                            dplyr::select(v, seq)
                          # Select means for all other variables   
                          df_m <- df_r %>% 
                            dplyr::select(ends_with("_mean"), seq) %>% 
                            dplyr::select(!v2)
                          # Combine data sets to predict from   
                          df_s <- df_v %>% 
                            full_join(df_m, by = "seq") %>% 
                            dplyr::select(!seq) %>% 
                            rename_with(~str_remove(., '_mean'))
                          # Predict data
                          df_p <- cbind(df_s,
                                        data.frame(predict(multinom.baseline1,
                                                           newdata = df_s,
                                                           type = "probs")))
                          
                          
                          # Format for figure interpretation   
                          df_c <- df_p %>% 
                            dplyr::select(v, starts_with("X")) %>% 
                            pivot_longer(cols = starts_with("X"), 
                                         names_to = "Codom", 
                                         values_to = "Probability")
                        }



# figure out what is going on here and if it can be included in fig 2
#ci <- as.data.frame(confint(multinom.baseline1, level = 0.95))



# Clarify column names 
df_combined <- df_predicted %>% 
  select(Codom...2, Aridity, GDiv, HumanFootprint, NDep, ANPP, cv_Precip,
         starts_with("Probability")) %>% 
  rename(Codom = Codom...2,
         Aridity = Aridity,
         "Gamma Diversity" = GDiv,
         "Human Footprint Index" = HumanFootprint, 
         "N Deposition" = NDep,
         ANPP = ANPP,
         Precip = cv_Precip,
         Prob_Aridity = Probability...3,
         Prob_gamma = Probability...6,
         Prob_Human = Probability...9,
         Prob_N = Probability...12,
         Prob_ANPP = Probability...15,
         Prob_Precip = Probability...18)

# Assign names
prob <- c("Prob_Aridity", "Prob_gamma", "Prob_ANPP", "Prob_Human", "Prob_N", "Prob_Precip")
named_var <- c("Aridity", "Gamma Diversity", "ANPP", "Human Footprint Index", "N Deposition", "Precip")





# # Determining differences in codominance across ecoregions ----------------------------------------------------------
# 
# modeSiteEcoregion <- modeSite %>% 
#   mutate(lump_mode_cat=ifelse(lumpMode==1, 'Monodominated',
#                        ifelse(lumpMode==2, 'Codominated', 'Even'))) 
# 
# 
# ecoregionModel <- xtabs(~ Formation + lump_mode_cat, data = modeSiteEcoregion)
# 
# print(chisq <- chisq.test(ecoregionModel))
# # X-squared = 36.811, df = 28, p-value = 0.1231
# 
# mosaicplot(ecoregionModel, shade = TRUE, las=2,
#            main = "modeSiteCodom")
# 
# summaryModeSiteCodom <- as.data.frame(ecoregionModel) %>% 
#   group_by(lump_mode_cat ) %>%
#   mutate(percent=Freq/sum(Freq)) %>% 
#   ungroup()
# 
# ggplot(summaryModeSiteCodom, aes(x=lump_mode_cat , y=Formation)) +
#   geom_tile(aes(fill=percent)) +
#   geom_text(aes(label=Freq), size=6, color='grey40') +
#   # geom_text(aes(label=round(100*percent, digits=0))) +
#   scale_fill_gradient(low='#F8FBF8', high='#031B88') +
#   scale_y_discrete(limits=rev) +
#   xlab('Control') + ylab('Treatment') + labs(fill='Column\nPercentage')
>>>>>>> dfae5fa6874c4c2f11a89e4069a1f367ece05c86

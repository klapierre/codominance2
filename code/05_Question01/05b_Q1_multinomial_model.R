
################################################################################
##  05b_Q1_site_modes.R: Calculate mode of codominance numbers in control plots and 
##  compare to environmental data.
##
##  Authors: Ashley LaRoque, Jordan Winter, Elise Grabda, Alyssa Young, Rachael Brenneman (modified K. Komatsu)
##  Date created: 
################################################################################

source("code/01_library.R")
source("code/02_functions.R")


# Read data ---------------------------------------------------------------

modePlot <- readRDS("data/modePlot.rds")

df_iap <- readRDS("data/modeSiteRichness.rds") %>% 
  rename(ANPP=anpp,
         GDiv=gamma_rich,
         NDep=NDeposition,
         PlotRichness=site_rich_mean) %>% 
  mutate(z_MAP=scale(MAP),
         z_MAT=scale(MAT),
         z_GDiv=scale(GDiv),
         z_ANPP=scale(ANPP),
         z_HumanFootprint=scale(HumanFootprint),
         z_NDep=scale(NDep),
         z_PlotRichness=scale(PlotRichness))

# df_iap10 <- readRDS("data/modeSite_cutoff10.rds")%>% 
#   rename(ANPP=anpp,
#          GDiv=gamma_rich,
#          NDep=NDeposition)
# df_iap15 <- readRDS("data/modeSite_cutoff15.rds")%>% 
#   rename(ANPP=anpp,
#          GDiv=gamma_rich,
#          NDep=NDeposition)
# df_iap25 <- readRDS("data/modeSite_cutoff25.rds")%>% 
#   rename(ANPP=anpp,
#          GDiv=gamma_rich,
#          NDep=NDeposition)
# df_iap30 <- readRDS("data/modeSite_cutoff30.rds")%>% 
#   rename(ANPP=anpp,
#          GDiv=gamma_rich,
#          NDep=NDeposition)


# Multinomial Analysis ----------------------------------------------------------

#multinom function comes from nnet
#mutinomial model using monodominance as the reference state

multinom.baselineScale <- multinom(factor(df_iap$lumpMode,
                                          levels = c(1,2,4)) ~ 
                                     z_ANPP*(z_MAP + z_MAT + z_GDiv + z_PlotRichness) +
                                     z_NDep + z_HumanFootprint,
                                   data = df_iap)

# #cutoff10 
# multinom.baseline10 <- multinom(factor(df_iap$lumpMode,
#                                       levels = c(1,2,4)) ~ 
#                                   ANPP*(MAP + MAT + GDiv) +
#                                   NDep + HumanFootprint,
#                                data = df_iap10)
# #cutoff15
# multinom.baseline15 <- multinom(factor(df_iap$lumpMode,
#                                       levels = c(1,2,4)) ~ 
#                                   ANPP*(MAP + MAT + GDiv) +
#                                   NDep + HumanFootprint,
#                                data = df_iap15)
# #cutoff25
# multinom.baseline25 <- multinom(factor(df_iap$lumpMode,
#                                       levels = c(1,2,4)) ~
#                                   ANPP*(MAP + MAT + GDiv) +
#                                   NDep + HumanFootprint,
#                                data = df_iap25)
# #cutoff30
# multinom.baseline30 <- multinom(factor(df_iap$lumpMode,
#                                       levels = c(1,2,4)) ~ 
#                                   ANPP*(MAP + MAT + GDiv) +
#                                   NDep + HumanFootprint,
#                                data = df_iap30)

#checking p-value manually
(coef <- summary(multinom.baselineScale)$coefficients)

(stderr <- summary(multinom.baselineScale)$standard.errors)

z <- coef/stderr
(p_values <- 2 * (1 - pnorm(abs(z))))


#PR2() from 'pscl' package provides many R2 estimates
pR2(multinom.baselineScale) #McFaddenR2 = 0.393; our current cutoff exhibits the most explanatory power
pR2(multinom.baseline10) #McFaddenR2 = 0.311
pR2(multinom.baseline15) #McFaddenR2 = 0.293
pR2(multinom.baseline25) #McFaddenR2 = 0.314
pR2(multinom.baseline30) #McFaddenR2 = 0.293

#Putting into more convenient table format, P-values and OR's match model output 
results <- tidy(multinom.baselineScale, conf.int = TRUE, exponentiate = TRUE) #tidy() is a 'broom' function

results$y.level[results$y.level == 2] <- "Codominated"
results$y.level[results$y.level == 4] <- "Even"

results_table <- results %>%
  mutate(OR_CI = sprintf("%.2f ", estimate),
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
  mutate(MAP_mean = mean(MAP),
         MAT_mean = mean(MAT),
         GDiv_mean = mean(GDiv),
         HumanFootprint_mean = mean(HumanFootprint),
         NDep_mean = mean(NDep),
         PlotRichness_mean = mean(PlotRichness),
         ANPP_mean = mean(ANPP))

anpp_levels <- quantile(df_iap$z_ANPP, probs = c(0.10, 0.50, 0.90),
                        na.rm = TRUE, names = FALSE)

# Create a MAP × ANPP prediction grid
df_interaction <- tidyr::expand_grid(MAP = seq(min(df_om$MAP, na.rm = TRUE),
                                               max(df_om$MAP, na.rm = TRUE),
                                               length.out = 100),
                                     ANPP = anpp_levels) %>%
  mutate(MAT = mean(df_om$MAT, na.rm = TRUE),
         GDiv = mean(df_om$GDiv, na.rm = TRUE),
         PlotRichness = mean(df_om$PlotRichness, na.rm = TRUE),
         NDep = mean(df_om$NDep, na.rm = TRUE),
         HumanFootprint = mean(df_om$HumanFootprint, na.rm = TRUE))

# Variables of interest
var <- c("MAP", "MAT", "GDiv", "HumanFootprint", "NDep", "ANPP", "PlotRichness")

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
var2 <- c("MAP_mean", "MAT_mean", "GDiv_mean", "HumanFootprint_mean", "NDep_mean", "ANPP_mean", "PlotRichness_mean") 

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
                                        data.frame(predict(multinom.baselineScale,
                                                           newdata = df_s,
                                                           type = "probs")))
                          
                          
                          # Format for figure interpretation   
                          df_c <- df_p %>% 
                            dplyr::select(v, starts_with("X")) %>% 
                            pivot_longer(cols = starts_with("X"), 
                                         names_to = "Codom", 
                                         values_to = "Probability")
                        }



# Clarify column names 
df_combined <- df_predicted %>% 
  select(Codom...2, MAP, MAT, GDiv, HumanFootprint, NDep, ANPP, PlotRichness, # cv_Precip,
         starts_with("Probability")) %>% 
  rename(Codom = Codom...2,
         "MAP" = MAP,
         "MAT" = MAT,
         "Gamma Diversity" = GDiv,
         "Alpha Diversity" = PlotRichness,
         "Human Footprint Index" = HumanFootprint, 
         "N Deposition" = NDep,
         ANPP = ANPP,
         Prob_MAP = Probability...3,
         Prob_MAT = Probability...6,
         Prob_gamma = Probability...9,
         Prob_alpha = Probability...21,
         Prob_Human = Probability...12,
         Prob_N = Probability...15,
         Prob_ANPP = Probability...18)

# saveRDS(df_combined, file = "data/multimodalModelPredictions.rds") # saving derived data for analyses

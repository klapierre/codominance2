
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
df_iap <- readRDS("data/modeSite.rds") %>% 
  rename(ANPP=anpp,
         GDiv=gamma_rich,
         NDep=NDeposition)

df_iap10 <- readRDS("data/modeSite_cutoff10.rds")%>% 
  rename(ANPP=anpp,
         GDiv=gamma_rich,
         NDep=NDeposition)
df_iap15 <- readRDS("data/modeSite_cutoff15.rds")%>% 
  rename(ANPP=anpp,
         GDiv=gamma_rich,
         NDep=NDeposition)
df_iap25 <- readRDS("data/modeSite_cutoff25.rds")%>% 
  rename(ANPP=anpp,
         GDiv=gamma_rich,
         NDep=NDeposition)
df_iap30 <- readRDS("data/modeSite_cutoff30.rds")%>% 
  rename(ANPP=anpp,
         GDiv=gamma_rich,
         NDep=NDeposition)


# Multinomial Analysis ----------------------------------------------------------

#multinom function comes from nnet
#mutinomial model using monodominance as the reference state

multinom.baseline1 <- multinom(factor(df_iap$lumpMode,
                                      levels = c(1,2,4)) ~ 
                                 ANPP*(MAP + MAT + GDiv) +
                                 NDep + HumanFootprint,
                               data = df_iap)

#cutoff10 
multinom.baseline10 <- multinom(factor(df_iap$lumpMode,
                                      levels = c(1,2,4)) ~ 
                                  ANPP*(MAP + MAT + GDiv) +
                                  NDep + HumanFootprint,
                               data = df_iap10)
#cutoff15
multinom.baseline15 <- multinom(factor(df_iap$lumpMode,
                                      levels = c(1,2,4)) ~ 
                                  ANPP*(MAP + MAT + GDiv) +
                                  NDep + HumanFootprint,
                               data = df_iap15)
#cutoff25
multinom.baseline25 <- multinom(factor(df_iap$lumpMode,
                                      levels = c(1,2,4)) ~
                                  ANPP*(MAP + MAT + GDiv) +
                                  NDep + HumanFootprint,
                               data = df_iap25)
#cutoff30
multinom.baseline30 <- multinom(factor(df_iap$lumpMode,
                                      levels = c(1,2,4)) ~ 
                                  ANPP*(MAP + MAT + GDiv) +
                                  NDep + HumanFootprint,
                               data = df_iap30)

#checking p-value manually
(coef <- summary(multinom.baseline1)$coefficients)

(stderr <- summary(multinom.baseline1)$standard.errors)

z <- coef/stderr
(p_values <- 2 * (1 - pnorm(abs(z))))



#PR2() from 'pscl' package provides many R2 estimates
pR2(multinom.baseline1) #McFaddenR2 = 0.325; our current cutoff exhibits the most explanatory power
pR2(multinom.baseline10) #McFaddenR2 = 0.311
pR2(multinom.baseline15) #McFaddenR2 = 0.293
pR2(multinom.baseline25) #McFaddenR2 = 0.314
pR2(multinom.baseline30) #McFaddenR2 = 0.293

#Putting into more convenient table format, P-values and OR's match model output 
results <- tidy(multinom.baseline1, conf.int = TRUE, exponentiate = TRUE) #tidy() is a 'broom' function

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

gt(results_table, caption = "") %>% 
  gtsave("multinomial_cutoff20.png")

#Results10#Resultscaption = 10
results10 <- tidy(multinom.baseline10, conf.int = TRUE, exponentiate = TRUE) #tidy() is a 'broom' function

results_table10 <- results10 %>%
  mutate(OR_CI = sprintf("%.2f ", estimate),
         p_value = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))) %>%
  select(y.level, term, OR_CI,p_value) %>%
  rename( "Outcome Category" = y.level,
          "Predictor" = term,
          "Odds Ratio (95% CI) vs. Monodominated" = OR_CI,
          "P-value" = p_value)

gt(results_table10, caption = "10% Cutoff (Base)") %>% 
  gtsave("C:/Users/elise/Downloads/Tables/cutoff10.png")

#Results15
results15 <- tidy(multinom.baseline15, conf.int = TRUE, exponentiate = TRUE) #tidy() is a 'broom' function

results_table15 <- results15 %>%
  mutate(OR_CI = sprintf("%.2f ", estimate),
         p_value = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))) %>%
  select(y.level, term, OR_CI,p_value) %>%
  rename( "Outcome Category" = y.level,
          "Predictor" = term,
          "Odds Ratio (95% CI) vs. Monodominated" = OR_CI,
          "P-value" = p_value)

gt(results_table15, caption = "15% Cutoff (Base)") %>% 
  gtsave("C:/Users/elise/Downloads/Tables/cutoff15.png")

#Results25
results25 <- tidy(multinom.baseline25, conf.int = TRUE, exponentiate = TRUE) #tidy() is a 'broom' function

results_table25 <- results25 %>%
  mutate(OR_CI = sprintf("%.2f ", estimate),
         p_value = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))) %>%
  select(y.level, term, OR_CI,p_value) %>%
  rename( "Outcome Category" = y.level,
          "Predictor" = term,
          "Odds Ratio (95% CI) vs. Monodominated" = OR_CI,
          "P-value" = p_value)

gt(results_table25, caption = "25% Cutoff (Base)") %>% 
  gtsave("C:/Users/elise/Downloads/Tables/cutoff25.png")

#Results30
results30 <- tidy(multinom.baseline30, conf.int = TRUE, exponentiate = TRUE) #tidy() is a 'broom' function

results_table30 <- results30 %>%
  mutate(OR_CI = sprintf("%.2f ", estimate),
         p_value = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))) %>%
  select(y.level, term, OR_CI,p_value) %>%
  rename( "Outcome Category" = y.level,
          "Predictor" = term,
          "Odds Ratio (95% CI) vs. Monodominated" = OR_CI,
          "P-value" = p_value)

gt(results_table30, caption = "30% Cutoff (Base)") %>% 
  gtsave("C:/Users/elise/Downloads/Tables/cutoff30.png")


# Format: for figure of multinomial model predictions----------------------------------------------------------

# Generate mean (later used in model predicted data)
# change this data file to represent the cutoff of interest
df_om <- df_iap %>% 
  na.omit() %>% 
  mutate(MAP_mean = mean(MAP),
         MAT_mean = mean(MAT),
         GDiv_mean = mean(GDiv),
         HumanFootprint_mean = mean(HumanFootprint),
         NDep_mean = mean(NDep),
         ANPP_mean = mean(ANPP))

# Variables of interest
var <- c("MAP", "MAT", "GDiv", "HumanFootprint", "NDep", "ANPP")

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
var2 <- c("MAP_mean", "MAT_mean", "GDiv_mean", "HumanFootprint_mean", "NDep_mean", "ANPP_mean") 

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
  select(Codom...2, MAP, MAT, GDiv, HumanFootprint, NDep, ANPP, # cv_Precip,
         starts_with("Probability")) %>% 
  rename(Codom = Codom...2,
         "MAP" = MAP,
         "MAT" = MAT,
         "Gamma Diversity" = GDiv,
         "Human Footprint Index" = HumanFootprint, 
         "N Deposition" = NDep,
         ANPP = ANPP,
         Prob_MAP = Probability...3,
         Prob_MAT = Probability...6,
         Prob_gamma = Probability...9,
         Prob_Human = Probability...12,
         Prob_N = Probability...15,
         Prob_ANPP = Probability...18)

saveRDS(df_combined, file = "data/multimodalModelPredictions.rds") # saving derived data for analyses
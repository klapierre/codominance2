
################################################################################
##  05_Q1_site_modes.R: Calculate mode of codominance numbers in control plots and 
##  compare to environmental data.
##
##  Authors: Ashley LaRoque, Jordan Winter, Elise Grabda (modified K. Komatsu)
##  Date created: 
################################################################################

source("code/01_library.R")
source("code/02_functions.R")


# Read data ---------------------------------------------------------------
 
numCodomPlotYear <- readRDS("data/numCodomPlotYear.rds") %>% 
   separate(exp_unit, into=c('site_code', 'project_name', 'community_type', 'plot_id', 'treatment', 'calendar_year'), sep='::', remove=F) %>% 
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
   summarise(plot_codom=round(mean(plot_codom), digits=0), .groups='drop') %>% # calculate mean for ties
   rbind(singletonCodomPlotYear)

# for plots with singleton ties for modes, calculate mean and round to nearest integer
multipleMode <- numCodomPlotYear %>% 
   select(database, site_code, project_name, community_type, plot_id, trt_type, treatment, num_group, calendar_year) %>% 
   unique() %>% 
   full_join(modePlotTrue) %>% 
   filter(is.na(plot_codom)) %>%
   group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>% 
   summarise(plot_codom=round(mean(num_group), digits=0), .groups='drop')
 
# bind dataframes for averaged ties and true modes at plot level
modePlot <- rbind(modePlotTrue, multipleMode)


# Calculate mode across all control plots for each experiment ------------------

# create a dataframe of codominant group for experiments with a single plot (because can't calculate mode of singleton)
singletonCodomPlot <- modePlot %>% 
  filter(trt_type=='control') %>% 
  group_by(database, site_code, project_name, community_type) %>% 
  mutate(length=length(community_type)) %>% 
  ungroup() %>% 
  filter(length==1) %>% 
  rename(mode_site=plot_codom) %>% 
  dplyr::select(database, site_code, project_name, community_type, mode_site) 

# calculate mode across plots for each experiment, dropping those with ties
modeSiteTrue <- modePlot %>%
   filter(trt_type=='control') %>%
   group_by(database, site_code, project_name, community_type) %>% # mode generated from these
   reframe(mode_site = DescTools::Mode(plot_codom)) %>%  
   ungroup() %>% 
   group_by(database, site_code, project_name, community_type) %>% 
   summarise(mode_site=round(mean(mode_site), digits=0), .groups='drop') %>% # calculate mean for ties
   filter(!is.na(mode_site)) %>% 
   rbind(singletonCodomPlot)

# for plots with ties for modes, calculate mean and round to nearest integer
multipleModeProj <- modePlot %>% 
  filter(trt_type=='control') %>% 
  select(database, site_code, project_name, community_type, plot_id, plot_codom) %>% 
  full_join(modeSiteTrue) %>% 
  filter(is.na(mode_site)) %>%
  group_by(database, site_code, project_name, community_type) %>% 
  summarise(mode_site=round(mean(plot_codom), digits=0), .groups='drop')

# bind dataframes for average ties and true modes at site level
modeSite <- rbind(modeSiteTrue, multipleModeProj) %>% 
  left_join(readRDS("data/envData.rds")) %>% 
  mutate(lumpMode = ifelse(mode_site == 3, 2, mode_site))

 
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

multinom.baseline1 <- multinom(factor(modeSite$lumpMode,
                     levels = c(1,2,4)) ~ MAP+gamma_rich+HumanDisturbance+N_Deposition+MAP*MAT*anpp+gamma_rich*N_Deposition*HumanDisturbance ,
              data=modeSite)
a <- multinom(factor(modeSite$lumpMode,
                        levels = c(1,2,4)) ~ MAP+gamma_rich+HumanDisturbance+N_Deposition+MAP*MAT*anpp+gamma_rich*N_Deposition*HumanDisturbance ,
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

naomit <- modeSite %>% 
  na.omit()

predict.data.map <- data.frame(tibble(
  MAP = seq(from = min(naomit$MAP), to = max(naomit$MAP), length.out = 100), 
  MAT = mean(naomit$MAT),
  gamma_rich = mean(naomit$gamma_rich),
  HumanDisturbance = mean(naomit$HumanDisturbance),
  N_Deposition = mean(naomit$N_Deposition),
  anpp = mean(naomit$anpp)
))
predict.data.mat <- data.frame(tibble(
  MAT = seq(from = min(naomit$MAT), to = max(naomit$MAT), length.out = 100), 
  MAP = mean(naomit$MAP),
  gamma_rich = mean(naomit$gamma_rich),
  HumanDisturbance = mean(naomit$HumanDisturbance),
  N_Deposition = mean(naomit$N_Deposition),
  anpp = mean(naomit$anpp)
))
predict.data.gamma <- data.frame(tibble(
  gamma_rich = seq(from = min(naomit$gamma_rich), to = max(naomit$gamma_rich), length.out = 100), 
  MAT = mean(naomit$MAT),
  MAP = mean(naomit$MAP),
  HumanDisturbance = mean(naomit$HumanDisturbance),
  N_Deposition = mean(naomit$N_Deposition),
  anpp = mean(naomit$anpp)
))
predict.data.hdisturb <- data.frame(tibble(
  HumanDisturbance = seq(from = min(naomit$HumanDisturbance), to = max(naomit$HumanDisturbance), length.out = 100), 
  MAT = mean(naomit$MAT),
  gamma_rich = mean(naomit$gamma_rich),
  MAP = mean(naomit$MAP),
  N_Deposition = mean(naomit$N_Deposition),
  anpp = mean(naomit$anpp)
))
predict.data.ndep <- data.frame(tibble(
  N_Deposition = seq(from = min(naomit$N_Deposition), to = max(naomit$N_Deposition), length.out = 100), 
  MAT = mean(naomit$MAT),
  gamma_rich = mean(naomit$gamma_rich),
  HumanDisturbance = mean(naomit$HumanDisturbance),
  MAP = mean(naomit$MAP),
  anpp = mean(naomit$anpp)
))
predict.data.anpp <- data.frame(tibble(
  anpp = seq(from = min(naomit$anpp), to = max(naomit$anpp), length.out = 100), 
  MAT = mean(naomit$MAT),
  gamma_rich = mean(naomit$gamma_rich),
  HumanDisturbance = mean(naomit$HumanDisturbance),
  N_Deposition = mean(naomit$N_Deposition),
  MAP = mean(naomit$MAP)
))

predict.map <- cbind(predict.data.map,data.frame(predict(multinom.baseline1, newdata = predict.data.map, type = "probs")))
predict.mat <- cbind(predict.data.mat,data.frame(predict(multinom.baseline1, newdata = predict.data.mat, type = "probs")))
predict.ndep <- cbind(predict.data.ndep,data.frame(predict(multinom.baseline1, newdata = predict.data.ndep, type = "probs")))
predict.hdisturb <- cbind(predict.data.hdisturb,data.frame(predict(multinom.baseline1, newdata = predict.data.hdisturb, type = "probs")))
predict.gdiv <- cbind(predict.data.gamma,data.frame(predict(multinom.baseline1, newdata = predict.data.gamma, type = "probs")))
predict.anpp <- cbind(predict.data.anpp,data.frame(predict(multinom.baseline1, newdata = predict.data.anpp, type = "probs")))

p.gather.map <- gather(predict.map, key= "codom", value = "Probability", 7:9)
p.gather.mat <- gather(predict.mat, key= "codom", value = "Probability", 7:9)
p.gather.ndep <- gather(predict.ndep, key= "codom", value = "Probability", 7:9)
p.gather.hdisturb <- gather(predict.hdisturb, key= "codom", value = "Probability", 7:9)
p.gather.gdiv <- gather(predict.gdiv, key= "codom", value = "Probability", 7:9)
p.gather.anpp <- gather(predict.anpp, key= "codom", value = "Probability", 7:9)

grid.arrange(
p.gather.map %>% 
  ggplot(aes(x= MAP, y= Probability, colour = codom))+
  geom_line()+
  labs(y= "Probability")+
  theme(legend.position = "none"),
p.gather.mat %>% 
  ggplot(aes(x= MAT, y= Probability, colour = codom))+
  geom_line()+
  labs(y= "Probability")+
  theme(legend.position = "none"),
p.gather.ndep %>% 
  ggplot(aes(x= N_Deposition, y= Probability, colour = codom))+
  geom_line()+
  labs(y= "Probability")+
  theme(legend.position = "none"),
p.gather.hdisturb %>% 
  ggplot(aes(x= HumanDisturbance, y= Probability, colour = codom))+
  geom_line()+
  labs(y= "Probability")+
  theme(legend.position = "none"),
p.gather.gdiv %>% 
  ggplot(aes(x= gamma_rich, y= Probability, colour = codom))+
  geom_line()+
  labs(y= "Probability")+
  theme(legend.position = "none"),
p.gather.anpp %>% 
  ggplot(aes(x= anpp, y= Probability, colour = codom))+
  geom_line()+
  labs(y= "Probability")+
  theme(legend.position = "none")
)

#### random forest modeling attempt ####
library(randomForest)
library(caret)

predictors <- c("MAP", "MAT", "gamma_rich",
                "anpp", "HumanDisturbance", "N_Deposition")

response <- "lumpMode"

select_predictors <- modeSite %>% 
  select(MAP, MAT, gamma_rich, anpp, HumanDisturbance, N_Deposition, lumpMode) %>% 
  na.omit()

# Split data into training (80%) and testing (20%) sets
set.seed(120)
trainIndex <- createDataPartition(select_predictors$lumpMode, p = 0.8, list = FALSE)
trainData <- select_predictors[trainIndex, ]
testData <- select_predictors[-trainIndex, ]

# Train the Random Forest model
rf_model <- randomForest(lumpMode ~ ., data = trainData[, c(response, predictors)], 
                         ntree = 1000, importance = TRUE)

# Print model summary
print(rf_model)

# Feature Importance Plot
importance_df <- data.frame(Feature = rownames(importance(rf_model)), 
                            Importance = importance(rf_model)[, 1])

ggplot(importance_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Variable Importance in Random Forest Model", 
       x = "Predictor Variable", y = "Importance Score") +
  theme_minimal()

# Predict on test set
testData$Predicted <- predict(rf_model, newdata = testData[, predictors])

# Compute R-squared
r2 <- cor(testData$lumpMode, testData$Predicted)^2

# Scatter Plot: Observed vs Predicted
ggplot(testData, aes(x = lumpMode, y = Predicted)) +
  geom_point(alpha = 0.7, color = "darkblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = paste("Observed vs. Predicted (R? =", round(r2, 2), ")"),
       x = "Observed Nodule Number", y = "Predicted Nodule Number") +
  theme_minimal()


# code clean up  from 154 to 251----------------------------------------------------------

df_om <- modeSite %>% 
  na.omit()

var <- c("MAP", "MAT", "gamma_rich", "HumanDisturbance", "N_Deposition", "anpp")

df_seq <- foreach(v = var, .combine = bind_cols) %do% {
  
  df_v <- df_om %>% 
    select(v)
  
  out <- seq(from = min(df_v), 
             to = max(df_v), 
             length.out = 100)
  
}
colnames(df_seq) <- c("MAP", "MAT", "gamma_rich", "HumanDisturbance", "N_Deposition", "anpp")

df_predicted <- foreach(v = var, .combine = bind_cols) %do% {
 
   df_s <- df_seq %>% 
    select(v)
   
   df_p <- cbind(df_s, 
                 data.frame(predict(multinom.baseline1, 
                                    newdata = df_seq, # values generated differ from previous code.... investigate
                                    type = "probs")))
  
}


df_comb <- df_predicted %>% 
  pivot_longer(cols = starts_with("X"), names_to = "Codom", values_to = "Probability") %>%
  mutate(Codom = case_when(Codom == "X1...2" ~ "X1",
                           Codom == "X2...3" ~ "X2",
                           Codom == "X4...4" ~ "X4",
                           Codom == "X1...6" ~ "X1",
                           Codom == "X2...7" ~ "X2",
                           Codom == "X4...8" ~ "X4",
                           Codom == "X1...10" ~ "X1",
                           Codom == "X2...11" ~ "X2",
                           Codom == "X4...12" ~ "X4",
                           Codom == "X1...14" ~ "X1",
                           Codom == "X2...15" ~ "X2",
                           Codom == "X4...16" ~ "X4",
                           Codom == "X1...18" ~ "X1",
                           Codom == "X2...19" ~ "X2",
                           Codom == "X4...20" ~ "X4",
                           Codom == "X1...22" ~ "X1",
                           Codom == "X2...23" ~ "X2",
                           Codom == "X4...24" ~ "X4")) # needs to be redone but i have tried for hours and cannot get it to work the same

df_comb <- df_comb %>% 
  pivot_longer(cols = c("MAP", "MAT", "gamma_rich", "HumanDisturbance", "N_Deposition", "anpp"),
                 names_to = "Variable", values_to = "Value")

ggplot(df_comb, aes(x = Value, y = Probability, color = Codom)) +
  geom_line() +  
  facet_wrap(~ Variable,
             scales = "free") +
  labs(y = "Probability" ) 
 # theme(legend.position = "none")

p.gather.map %>% 
  ggplot(aes(x= MAP, y= Probability, colour = codom))+
  geom_line()+
  labs(y= "Probability")+
  theme(legend.position = "none")

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


################################################################################
##  05_Q1_site_modes.R: Calculate mode of codominance numbers in control plots and 
##  compare to environmental data.
##
##  Authors: Ashley LaRoque, Jordan Winter, Elise Grabda, Alyssa Young, Rachael Brenneman (modified K. Komatsu)
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


# Histogram- count of codom level per variable ----------------------------

df_hist <- modeSite %>% 
  select(mode_site, lumpMode, MAP, MAT, gamma_rich, anpp, HumanDisturbance, N_Deposition) %>% 
  pivot_longer(cols = c("MAP", "MAT", "gamma_rich", "anpp", "HumanDisturbance", "N_Deposition"), 
               names_to = "variable", values_to = "value") %>% 
  mutate(lumpMode = ifelse(lumpMode == 1, "Monodominated", 
                           ifelse(lumpMode == 2, "Codominated", "Even")),
         mode.levels = factor(lumpMode, levels = c("Monodominated", "Codominated", "Even")))

mode.labs <- c("monodominated", "codominated", "even")
names(mode.labs) <- c("Monodominated", "Codominated", "Even")

# bar counts are all the same per variable ...?
# this is at the site level 
ggplot(df_hist1, aes(mode.levels)) +
  geom_bar(aes(fill = mode.levels)) +
  facet_wrap(~ variable1,
             labeller = labeller(mode.levels = mode.labs)) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_manual(name = "",
                    labels = c("Monodominated", "Codominated", "Even"),
                    values = c("#02385A", "#A63922", "#D8B573")) +
  labs(x = "Codominance Level",
       y = "Count")

# histogram of variable values
ggplot(df_hist1, aes(value)) + # df_hist comes from formatted df above 
  geom_histogram(aes(fill = mode.levels), bins = 50) +
  facet_wrap(~ variable1,
             scales = "free") +
  theme_bw() +
  theme(legend.position = "top") +
  scale_fill_manual(name = "",
                    labels = c("Monodominated", "Codominated", "Even"),
                    values = c("#02385A", "#A63922", "#D8B573")) +
  labs(x = "Value",
       y = "Count")

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

##########################################################################
## DON'T USE THIS BLOCK OF CODE ANYMORE ####
# naomit <- modeSite %>% 
#   na.omit()
# 
# predict.data.map <- data.frame(tibble(
#   MAP = seq(from = min(naomit$MAP), to = max(naomit$MAP), length.out = 100), 
#   MAT = mean(naomit$MAT),
#   gamma_rich = mean(naomit$gamma_rich),
#   HumanDisturbance = mean(naomit$HumanDisturbance),
#   N_Deposition = mean(naomit$N_Deposition),
#   anpp = mean(naomit$anpp)
# ))
# predict.data.mat <- data.frame(tibble(
#   MAT = seq(from = min(naomit$MAT), to = max(naomit$MAT), length.out = 100), 
#   MAP = mean(naomit$MAP),
#   gamma_rich = mean(naomit$gamma_rich),
#   HumanDisturbance = mean(naomit$HumanDisturbance),
#   N_Deposition = mean(naomit$N_Deposition),
#   anpp = mean(naomit$anpp)
# ))
# predict.data.gamma <- data.frame(tibble(
#   gamma_rich = seq(from = min(naomit$gamma_rich), to = max(naomit$gamma_rich), length.out = 100), 
#   MAT = mean(naomit$MAT),
#   MAP = mean(naomit$MAP),
#   HumanDisturbance = mean(naomit$HumanDisturbance),
#   N_Deposition = mean(naomit$N_Deposition),
#   anpp = mean(naomit$anpp)
# ))
# predict.data.hdisturb <- data.frame(tibble(
#   HumanDisturbance = seq(from = min(naomit$HumanDisturbance), to = max(naomit$HumanDisturbance), length.out = 100), 
#   MAT = mean(naomit$MAT),
#   gamma_rich = mean(naomit$gamma_rich),
#   MAP = mean(naomit$MAP),
#   N_Deposition = mean(naomit$N_Deposition),
#   anpp = mean(naomit$anpp)
# ))
# predict.data.ndep <- data.frame(tibble(
#   N_Deposition = seq(from = min(naomit$N_Deposition), to = max(naomit$N_Deposition), length.out = 100), 
#   MAT = mean(naomit$MAT),
#   gamma_rich = mean(naomit$gamma_rich),
#   HumanDisturbance = mean(naomit$HumanDisturbance),
#   MAP = mean(naomit$MAP),
#   anpp = mean(naomit$anpp)
# ))
# predict.data.anpp <- data.frame(tibble(
#   anpp = seq(from = min(naomit$anpp), to = max(naomit$anpp), length.out = 100), 
#   MAT = mean(naomit$MAT),
#   gamma_rich = mean(naomit$gamma_rich),
#   HumanDisturbance = mean(naomit$HumanDisturbance),
#   N_Deposition = mean(naomit$N_Deposition),
#   MAP = mean(naomit$MAP)
# ))
# 
# predict.map <- cbind(predict.data.map,data.frame(predict(multinom.baseline1, newdata = predict.data.map, type = "probs")))
# predict.mat <- cbind(predict.data.mat,data.frame(predict(multinom.baseline1, newdata = predict.data.mat, type = "probs")))
# predict.ndep <- cbind(predict.data.ndep,data.frame(predict(multinom.baseline1, newdata = predict.data.ndep, type = "probs")))
# predict.hdisturb <- cbind(predict.data.hdisturb,data.frame(predict(multinom.baseline1, newdata = predict.data.hdisturb, type = "probs")))
# predict.gdiv <- cbind(predict.data.gamma,data.frame(predict(multinom.baseline1, newdata = predict.data.gamma, type = "probs")))
# predict.anpp <- cbind(predict.data.anpp,data.frame(predict(multinom.baseline1, newdata = predict.data.anpp, type = "probs")))
# 
# p.gather.map <- gather(predict.map, key= "codom", value = "Probability", 7:9)
# p.gather.mat <- gather(predict.mat, key= "codom", value = "Probability", 7:9)
# p.gather.ndep <- gather(predict.ndep, key= "codom", value = "Probability", 7:9)
# p.gather.hdisturb <- gather(predict.hdisturb, key= "codom", value = "Probability", 7:9)
# p.gather.gdiv <- gather(predict.gdiv, key= "codom", value = "Probability", 7:9)
# p.gather.anpp <- gather(predict.anpp, key= "codom", value = "Probability", 7:9)
# 
# grid.arrange(
# p.gather.map %>% 
#   ggplot(aes(x= MAP, y= Probability, colour = codom))+
#   geom_line()+
#   labs(y= "Probability")+
#   theme(legend.position = "none"),
# p.gather.mat %>% 
#   ggplot(aes(x= MAT, y= Probability, colour = codom))+
#   geom_line()+
#   labs(y= "Probability")+
#   theme(legend.position = "none"),
# p.gather.ndep %>% 
#   ggplot(aes(x= N_Deposition, y= Probability, colour = codom))+
#   geom_line()+
#   labs(y= "Probability")+
#   theme(legend.position = "none"),
# p.gather.hdisturb %>% 
#   ggplot(aes(x= HumanDisturbance, y= Probability, colour = codom))+
#   geom_line()+
#   labs(y= "Probability")+
#   theme(legend.position = "none"),
# p.gather.gdiv %>% 
#   ggplot(aes(x= gamma_rich, y= Probability, colour = codom))+
#   geom_line()+
#   labs(y= "Probability")+
#   theme(legend.position = "none"),
# p.gather.anpp %>% 
#   ggplot(aes(x= anpp, y= Probability, colour = codom))+
#   geom_line()+
#   labs(y= "Probability")+
#   theme(legend.position = "none")
# )
#######################################################################



# Code to generate figure of multinomial model----------------------------------------------------------

# Generate mean (later used in model predicted data)
df_om <- modeSite %>% 
  na.omit() %>% 
  mutate(MAP_mean = mean(MAP),
         MAT_mean = mean(MAT),
         gamma_rich_mean = mean(gamma_rich),
         HumanDisturbance_mean = mean(HumanDisturbance),
         N_Deposition_mean = mean(N_Deposition),
         anpp_mean = mean(anpp)) 

# Variables of interest
var <- c("MAP", "MAT", "gamma_rich", "HumanDisturbance", "N_Deposition", "anpp") 

# Generate sequence of values looped for each variable 
df_seq <- foreach(v = var, .combine = bind_cols) %do% {
  
  df_v <- df_om %>% 
    select(v)
  
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
  select(ends_with("_mean")) %>% 
  slice(1:100) %>% 
  add_column(seq = seq(from = 1, to = 100, length.out = 100)) %>% 
  full_join(df_seq, by = "seq") 

# Mean of variables of interest
var2 <- c("MAP_mean", "MAT_mean", "gamma_rich_mean", "HumanDisturbance_mean", "N_Deposition_mean", "anpp_mean") 

# Predict data using model, sequence, and mean
df_predicted <- foreach(v = var, 
                        v2 = var2,
                        .combine = bind_cols) %do% {
 # Select variable sequence
   df_v <- df_r %>% 
    select(v, seq)
 # Select means for all other variables   
   df_m <- df_r %>% 
     select(ends_with("_mean"), seq) %>% 
     select(!v2)
 # Combine data sets to predict from   
   df_s <- df_v %>% 
     full_join(df_m, by = "seq") %>% 
     select(!seq) %>% 
     rename_with(~str_remove(., '_mean'))
 # Predict data
   df_p <- cbind(df_s,
                 data.frame(predict(multinom.baseline1,
                                    newdata = df_s,
                                    type = "probs")))
 # Format for figure interpretation   
   df_c <- df_p %>% 
     select(v, starts_with("X")) %>% 
     pivot_longer(cols = starts_with("X"), 
                  names_to = "Codom", 
                  values_to = "Probability")
}

# Clarify column names 
df_combined <- df_predicted %>% 
  select(Codom...2, MAP, MAT, gamma_rich, HumanDisturbance, N_Deposition, anpp,
         starts_with("Probability")) %>% 
  rename(Codom = Codom...2,
         "Gamma Diversity" = gamma_rich,
         "Human Disturbance" = HumanDisturbance, 
         "N Deposition" = N_Deposition,
         ANPP = anpp,
         Prob_MAP = Probability...3,
         Prob_MAT = Probability...6,
         Prob_gamma = Probability...9,
         Prob_Human = Probability...12,
         Prob_N = Probability...15,
         Prob_anpp = Probability...18)

# Assign names
prob <- c("Prob_MAP", "Prob_MAT", "Prob_gamma", "Prob_anpp", "Prob_Human", "Prob_N")
named_var <- c("MAP", "MAT", "Gamma Diversity", "ANPP", "Human Disturbance", "N Deposition")

# Pre-allocate list
output <- list()

# Generate figures using loop across variables and their probabilities 
output <- foreach(v = named_var, p = prob, 
                  .combine = 'c') %do% {
    
 # Select respective variable and its probability 
  df_combo <- df_combined %>% 
    select(Codom, v, p) %>% 
    rename(v1 = v, p1 = p)
  
 # Generate figure 
  fig <- ggplot(df_combo, aes(x = v1, y = p1, color = Codom)) +
    geom_point(size = 2) + # ggMarginal must use geom_point
    labs(y = "Probability",
         x = v) +
    ylim(0.0, 0.8) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          text = element_text(size = 25),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          legend.position = "right",
          legend.text = element_text(size = 16),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.box.background = element_rect(fill = "transparent", color = NA),
          legend.key = element_rect(fill = "transparent", color = NA)) +
    scale_color_manual(name = "",
                       labels = c("Monodominated", "Codominated", "Even"),
                       values = c("#02385A", "#A63922", "#D8B573"))
  
  # Extract legend
  legend <- get_legend(fig) 
  
  # Remove legend so density can be added 
  fig2 <- fig + theme(legend.position = "none")
  
 # Add density plot on y-axis
  fig_q1 <- ggExtra::ggMarginal(fig2, type = 'density', margins = 'y',
                     size = 5, groupColour = TRUE, groupFill = TRUE)
 
  list(fig_q1)  # Return plot as list item 
} 

# generate and add legend
# get_legend <- function(myplot) {
#   tmp <- ggplotGrob(myplot)
#   gtable::gtable_filter(tmp, "guide-box")
# }

# Create standalone plot to extract legend using the first variable pair
# example_df <- df_combined %>%
#   select(Codom, all_of(named_var[1]), all_of(prob[1])) %>%
#   rename(v1 = all_of(named_var[1]), p1 = all_of(prob[1]))
# 
# legend_plot <- ggplot(example_df, aes(x = v1, y = p1, color = Codom)) +
#   geom_point() +
#   theme(legend.position = "right",
#         legend.text = element_text(size = 20),    # make labels smaller
#         legend.key.size = unit(0.4, "cm"),        # make color boxes smaller
#         legend.spacing.y = unit(0.2, "cm"),        # reduce vertical spacing
#         # Make legend background clear
#         legend.background = element_rect(fill = "transparent", color = NA),
#         legend.box.background = element_rect(fill = "transparent", color = NA),
#         legend.key = element_rect(fill = "transparent", color = NA)
#            ) +
#   # make legend dots bigger
#   guides(color = guide_legend(override.aes = list(size = 5))) +
#   scale_color_manual(name = "",
#                      labels = c("Monodominated", "Codominated", "Even"),
#                      values = c("#02385A", "#A63922", "#D8B573"))
# 
# legend <- get_legend(legend_plot)

# Overlay the legend on the top-right of the first plot
plot_with_legend <- grobTree(
  output[[1]],  # already a gtable/grob from ggMarginal
  grobTree(legend, vp = viewport(x = 0.87, y = 1.33, just = c("right", "top"))))

# Replace only the first plot with overlaid legend
output[[1]] <- plot_with_legend

# Show all plots
grid.arrange(grobs = output, nrow = 2, ncol = 3)


out_hist <- foreach(h = named_var) %do% {
  
  # Prepare dataframe
  df_o <- df_hist1 %>%
    select(variable1, mode.levels, value) %>% 
    filter(variable1 == h) %>% 
    na.omit()
  
  # Generate histogram
  fig_h <- ggplot(df_o, aes(value)) +
    geom_histogram(aes(fill = mode.levels), bins = 50) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0, 2.2, 0, 1.3), "cm"),
          plot.background = element_rect(fill = "white", colour = "white")) +
    scale_fill_manual(name = "",
                      labels = c("Monodominated", "Codominated", "Even"),
                      values = c("#02385A", "#A63922", "#D8B573")) +
    labs(x = "", y = "")
}


print(grid.arrange(out_hist[[1]], out_hist[[2]], out_hist[[3]],  
                   output[[1]], output[[2]], output[[3]], 
                   out_hist[[4]], out_hist[[5]], out_hist[[6]], 
                   output[[4]], output[[5]], output[[6]], 
                   nrow = 4, ncol = 3, heights = c(1, 2, 1, 2)))




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
       x = "Observed", y = "Predicted") +
  theme_minimal()




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

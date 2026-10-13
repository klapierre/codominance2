################################################################################
##  05b_Q1_site_modes.R: Calculate mode of codominance numbers in control plots and 
##  compare to environmental data.
##
##  Authors: Ashley LaRoque, Jordan Winter, Elise Grabda, Alyssa Young, Rachael Brenneman (modified K. Komatsu)
##  Date created: 
################################################################################

source("code/01_library.R")
source("code/02_functions.R")


# Read data --------------------------------------------------------------------

modePlot <- readRDS("data/modePlot.rds")

df_iap <- readRDS("data/modeSiteRichness.rds") %>% 
  rename(ANPP=anpp,
         GDiv=gamma_rich,
         NDep=NDeposition,
         PlotRichness=site_rich_mean) %>% 
  mutate(z_MAP=as.numeric(scale(MAP)),
         z_MAT=as.numeric(scale(MAT)),
         z_GDiv=as.numeric(scale(GDiv)),
         z_ANPP=as.numeric(scale(ANPP)),
         z_HumanFootprint=as.numeric(scale(HumanFootprint)),
         z_NDep=as.numeric(scale(NDep)),
         z_PlotRichness=as.numeric(scale(PlotRichness)))

plotRichness <- df_iap %>% 
  select(database, site_code, project_name, community_type, PlotRichness, ANPP, GDiv, NDep, 
         z_MAP, z_MAT, z_GDiv, z_ANPP, z_PlotRichness, z_HumanFootprint, z_NDep)

df_iap15 <- readRDS("data/modeSite_cutoff15.rds")%>%
  left_join(plotRichness)
df_iap25 <- readRDS("data/modeSite_cutoff25.rds")%>%
  left_join(plotRichness)
df_iap30 <- readRDS("data/modeSite_cutoff30.rds")%>%
  left_join(plotRichness)

# Colors
dominance_colors <- c("Monodominated" = "#007BA7", "Codominated" = "#A63922", "Even" = "#D8B573")
dominance_colors2 <- c("monodominated" = "#007BA7", "codominated" = "#A63922", "even" = "#D8B573")


# Figure set-up ----------------------------------------------------------------
theme_set(theme_bw())
theme_update(text = element_text(size = 35),
             axis.text.x = element_text(size = 30),
             axis.text.y = element_text(size = 30),
             legend.position = "right",
             legend.text = element_text(size = 30),
             legend.background = element_rect(fill = "transparent", color = NA),
             legend.box.background = element_rect(fill = "transparent", color = NA),
             legend.key = element_rect(fill = "transparent", color = NA),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank())


# Fig: Distribution of codominance across sites ----------------------------------------------------------

mapData <- df_iap %>% 
  mutate(codom_category = ifelse(lumpMode == 1, "Mono", 
                                 ifelse(lumpMode == 2, "Codom", "Even"))) %>% 
  filter(!is.na(NDep))

mapData$codom_category <- factor(mapData$codom_category, levels = c("Mono", "Codom", "Even"))

cc_bins <- mapData %>% 
  group_by(codom_category) %>% 
  summarise(count = n()) 
cc_binsdf <- as.data.frame(cc_bins)  

cc_binsdf %>% 
  ggplot(aes(x = factor(codom_category, levels = c("Mono", "Codom", "Even"))))+
  geom_bar(aes(y=count),stat = "identity", fill = dominance_colors)+
  geom_text(aes(y=count,label = count),vjust = -0.1, size = 20)+
  labs(x="",y= "Count")+
  theme_classic2()+
  theme(element_text(size = 20),
        axis.text.x = element_text(size = 40, angle = 0, vjust = 0.5, hjust =0.5),
        axis.title.y = element_text(size = 40),
        axis.text.y  = element_text(size = 40),
        plot.margin = margin(1,.4,.4,.4,"cm"))

# ggsave("Fig1_hist.png", width = 10, height = 10, dpi = 400)


# Fig: Correlation matrix of all the response variables ------------------------

# shows that correlations are statistically significant but weak-moderately correlated with one another
(fig_corr <- chart.Correlation(df_iap[, c("MAP", 
                                          "MAT", 
                                          # "Aridity",
                                          # "cv_Precip",
                                          "ANPP", 
                                          "GDiv", 
                                          "PlotRichness",
                                          "HumanFootprint", 
                                          "NDep"
)],
method = "pearson",
histogram = TRUE))


# Multinomial Analysis (using monodominance as reference state) ----------------

multinom.interaction <- multinom(factor(df_iap$lumpMode, levels = c(1,2,4)) ~ z_ANPP*(z_MAP + z_MAT + z_GDiv + z_PlotRichness) + z_NDep + z_HumanFootprint,
                                 data = df_iap)

multinom.15 <- multinom(factor(df_iap15$lumpMode, levels = c(1,2,4)) ~ z_ANPP*(z_MAP + z_MAT + z_GDiv + z_PlotRichness) + z_NDep + z_HumanFootprint,
                                 data = df_iap15)

multinom.25 <- multinom(factor(df_iap25$lumpMode, levels = c(1,2,4)) ~ z_ANPP*(z_MAP + z_MAT + z_GDiv + z_PlotRichness) + z_NDep + z_HumanFootprint,
                                 data = df_iap25)

multinom.30 <- multinom(factor(df_iap30$lumpMode, levels = c(1,2,4)) ~ z_ANPP*(z_MAP + z_MAT + z_GDiv + z_PlotRichness) + z_NDep + z_HumanFootprint,
                                 data = df_iap30)


AIC(multinom.interaction, multinom.15, multinom.25, multinom.30) # lowest AIC is multinom.interaction

pR2(multinom.interaction) #McFadden R2 = 0.393
pR2(multinom.15) #McFadden R2 = 0.349
pR2(multinom.25) #McFadden R2 = 0.417
pR2(multinom.30) #McFadden R2 = 0.402

# Model Output in table form 
results <- tidy(multinom.interaction, conf.int = TRUE, exponentiate = TRUE) 

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


# Multinomial model predictions ------------------------------------------------

# Interactive Responses

# Model predictors
z_vars <- c("z_ANPP", "z_MAP", "z_MAT", "z_GDiv", "z_PlotRichness", "z_NDep", "z_HumanFootprint")

# non-focal predictors held at mean values
z_means <- df_iap %>%
  summarise(across(all_of(z_vars), ~ mean(.x, na.rm = TRUE)))

# Back-transform z-scores to raw units
raw_name <- c(z_ANPP = "ANPP", z_MAP = "MAP", z_GDiv = "GDiv")

backtransform <- function(z_value, z_variable) {
  raw_variable <- raw_name[[z_variable]]
  z_value * sd(df_iap[[raw_variable]], na.rm = TRUE) + mean(df_iap[[raw_variable]], na.rm = TRUE)
}

interaction_predictions <- function(model, panel_z, panel_name) {
  
  # Low, median, and high values for the panel variable
  panel_values <- quantile(df_iap[[panel_z]],
                           probs = c(0.25, 0.50, 0.75), na.rm = TRUE, names = FALSE)
  
  panel_key <- tibble(panel_z_value = panel_values,
                      panel_level = c("Low", "Median", "High"),
                      panel_raw_value = backtransform(panel_values, panel_z)) %>%
    mutate(panel_label = paste0(panel_level, " ", panel_name, " (", round(panel_raw_value, 0), ")"))
  
  # ANPP varies continuously along the x-axis
  anpp_values <- seq(min(df_iap$z_ANPP, na.rm = TRUE),
                     max(df_iap$z_ANPP, na.rm = TRUE),
                     length.out = 200)
  
  prediction_grid <- expand_grid(z_ANPP = anpp_values,
                                 panel_z_value = panel_values)
  
  # Set all predictors to their mean, then replace ANPP and panel variable
  newdata <- as_tibble(z_means[rep(1, nrow(prediction_grid)), , drop = FALSE])
  
  newdata$z_ANPP <- prediction_grid$z_ANPP
  newdata[[panel_z]] <- prediction_grid$panel_z_value
  newdata$panel_z_value <- prediction_grid$panel_z_value
  
  probabilities <- predict(model, newdata = newdata, type = "probs") %>%
    as.data.frame()
  
  # Outcome levels: 1 = monodominated; 2 = codominated; 4 = even
  names(probabilities) <- c("1", "2", "4")
  
  bind_cols(newdata, probabilities) %>%
    left_join(panel_key, by = "panel_z_value") %>%
    mutate(ANPP_raw = backtransform(z_ANPP, "z_ANPP"),
           panel_label = factor(panel_label, levels = panel_key$panel_label[c(3, 2, 1)])) %>%
    pivot_longer(cols = c("1", "2", "4"), names_to = "dominance_code", values_to = "Probability") %>%
    mutate(Community = recode(dominance_code, "1" = "Monodominated", "2" = "Codominated", "4" = "Even"),
           Community = factor(Community, levels = c("Monodominated", "Codominated", "Even")))
}

# Predictions for interactions
pred_anpp_map <- interaction_predictions(model = multinom.interaction, panel_z = "z_MAP", panel_name = "MAP")

pred_anpp_gdiv <- interaction_predictions(model = multinom.interaction, panel_z = "z_GDiv", panel_name = "Gamma Diversity")


# ANPP v MAP scatter
guidesMAP <- pred_anpp_map %>%
    distinct(panel_label, panel_raw_value)
  
p_map_scatter <- ggplot(data=df_iap, aes(x=ANPP, y=.data[['MAP']], color=group)) +
  geom_point(alpha=0.35, size=2) +
  geom_hline(data=guidesMAP, aes(yintercept=panel_raw_value), inherit.aes=F, linetype="dashed", color="black") +
  labs(x=NULL, y='MAP') +
  scale_color_manual(values=dominance_colors2) +
  theme(legend.position='none')
  
# ANPP v MAP scatter
guidesGDiv <- pred_anpp_gdiv %>%
  distinct(panel_label, panel_raw_value)

p_gdiv_scatter <- ggplot(data=df_iap, aes(x=ANPP, y=.data[['GDiv']], color=group)) +
  geom_point(alpha=0.35, size=2) +
  geom_hline(data=guidesGDiv, aes(yintercept=panel_raw_value), inherit.aes=F, linetype="dashed", color="black") +
  labs(x=NULL, y='Gamma Diversity') +
  scale_color_manual(values=dominance_colors2) +
  theme(legend.position='none')

# ANPP:MAP probability panels
p_map_lines <- ggplot(data=pred_anpp_map, aes(x = ANPP_raw, y = Probability, color = Community, group = Community)) +
  geom_line(linewidth = 1) +
  facet_grid(panel_label ~ .) +
  scale_color_manual(values = dominance_colors) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(x = 'Aboveground Net Primary Productivity', y = 'Predicted Probability', color = NULL) +
  theme(legend.position = "top")

# ANPP:GDiv probability panels
p_gdiv_lines <- ggplot(data=pred_anpp_gdiv, aes(x = ANPP_raw, y = Probability, color = Community, group = Community)) +
  geom_line(linewidth = 1) +
  facet_grid(panel_label ~ .) +
  scale_color_manual(values = dominance_colors) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(x = 'Aboveground Net Primary Productivity', y = 'Predicted Probability', color = NULL) +
  theme(legend.position = "top")


# Combine panels into figure
map_block <- p_map_scatter / p_map_lines +
  plot_layout(heights = c(0.8, 3))

gdiv_block <- p_gdiv_scatter / p_gdiv_lines +
  plot_layout(heights = c(0.8, 3))

(map_block | gdiv_block) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

# ggsave("Fig2_model_4.0.png", width = 20, height = 25, dpi = 400)



# Mean Responses ---------------------------------------------------------------

z_vars <- c("z_ANPP", "z_MAP", "z_MAT", "z_GDiv", "z_PlotRichness", "z_NDep", "z_HumanFootprint")
raw_vars <- c("ANPP", "MAP", "MAT", "GDiv", "PlotRichness", "NDep", "HumanFootprint")

model_df <- df_iap %>%
  filter(if_all(all_of(c("lumpMode", z_vars, raw_vars)), ~ !is.na(.x)))

z_reference <- model_df %>%
  summarise(across(all_of(z_vars), ~ mean(.x, na.rm = TRUE)))

raw_to_z <- function(raw_values, raw_variable) {
  (raw_values - mean(df_iap[[raw_variable]], na.rm = TRUE)) / sd(df_iap[[raw_variable]], na.rm = TRUE)
}

main_effect_predictions <- function(raw_variable, z_variable, figure_label, model = multinom.interaction) {
  
  # Raw predictor values for the x-axis
  raw_sequence <- seq(
    min(model_df[[raw_variable]], na.rm = TRUE),
    max(model_df[[raw_variable]], na.rm = TRUE),
    length.out = 200)
  
  # predictors at means
  newdata <- as.data.frame(
    z_reference[rep(1, length(raw_sequence)), , drop = FALSE]
  )
  
  # replace z with raw data
  newdata[[z_variable]] <- raw_to_z(raw_sequence, raw_variable)
  
  probabilities <- as.data.frame(predict(model, newdata = newdata, type = "probs"))
  
  # fix names
  names(probabilities) <- c("1", "2", "4")
  
  bind_cols(tibble(x_raw = raw_sequence, Predictor = figure_label), probabilities) %>%
    pivot_longer(cols = all_of(c("1", "2", "4")), names_to = "dominance_code", values_to = "Probability") %>%
    mutate(Community = recode(dominance_code, "1" = "Monodominated", "2" = "Codominated", "4" = "Even"),
           Community = factor(Community, levels = c("Monodominated", "Codominated", "Even")))
}

pred_main_effects <- bind_rows(main_effect_predictions(raw_variable = "MAT",
                                                       z_variable = "z_MAT",
                                                       figure_label = "Mean Annual Temperature"),
                               main_effect_predictions(raw_variable = "PlotRichness",
                                                       z_variable = "z_PlotRichness",figure_label = "Alpha Diversity"),
                               main_effect_predictions(raw_variable = "HumanFootprint",
                                                       z_variable = "z_HumanFootprint",
                                                       figure_label = "Human Footprint Index"),
                               main_effect_predictions(raw_variable = "NDep",
                                                       z_variable = "z_NDep",
                                                       figure_label = "Nitrogen Deposition"))

support_data <- model_df %>%
  mutate(Community = recode(as.character(lumpMode), "1" = "Monodominated", "2" = "Codominated", "4" = "Even"),
         Community = factor(Community, levels = c("Monodominated", "Codominated", "Even")))


# Figure of main effects -------------------------------------------------------
main_effect_block <- function(raw_variable, figure_label, bins = 30) {
  
  x_limits <- range(model_df[[raw_variable]], na.rm = TRUE)
  
  # Histograms 
  p_hist <- ggplot(support_data, aes(x = .data[[raw_variable]], fill = Community)) +
    geom_histogram(bins = bins, position = "stack", alpha = 0.75) +
    scale_fill_manual(values = dominance_colors, guide = "none") +
    coord_cartesian(xlim = x_limits) +
    ylab('Count') +
    theme(axis.title.x=element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  # Predicted probabilities
  p_lines <- pred_main_effects %>%
    filter(Predictor == figure_label) %>%
    ggplot(aes(x = x_raw, y = Probability, color = Community, group = Community)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = dominance_colors) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    coord_cartesian(xlim = x_limits, ylim = c(0, 1)) +
    labs(x = figure_label, y = "Predicted Probability", color = NULL) +
    theme(legend.position = "top")
  
  p_hist / p_lines +
    plot_layout(heights = c(0.3, 1))
}

# Make figures and combine
p_MAT <- main_effect_block(raw_variable = "MAT", figure_label = "Mean Annual Temperature")
p_PlotRichness <- main_effect_block(raw_variable = "PlotRichness", figure_label = "Alpha Diversity")
p_HumanFootprint <- main_effect_block(raw_variable = "HumanFootprint", figure_label = "Human Footprint Index")
p_NDep <- main_effect_block(raw_variable = "NDep", figure_label = "Nitrogen Deposition")

wrap_plots(p_MAT, p_PlotRichness,
           p_HumanFootprint, p_NDep,
           ncol = 2,
           guides = "collect") &
  theme(legend.position = "top")

# ggsave("FigS5_modelExtra.png", width = 20, height = 20, dpi = 400)



ggplot(
  pred_main_effects,
  aes(
    x = x_raw,
    y = Probability,
    color = Community,
    group = Community
  )
) +
  geom_line(linewidth = 1) +
  facet_wrap(
    ~ Predictor,
    ncol = 2,
    scales = "free_x"
  ) +
  scale_color_manual(values = dominance_colors) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  labs(
    x = "Predictor value",
    y = "Predicted probability",
    color = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold")
  )


# saveRDS(df_combined, file = "data/multimodalModelPredictions.rds") # saving derived data for analyses

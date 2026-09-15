
################################################################################
##  05d_Q1_continuousDominance.R: Calculate dominance as integers in control plots and 
##  compare to environmental data.
##
##  Authors: Kimberly Komatsu
##  Date created: Sept 15, 2026
################################################################################

source("code/01_library.R")
source("code/02_functions.R")


# Read data ---------------------------------------------------------------

numCodomPlotYear <- readRDS("data/numCodomPlotYear.rds") %>% 
  separate(exp_unit, into=c('site_code', 'project_name', 'community_type', 
                            'plot_id', 'treatment', 'calendar_year'), 
           sep='::', remove=F) %>% 
  left_join(readRDS("data/expInfo.rds")) %>% 
  left_join(readRDS("data/envData.rds")) %>% 
  filter(trt_type=='control')


# # Determine structure of data within each experiment
# projects <- unique(numCodomPlotYear$project_name)
# 
# for (project in projects) {
#   project_data <- subset(numCodomPlotYear, project_name == project)
#   
#   p <- ggplot(project_data, aes(x = num_codominants_continuous)) +
#     geom_bar() +
#     facet_wrap(~ calendar_year) +
#     theme_minimal() +
#     labs(title = project)
#   
#   ggsave(filename = paste0("figs/histogram_", project, ".png"),
#          plot = p,
#          width = 10,
#          height = 4,
#          dpi = 300)
# }


# Mixed-Effects Model ----------------------------------------------------

model <- glmmTMB(
  num_codominants_continuous ~ scale(anpp)*(scale(MAP) + scale(MAT) + scale(gamma_rich)) + scale(NDeposition) + scale(HumanFootprint) +
    (1 | site_code / project_name / community_type / plot_id),
  family = nbinom2,
  data = numCodomPlotYear
)

summary(model)

# ANPP * MAP figure
pred_MAP <- ggpredict(model, terms = c("MAP [all]", "anpp [quartiles]"))

ggplot(pred_MAP, aes(x = x, y = predicted, color=group, fill=group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.2
  ) +
  labs(
    x = "MAP",
    y = "Predicted number of dominant species",
    color = 'ANPP', fill = 'ANPP'
  ) +
  theme_minimal()


# MAT figure

pred_MAT <- ggpredict(model, terms = "MAT [all]")

ggplot(pred_MAT, aes(x = x, y = predicted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.2
  ) +
  labs(
    x = "MAT",
    y = "Predicted number of dominant species"
  ) +
  theme_minimal()


# gamma div figure

pred_gammaDiv <- ggpredict(model, terms = "gamma_rich [all]")

ggplot(pred_gammaDiv, aes(x = x, y = predicted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.2
  ) +
  labs(
    x = "Gamma Diversity",
    y = "Predicted number of dominant species"
  ) +
  theme_minimal()

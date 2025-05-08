
rm(list = ls(all = TRUE))
# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(cowplot)

# PATHS -------------------------------------------------------------------

BASE = here::here()
BASE = "/Users/hanleyb/Documents/GitHub/VAULT_Breast"
OUT_DIR = file.path(BASE, "analysis", "figures")
FLOW_PATH = file.path(BASE, "data", "fixed_flow", "fixed_flow.txt")
KI67_PATH = file.path(BASE, "data", "image_analysis", "Ki67_pathologist_scores.txt")
CLINDATA = file.path(BASE, "data", "metadata", "clinical_data.txt")


# LOAD DATA ---------------------------------------------------------------
flow = read.delim(FLOW_PATH)
ki67 = read.delim(KI67_PATH)
clinical_data = read.delim(CLINDATA)

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_features.R"))



# GRAPHICAL PARAMETERS ----------------------------------------------------
point_size = 5
base_size = 20
font_size = 8


# ORGANISE DATA -----------------------------------------------------------
flow$Trial_ID = flow$Trial_id
mat_out = left_join(ki67, flow, by = "Trial_ID")
clinical_data = clinical_data%>%select(everything(), Trial_ID = Trial.ID)
flow_out = left_join(flow, clinical_data, by = "Trial_ID")



# DRAW PLOTS --------------------------------------------------------------
plot_1 = 
mat_out%>%
  group_by(Trial_ID)%>%
  mutate(distinct = n_distinct(filenames))%>%
  group_by(filenames)%>%
  mutate(order = Surface_Area*exp_perc,
         region_mean = sum(order))%>%
  group_by(Trial_ID)%>%
  mutate(max = max(region_mean),
         min = min(region_mean),
         mean = mean(region_mean),
         range = max-min)%>%
  ggplot(aes(mean, KI67.ALL, label = Trial_ID))+
  #geom_text()+
  geom_point(size = point_size)+
  geom_smooth(method = "lm", col = "black")+
  stat_cor(size = font_size)+
  theme_classic(base_size = base_size)+
  ylab("Ki67 (%) Fixed FACS")+
  xlab("Ki67 (%) IHC - Tumour Mean")



plot_2 =
flow_out%>%
  filter(!is.na(Pathology.Grading.Components_M))%>%
  ggplot(aes(as.character(Pathology.Grading.Components_M), KI67.ALL))+
  #geom_boxplot(outlier.alpha = 0)+
  #geom_jitter(size = point_size)+
  geom_violin(fill = "lightgrey", width= 1.2)+
  geom_boxplot(width = 0.5, outlier.alpha = 0)+
  geom_jitter(size = point_size/2, width = 0.11)+
  theme_classic(base_size = base_size)+
  xlab("Mitotic Rate - H&E")+
  ylab("Ki67 (%) Fixed FACS")+
  stat_compare_means(comparisons = list(c("1", "2"),
                                        c("1", "3"),
                                        c("3", "2")),label = "p.signif", size = font_size)+
  scale_y_log10()


plot_3 = 
flow_out%>%
  mutate(Recurrence.score_GEP = if_else(Recurrence.score_GEP == 0, NA, Recurrence.score_GEP))%>%
  mutate(label = if_else(Trial_ID == "HF007", "HF007", ""))%>%
  filter(!is.na(Recurrence.score_GEP))%>%
  ggplot(aes(KI67.ALL, Recurrence.score_GEP, label = label))+
  geom_point(size = point_size, aes(col = GEP_performed))+
  geom_smooth(method = "lm", col = "black", se = F)+
  theme_classic(base_size = base_size)+
  geom_text(nudge_x = 5, fontface = "bold", size = font_size*0.75)+
  stat_cor(size = font_size)+
  scale_color_manual(values = colours[3:4])+
  theme(legend.title = element_blank(),
        legend.position = c(0.8,0.1),
        legend.text = element_text(size = font_size*2))+
  xlab("Ki67 (%) Fixed FACS")+
  ylab("Recurrence Score")

plot_grid(plot_1, plot_2, plot_3, ncol = 3, labels = LETTERS[1:3], label_size = 25)


# WRITE DATA --------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Extended_Data_Figure10_correlations_w_proliferation.png"))
       
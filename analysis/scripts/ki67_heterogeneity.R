# Create a figure showing the selected genes by CCF

rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------
library(tidyverse)


# PATHS -------------------------------------------------------------------
BASE = "/Users/hanleyb/Dropbox (The Francis Crick)/HoLSTF_Breast/Github_Repo"
BASE = here::here()
OUT_DIR = file.path(BASE, "analysis", "figures")
KI67_LOC = file.path(BASE, "data", "image_analysis", "Ki67_pathologist_scores.txt")


# LOAD DATA ---------------------------------------------------------------
mat_out = read.delim(KI67_LOC)


# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "rot.lab.R"))

# GRAPHICAL PARAMETERS ----------------------------------------------------
sr_col = "#33A02Cb7"
rs_col = "#6A3D9Aff"


# DRAW PLOT ---------------------------------------------------------------
# get x-axis order
order = 
  mat_out%>%
  group_by(filenames)%>%
  mutate(order = Surface_Area*exp_perc,
         order = sum(order))%>%
  arrange(Trial_ID, order)%>%
  filter(!duplicated(filenames))%>%
  ungroup()%>%pull(filenames)

mat_out%>%  
  group_by(filenames)%>%
  mutate(order = Surface_Area*exp_perc,
         region_mean = sum(order))%>%
  group_by(Trial_ID)%>%
  mutate(max = max(region_mean),
         min = min(region_mean),
         mean = mean(region_mean),
         range = max-min)%>%
  ggplot()+
  geom_jitter(aes(reorder(Trial_ID, mean), exp_perc*100, col = "sub-regional"), alpha = .5, size = 2)+
  geom_point(aes(Trial_ID, region_mean, col = "regional_mean"), alpha = 0.6, size = 3)+
  geom_point(aes(Trial_ID, mean, col = "tumour_mean"), alpha = 0.9, size = 4)+
  scale_color_manual(values = c("tumour_mean" = rs_col, "regional_mean" = sr_col, "sub-regional" = "grey"))+
  theme_classic(base_size = 20)+
  ylab("Ki67 Positivity (%)")+
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.15,0.8))+
  rot.lab()



# WRITE PLOT --------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Figure5A_Ki67_heterogeneity.png"))



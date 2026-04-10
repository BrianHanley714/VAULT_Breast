# draw extended data plots for time in formalin and time for repsamp

# Consort Diagram

rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(cowplot)


# LOAD DATA ---------------------------------------------------------------
BASE = here::here()
BASE = "/Users/hanleyb/Documents/GitHub/VAULT_Breast"
OUT_DIR = file.path(BASE, "analysis", "figures")
CLINDATA = file.path(BASE, "data", "metadata", "clinical_data.txt")
ENDPOINT = file.path(BASE, "data", "metadata", "VAULT_endpoint_data.txt")
TUM_COUNTS = file.path(BASE, "data", "image_analysis", "tumour_cell_counts.tsv")


# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_features.R"))
get_lighter = function(col){colorRampPalette(c("white", col))(3)[2]} # blend colours with white

# LOAD DATA ---------------------------------------------------------------
clinical_data = read.delim(CLINDATA)
endpoint_data = read.delim(ENDPOINT)
annotations = read.delim(TUM_COUNTS)
# get clinical data
clinical_data = left_join(clinical_data, endpoint_data, by = "Trial.ID")


# GRAPHICAL PARAMETERS ----------------------------------------------------
rs_col = "#54278F"
sr_col = "#33A02Cb7"
dark_grey  <- "#4D4D4D"
light_grey <- "#D9D9D9"
base_size = 30
labels_x = 1.5 #control location of RCPath equivalent workload points on the x-axis

# Blend colours with white to lighten
rs_col = get_lighter(rs_col)
sr_col = get_lighter(sr_col)
col_palette <- colorRampPalette(c(rs_col, sr_col))
col_blends <- col_palette(5)



# DRAW TIME IN FORMALIN PLOT ----------------------------------------------
plot_time_formalin = clinical_data%>%
  dplyr::select(Time..Buffer.to.Homogenisation..days., Time..Surgery.to.Homogenisation..days., Time..in.formalin..days., Time..Surgery.to.Grossing..days., Trial.ID)%>%
  pivot_longer(values_to = "days", cols = c(Time..Buffer.to.Homogenisation..days., Time..in.formalin..days.))%>%
  dplyr::filter(!is.na(days))%>%
  group_by(Trial.ID)%>%
  mutate(Time_to_Homogenisation = sum(days))%>%
  ungroup()%>%
  mutate(name = if_else(name == "Time..Buffer.to.Homogenisation..days.", "PBS", "Formalin"),
         Solution = factor(name, levels = c("PBS", "Formalin")))%>%
  ggplot()+
  geom_col(aes(0-days, reorder(Trial.ID, -days), fill = Solution))+
  geom_vline(xintercept = -1, col = "black")+
  geom_vline(xintercept = -3, col = "black")+
  geom_tile(aes(0-Time_to_Homogenisation,  reorder(Trial.ID, days)), fill = "#223E8B", width =1, height = 1)+
  
  geom_tile(aes(0-Time..Surgery.to.Homogenisation..days., reorder(Trial.ID, days)), fill = "#FFDB58", width =1, height = 1)+

  scale_x_continuous(position = "bottom",
                     labels = function(x) abs(x)) +
  scale_y_discrete(position = "right") +
  geom_tile(aes(0-Time..Surgery.to.Grossing..days., reorder(Trial.ID, days)), fill = "#30D5C8", width =1, height = 1)+
  ylab("Trial ID")+
  xlab("Time (days)")+
  theme_classic(base_size = base_size)+
  theme(text = element_text(size = 30), 
        legend.text = element_text(size = 20),
        legend.title = element_text(size  = 20),
        axis.text.y = element_text(size = 13),
        
        legend.position = c(0.75, 0.05))+
  guides(fill=guide_legend(nrow=1,byrow=TRUE, title = element_blank()),
   #      shape=guide_legend(nrow=1,byrow=TRUE, title = element_blank(), override.aes = list(size = 5))
  )+
  scale_fill_manual(values = c("PBS" = light_grey, "Formalin" = dark_grey))+
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background  = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA)
  )


plot_time_formalin



# WRITE PLOT --------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Figure2a_Time_in_formalin.png"), height = 12, width = 10)


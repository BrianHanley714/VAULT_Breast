# draw extended data plots for time in formalin and time for repsamp

# Consort Diagram

rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(cowplot)


# LOAD DATA ---------------------------------------------------------------
BASE = here::here()
BASE = "/Users/hanleyb/Documents/GitHub/VAULT_Breast/"
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
  filter(!is.na(days))%>%
  group_by(Trial.ID)%>%
  mutate(Time_to_Homogenisation = sum(days))%>%
  ungroup()%>%
  mutate(name = if_else(name == "Time..Buffer.to.Homogenisation..days.", "PBS", "Formalin"),
         Solution = factor(name, levels = c("PBS", "Formalin")))%>%
  ggplot()+
  geom_col(aes(days, reorder(Trial.ID, days), fill = Solution))+
  geom_point(aes(Time_to_Homogenisation,  reorder(Trial.ID, days), shape = "Homogenisation"))+
  
  geom_point(aes(Time..Surgery.to.Homogenisation..days., reorder(Trial.ID, days), shape = "Dissection"))+
  geom_vline(xintercept = 1, col = "red")+
  geom_vline(xintercept = 3, col = "red")+
  geom_point(aes(Time..Surgery.to.Grossing..days., reorder(Trial.ID, days), shape = "Routine Grossing"))+
  ylab("Trial ID")+
  xlab("Time (days)")+
  theme_classic(base_size = base_size)+
  theme(text = element_text(size = 30), 
        legend.text = element_text(size = 20),
        legend.title = element_text(size  = 20),
        axis.text.y = element_text(size = 10),

        legend.position = c(0.7, 0.3))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE),
         shape=guide_legend(nrow=3,byrow=TRUE, override.aes = list(size = 5))
  )+
  scale_fill_manual(values = c("PBS" = rs_col, "Formalin" = sr_col))



# DRAW TIME FOR REPSAMP PLOT ----------------------------------------------
plot_time_repseq = clinical_data%>%
  dplyr::select(Time.for.homogenisation..minutes., Time.for.dissection..minutes., Trial.ID)%>%
  pivot_longer(cols = c(Time.for.homogenisation..minutes., Time.for.dissection..minutes.), names_to = "Task", values_to = "Time")%>%
  mutate(Task = if_else(Task == "Time.for.homogenisation..minutes.", "Homogenisation", "Dissection"))%>%
  ggplot(aes(Task, Time))+
  geom_rect(xmin = 0, xmax= 4, ymin = 0, ymax = 5, fill = col_blends[1])+
  geom_rect(xmin = 0, xmax= 4, ymin = 5, ymax = 10, fill =  col_blends[2])+
  geom_rect(xmin = 0, xmax= 4, ymin = 10, ymax = 20, fill =  col_blends[3])+
  geom_rect(xmin = 0, xmax= 4, ymin = 20, ymax = 30, fill  = col_blends[4])+
  geom_rect(xmin = 0, xmax= 4, ymin = 30, ymax = 55, fill = col_blends[5])+
  geom_boxplot(width = 0.5, outlier.alpha = 0)+
  geom_jitter(width = 0.25)+
  ylab("Time (Minutes)")+
  stat_compare_means(label = "p.signif", label.x = 1.475, label.y = 50, fontface = "bold", size = 6)+
  geom_label(aes(labels_x, 45, label = "RCPath\nWorkload Points\nEquivalent"), size = 5, fontface = "bold")+
  geom_label(aes(labels_x, 7.5, label = "Two\nPoints"), size = 5)+
  geom_label(aes(labels_x, 15, label = "Three\nPoints"), size = 5)+
  geom_label(aes(labels_x, 25, label = "Five\nPoints"), size = 5)+
  geom_label(aes(labels_x, 35, label = "Eight\nPoints"), size = 5)+
  theme_classic(base_size = base_size)



# PLOT GRID ---------------------------------------------------------------
plot_grid(plot_time_formalin, plot_time_repseq, nrow = 1, labels = LETTERS[1:2],
          label_size = 25)


# WRITE PLOT --------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Extended_Data_Figure2_Secondary_Endpoints.png"))

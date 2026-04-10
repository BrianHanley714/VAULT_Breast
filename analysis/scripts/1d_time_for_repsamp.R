# draw extended data plots for time in formalin and time for repsamp

# Consort Diagram

rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(cowplot)


# LOAD DATA ---------------------------------------------------------------
BASE = here::here()
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
labels_x = 2.5 #control location of RCPath equivalent workload points on the x-axis

# Blend colours with white to lighten
rs_col = get_lighter(rs_col)
sr_col = get_lighter(sr_col)
col_palette <- colorRampPalette(c(rs_col, sr_col))
col_blends <- col_palette(5)





# DRAW TIME FOR REPSAMP PLOT ----------------------------------------------
plot_time_repseq = clinical_data%>%
  dplyr::select(Time.for.homogenisation..minutes., Time.for.dissection..minutes., Trial.ID)%>%
  pivot_longer(cols = c(Time.for.homogenisation..minutes., Time.for.dissection..minutes.), names_to = "Task", values_to = "Time")%>%
  mutate(Task = if_else(Task == "Time.for.homogenisation..minutes.", "Homogenisation", "Dissection"))%>%
  ggplot(aes(Task, Time))+
  geom_hline(yintercept = 5, linetype = "dashed")+
  geom_hline(yintercept = 10, linetype = "dashed")+
  geom_hline(yintercept = 20, linetype = "dashed")+
  geom_hline(yintercept = 30, linetype = "dashed")+
  geom_violin(width = 0.8, aes(fill = Task))+
  geom_boxplot(width = 0.2, outlier.alpha = 0)+
  geom_jitter(width = 0.05, size = 1)+
  ylab("Time (Minutes)")+
  scale_fill_manual(values =c(rs_col, sr_col))+
  stat_compare_means(label = "p.signif", label.x = 1.475, label.y.npc = 0.9, fontface = "bold", size = 10)+
  geom_point(aes(2.8, 1), alpha = 0)+
  #geom_label(aes(labels_x, 45, label = "RCPath\nWorkload Points\nEquivalent"), size = 5, fontface = "bold")+
  geom_label(aes(labels_x, 7.5, label = "Two\nPoints"), size = 5)+
  geom_label(aes(labels_x, 15, label = "Three\nPoints"), size = 5)+
  geom_label(aes(labels_x, 25, label = "Five\nPoints"), size = 5)+
  geom_label(aes(labels_x, 35, label = "Eight\nPoints"), size = 5)+
  theme_classic(base_size = base_size)+
  theme(legend.position = "none")


plot_time_repseq

# Statistical Tests -------------------------------------------------------



# dissection sigificantly less time than homogenisation
clinical_data%>%
  dplyr::select(Time.for.homogenisation..minutes., Time.for.dissection..minutes., Trial.ID)%>%
  pivot_longer(cols = c(Time.for.homogenisation..minutes., Time.for.dissection..minutes.), names_to = "Task", values_to = "Time")%>%
  mutate(Task = if_else(Task == "Time.for.homogenisation..minutes.", "Homogenisation", "Dissection"))%>%
  wilcox.test(Time ~ Task, data = ., conf.int = T, alternative = "two.sided")



# WRITE PLOT --------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Figure1g_time_for_RepSamp.png"))

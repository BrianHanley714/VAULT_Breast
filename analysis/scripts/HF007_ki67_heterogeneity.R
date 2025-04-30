# Demonstrate the regional heterogeneity in HF007

rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(ggpubr)

# PATHS -------------------------------------------------------------------
BASE = here::here()
BASE = "/Users/hanleyb/Dropbox (The Francis Crick)/HoLSTF_Breast/Github_Repo"
OUT_DIR = file.path(BASE, "analysis", "figures")
HF007_LOC = file.path(BASE, "data", "image_analysis", "HF007_image_analysis.txt")


# LOAD DATA ---------------------------------------------------------------
HF007 = read.delim(HF007_LOC)


# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "rot.lab.R"))

# GRAPHICAL PARAMETERS ----------------------------------------------------
sr_col = "#33A02Cb7"
rs_col = "#6A3D9Aff"


# DRAW PLOT ---------------------------------------------------------------

HF007%>%
  
  filter(Name == "Tumor")%>%
  mutate(OncotypeDX_performed = if_else(block == "BA015", TRUE, FALSE),
         slice = case_when(block == "BA010"~"Slice_7",
                           block == "BA012"~"Slice_8",
                           block == "BA014"~"Slice_10",
                           block == "BA015"~"Slice_11",
                           block == "BA016"~"Slice_12",
                           block == "BA018"~"Slice_13",
         ),
         order = as.numeric(sub("BA0", "", block)))%>%
  ggplot(aes(reorder(slice, order), Tumour_cells..Positive..))+
  geom_boxplot(outlier.alpha = 0, aes(fill = OncotypeDX_performed))+
  geom_jitter(size = 3)+
  theme_classic(base_size = 20)+
  stat_compare_means(label.x = 3, label.y = 40, label = "p.signif", size = 10)+
  scale_fill_manual(values = c("TRUE" = sr_col, "FALSE" = rs_col))+
  ylab("Ki67+ Cancer Cells (%)")+
  theme(axis.title.x = element_blank(), 
        legend.position = "top")


# WRITE PLOT --------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Figure5B_HF007_Ki67_heterogeneity.png"))






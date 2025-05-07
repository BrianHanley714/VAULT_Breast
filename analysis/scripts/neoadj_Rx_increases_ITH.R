# boxplot showing ITH index increases with neoadj treatment

rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(ggpubr)

# PATHS -------------------------------------------------------------------
BASE = here::here()
BASE = "/Users/hanleyb/Documents/GitHub/VAULT_Breast/"
OUT_DIR = file.path(BASE, "analysis", "figures")
VARIANTS_VAULT = file.path(BASE, "data","variants", "variant_calls_VAULT.txt")
INCLUDED_PATIENTS = file.path(BASE, "data","metadata", "cases_included.xlsx")
CLINDATA = file.path(BASE, "data", "metadata", "clinical_data.txt")


# LOAD DATA ---------------------------------------------------------------

rs_patients = read.delim(INCLUDED_PATIENTS)[,1]
vault = read.delim(VARIANTS_VAULT)
clinical_data = read.delim(CLINDATA)
source(file.path(BASE, "src", "custom_filters.R"))



# DRAW PLOT AND CALCULATE ITH INDEX ---------------------------------------
vault%>%
  mutate(study = "VAULT")%>%
  filter_rs()%>%
  group_by(Tumor_Sample_Barcode, clonality)%>%
  count()%>%
  filter(grepl("CLONAL", clonality))%>%
  pivot_wider(values_from = n, names_from = clonality)%>%
  mutate(ITH_Index = SUBCLONAL/CLONAL)%>%
  mutate(Trial.ID = substr(Tumor_Sample_Barcode, 1, 5))%>%
  left_join(clinical_data, by = "Trial.ID")%>%
  mutate(neoadj_RX = if_else(is.na(neo_adj_treatment_class), "naive", "treated"))%>%
  ggplot(aes(neoadj_RX, ITH_Index))+
  geom_boxplot(outlier.alpha = 0)+
  geom_jitter(size = 5)+
  stat_compare_means(label = "p.signif", size = 10, label.x.npc = 0.5)+
  scale_y_log10()+
  theme_classic(base_size = 30)



# WRITE PLOT --------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Extended_Data_Figure5_neoadjuvantRX_increases_ITH.png"))

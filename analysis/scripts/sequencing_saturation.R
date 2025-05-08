
rm(list = ls(all = TRUE))
# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
# PATHS -------------------------------------------------------------------

BASE = here::here()
BASE = "/Users/hanleyb/Documents/GitHub/VAULT_Breast/"
OUT_DIR = file.path(BASE, "analysis", "figures")
BAM_PATH = file.path(BASE, "data", "metadata", "bam_qc.txt")

# LOAD DATA ---------------------------------------------------------------
bamqc_data = read.delim(BAM_PATH)



# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_features.R"))


# ANNOTATE PANEL CASES ----------------------------------------------------
# These are the cases where the custom panel was performed. This is to differentiate them from cases where whole exome sequencing was performed. 
panel_cases = c(
  "HF016_ET001_D02",
  "HF036_ET001_D02",
  "HF082_ET001_D02",
  "HF105_ET001_D02",
  "HF132_ET001_D02",
  "HF170_ET001_D02",
  "HF195_ET001_D02",
  "HF199_ET001_D02",
  "HF170_BB001_D01",
  "HF170_BA006_D01",
  "HF170_BA009_D01",
  "HF170_BA019_D01",
  "HF170_BA007_D01",
  "HF057_BB008_D01",
  "HF036_BB017_D01",
  "HF105_BB001_D01"
  )


# DRAW PLOT ---------------------------------------------------------------
bamqc_data%>%
  mutate(panel = if_else(Tumor_Sample_Barcode %in% panel_cases, "panel", "wes"))%>%
  group_by(panel, Tumor_Sample_Barcode, Recalibrated)%>%
  reframe(sum = sum(mean_coverage, na.rm = T))%>%
  pivot_wider(names_from = Recalibrated, values_from = sum)%>%
  ggplot(aes(No, Yes, col = panel, label = Tumor_Sample_Barcode))+
  geom_point(size = 5)+
  scale_x_log10()+
  scale_y_log10(limits = c(1,80000))+
  geom_smooth(method = "gam", col = "black")+
  scale_color_manual(values = colours[3:4])+
  xlab("log10(Total Coverage)")+
  ylab("log10(De-duplicated Coverage)")+
  theme_classic(base_size = 30)


# WRITE PLOT --------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Extended_Data_Figure9_saturation_of_highdepthseq.png"))


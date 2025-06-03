# Create a plot comparing the known and predicted drivers in TCGA and VAULT

rm(list = ls(all = TRUE))
# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(ggpubr)
# PATHS -------------------------------------------------------------------

BASE = here::here()
OUT_DIR = file.path(BASE, "analysis", "figures")
VARIANTS_VAULT = file.path(BASE, "data","variants", "variant_calls_VAULT.txt")
VARIANTS_MBTCGA = file.path(BASE, "data","variants", "variant_calls_TCGA_MB.txt")
MATCHED_PT_CHAR_MBTCGA = file.path(BASE, "data","metadata", "matched_patients_characteristics.txt")
IDMAP = file.path(BASE, "data","metadata", "tumouridmap_MB.txt")
INCLUDED_PATIENTS = file.path(BASE, "data","metadata", "cases_included.xlsx")
CLINDATA = file.path(BASE, "data", "metadata", "clinical_data.txt")


# LOAD DATA ---------------------------------------------------------------

rs_patients = read.delim(INCLUDED_PATIENTS)[,1]
vault = read.delim(VARIANTS_VAULT)
mbtcga = read.delim(VARIANTS_MBTCGA)
matched_patients_char = read.delim(MATCHED_PT_CHAR_MBTCGA)
tumour_IDs = read.delim(IDMAP)
matched_patients = c(matched_patients_char$PATIENT_ID, tumour_IDs$sample[match(matched_patients_char$PATIENT_ID, tumour_IDs$metabricId)])
matched_patients = matched_patients[!is.na(matched_patients)]
clinical_data = read.delim(CLINDATA)

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "custom_filters.R"))
source(file.path(BASE, "src", "rot.lab.R"))
source(file.path(BASE, "src", "annotate_driver_summary.R"))

rs_col = "#54278F"
sr_col = "#33A02Cb7"
sr_col_1 = "#33A02Cb7"
sr_col_2 = "#B2DF8Aff"
gt_col = "grey"
clonal_col = "#33A02Cb7"
subclonal_col= "#54278F"



# COMBINED DATAFRAME ------------------------------------------------------
vault$study = "VAULT"
common_names = Reduce(intersect,list(names(mbtcga), names(vault)))


combined_df = bind_rows(mbtcga%>%dplyr::select(common_names),
                        vault%>%dplyr::select(common_names)
)


# DRAW PLOT ---------------------------------------------------------------
combined_df%>%
  annotate_driver_summary()%>%  
  filter_rs()%>%
  filter(study != "METABRIC")%>%
  filter_matched()%>%
  mutate(Consensus_annotated = if_else(trial_specific =="trial_inclusion_criterion" & 
                                         grepl("pred", CGI.Oncogenic.Summary) & 
                                         Consensus_annotated == "Non_Driver", "pDriver;CTI" , 
                                       if_else(trial_specific !="trial_inclusion_criterion" & 
                                                 grepl("pred", CGI.Oncogenic.Summary) & 
                                                 Consensus_annotated == "Non_Driver", "pDriver;other", Consensus_annotated)))%>%
  mutate(Consensus_annotated = if_else(grepl("Biomarker|Driver_Other", Consensus_annotated), "kDriver", Consensus_annotated))%>%
  mutate(Consensus_annotated = factor(Consensus_annotated, levels = c("kDriver", "pDriver;CTI", "pDriver;other", "Non_Driver")))%>%
  ggplot(aes(Consensus_annotated, ccf_expected_copies))+
  geom_boxplot(aes(fill = study))+
  #geom_jitter()+
  stat_compare_means(comparisons = list(c("kDriver", "pDriver;CTI"),
                                        c("kDriver", "pDriver;other"),
                                        c("kDriver", "Non_Driver")
  ), label = "p.signif")+
  facet_wrap(~study, ncol = 1)+
  rot.lab()+
  theme_classic(base_size = 20)+
  ylim(0, 1.3)+
  ylab("Cancer Cell Fraction")+
  scale_fill_manual(values = c("TCGA" = sr_col, "VAULT" = rs_col))+
  rot.lab()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15),
        legend.position = "none")


# WRITE PLOT --------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Figure3H_known_vs_predicted_drivers.png"))


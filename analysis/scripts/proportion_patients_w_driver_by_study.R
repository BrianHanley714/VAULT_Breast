# Create a figure showing the proportion of patients with a driver by study and clonality

rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------
library(tidyverse)

# PATHS -------------------------------------------------------------------

BASE = here::here()
OUT_DIR = file.path(BASE, "analysis", "figures")
VARIANTS_VAULT = file.path(BASE, "data","variants", "variant_calls_VAULT.txt")
VARIANTS_MBTCGA = file.path(BASE, "data","variants", "variant_calls_TCGA_MB.txt")
MATCHED_PT_MBTCGA = file.path(BASE, "data","metadata", "matched_patients.txt")
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
matched_patients = c(matched_patients_char$PATIENT_ID, unique(tumour_IDs$sample[match(matched_patients_char$PATIENT_ID, tumour_IDs$metabricId)]))
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



# RUN CHI-SQUARE TEST -----------------------------------------------------

chisquare = combined_df%>%
  annotate_driver_summary()%>%  
  filter_rs()%>%
  filter(study != "METABRIC")%>%
  filter_matched()%>%
  group_by(study, Tumor_Sample_Barcode, Consensus_annotated, clonality)%>%
  reframe(count = n())%>%
  filter(grepl("CLONAL", clonality))%>%
  tidyr::complete(Consensus_annotated, clonality, Tumor_Sample_Barcode)%>%mutate(study = if_else(substr(Tumor_Sample_Barcode, 1, 2) == "HF", "VAULT", "TCGA"))%>%
  pivot_wider(names_from = clonality, values_from = count)%>%
  pivot_longer(cols = c(CLONAL, SUBCLONAL), names_to = "clonality", values_to = "count")%>%
  mutate(count = if_else(is.na(count), 0, count))%>%#filter(study == "VAULT")%>%pull(Tumor_Sample_Barcode)%>%unique
  filter(Consensus_annotated != "Non_Driver")%>%
  group_by(study, Tumor_Sample_Barcode, clonality)%>%
  reframe(count = sum(count))%>%
  group_by(study, Tumor_Sample_Barcode)%>%
  mutate(sum = sum(count))%>%
  pivot_wider(names_from = clonality, values_from = count)%>%filter(!is.na(study))%>%
  group_by(study)%>%
  reframe(cases = n(),
          any_driver = sum(sum>0),
          clonal = sum(CLONAL>0),
          subclonal = sum(SUBCLONAL>0)
  )%>%
  pivot_longer(cols = c(any_driver, clonal, subclonal), names_to = "clonality", values_to = "count")%>%
  pivot_wider(names_from = study, values_from = c(count, cases))%>%
  rowwise()%>%
  mutate(p_value = (chisq.test(matrix(c(count_TCGA,cases_TCGA-count_TCGA, count_VAULT, cases_VAULT-count_VAULT), ncol = 2))$p.value),
         p_label =  ifelse(p_value < 0.001, "p < 0.001", paste0("p = ", round(p_value, 3))))

# DRAW PLOT ---------------------------------------------------------------
combined_df%>%
  #filter(!is.na(clonality)& clonality != "INDETERMINATE")%>%
  annotate_driver_summary()%>% #filter(study == "TCGA") %>%filter(Consensus_annotated == "Non_Driver")%>%pull(Hugo_Symbol)%>%table()%>%View()
  filter_rs()%>%
  filter(study != "METABRIC")%>%
  filter_matched()%>%
  group_by(study, Tumor_Sample_Barcode, Consensus_annotated, clonality)%>%
  reframe(count = n())%>%
  filter(grepl("CLONAL", clonality))%>%
  tidyr::complete(Consensus_annotated, clonality, Tumor_Sample_Barcode)%>%mutate(study = if_else(substr(Tumor_Sample_Barcode, 1, 2) == "HF", "VAULT", "TCGA"))%>%
  pivot_wider(names_from = clonality, values_from = count)%>%
  pivot_longer(cols = c(CLONAL, SUBCLONAL), names_to = "clonality", values_to = "count")%>%
  mutate(count = if_else(is.na(count), 0, count))%>%#filter(study == "VAULT")%>%pull(Tumor_Sample_Barcode)%>%unique
  filter(Consensus_annotated != "Non_Driver")%>%
  group_by(study, Tumor_Sample_Barcode, clonality)%>%
  reframe(count = sum(count))%>%
  group_by(study, Tumor_Sample_Barcode)%>%
  mutate(sum = sum(count))%>%
  pivot_wider(names_from = clonality, values_from = count)%>%filter(!is.na(study))%>%
  group_by(study)%>%
  reframe(cases = n(),
          any_driver = sum(sum>0),
          clonal = sum(CLONAL>0),
          subclonal = sum(SUBCLONAL>0)
          )%>%
  pivot_longer(cols = c(any_driver, clonal, subclonal), names_to = "clonality", values_to = "count")%>%
  mutate(prop = count/cases,
         p_label = chisquare$p_label[match(clonality, chisquare$clonality)],
         p_label = if_else(study == "VAULT", "", p_label))%>%
  ggplot(aes(study, prop, fill = clonality))+
  geom_col()+
  geom_text(aes(x = study, y = 1.1, label = p_label), nudge_x = 0.5)+
  facet_wrap(~clonality)+
  theme_classic(base_size = 25)+
  scale_fill_manual(values = c("subclonal" = subclonal_col, "clonal" = clonal_col))+
  rot.lab()+
  theme(plot.title = element_text(hjust = .5),
        legend.position = "none",
        axis.title.x = element_blank(),strip.text = element_text(size = 12),
        axis.text.x = element_text(size = 15))+
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  ylab("Patients (Proportion)")+
  ggtitle("Driver Detected")


# WRITE PLOT --------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Figure3G_proportion_patients_w_driver_by_study.png"))


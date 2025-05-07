# For VAULT, get the proportion of tumours with drivers

rm(list = ls(all = TRUE))
# LIBRARIES ---------------------------------------------------------------
library(tidyverse)

# PATHS -------------------------------------------------------------------

BASE = here::here()
BASE = "/Users/hanleyb/Documents/GitHub/VAULT_Breast/"
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
source(file.path(BASE, "src", "annotate_driver_summary.R"))

rs_col = "#54278F"
sr_col = "#33A02Cb7"
sr_col_1 = "#33A02Cb7"
sr_col_2 = "#B2DF8Aff"
gt_col = "grey"
clonal_col = "#33A02Cb7"
subclonal_col= "#54278F"
mix = colorRampPalette(c(rs_col, sr_col))
mix = mix(3)[2]  



# COMBINED DATAFRAME ------------------------------------------------------
vault$study = "VAULT"
common_names = Reduce(intersect,list(names(mbtcga), names(vault)))


combined_df = bind_rows(mbtcga%>%dplyr::select(common_names),
                        vault%>%dplyr::select(common_names)
)



# DRAW PLOT ---------------------------------------------------------------
combined_df%>%
  annotate_driver_summary()%>%
  filter(study == "VAULT")%>%
  filter_rs()%>%
  group_by(Tumor_Sample_Barcode, clonality, Consensus_annotated)%>%
  reframe(count = n())%>%
  pivot_wider(names_from = Consensus_annotated, values_from = count)%>%
  dplyr::select(-Non_Driver)%>%
  rowwise()%>%
  mutate(driver_sum_by_clonality = sum(Driver_Other, Biomarker, na.rm = T))%>%
  dplyr::select(-c(Driver_Other, Biomarker))%>%
  ungroup()%>%
  group_by(Tumor_Sample_Barcode)%>%
  pivot_wider(names_from = clonality, values_from = driver_sum_by_clonality)%>%
  dplyr::select(CLONAL, SUBCLONAL, INDETERMINATE, Tumor_Sample_Barcode)%>%
  mutate(driver = case_when(CLONAL>0 & SUBCLONAL ==0 ~ "Clonal_only",
                            CLONAL>0 & SUBCLONAL >0 ~ "Clonal&subclonal",
                            CLONAL==0 & SUBCLONAL >0 ~ "Subclonal_only",
                            CLONAL==0 & SUBCLONAL ==0 ~ "No driver identified"
  ))%>%
  group_by(driver)%>%
  reframe(count = n())%>%
  mutate(driver = factor(driver, levels = c("Clonal_only","Subclonal_only", "Clonal&subclonal","No driver identified")))%>%
  ggplot(aes("x", count, fill = driver))+
  geom_col()+
  theme_classic(base_size = 20)+
  scale_fill_manual(values = c(
    "Clonal_only" = sr_col,
    "Clonal&subclonal" = mix,
    "Subclonal_only" = rs_col,
    "No driver identified" = "lightgrey"
  ))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "bottom")+
  ylab("VAULT patients - Count")+
  guides(fill =  guide_legend(nrow = 4))


# WRITE PLOT --------------------------------------------------------------

ggsave(file.path(OUT_DIR, "Figure3A_prop_cases_w_drivers.png"))


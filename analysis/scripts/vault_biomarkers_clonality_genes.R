# For VAULT biomarkers - identify which genes they are in

rm(list = ls(all = TRUE))
# LIBRARIES ---------------------------------------------------------------
library(tidyverse)

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
source(file.path(BASE, "src", "annotate_driver_summary.R"))
source(file.path(BASE, "src", "plotting_features.R"))


# GRAPHICAL PARAMETERS ----------------------------------------------------
rs_col = "#54278F"
sr_col = "#33A02Cb7"
sr_col_1 = "#33A02Cb7"
sr_col_2 = "#B2DF8Aff"
gt_col = "grey"
clonal_col = "#33A02Cb7"
subclonal_col= "#54278F"
mix = colorRampPalette(c(rs_col, sr_col))
mix = mix(3)[2]  



# COMBINED DATASET --------------------------------------------------------
vault$study = "VAULT"
common_names = Reduce(intersect,list(names(mbtcga), names(vault)))
combined_df = bind_rows(mbtcga%>%dplyr::select(all_of(common_names)),
                        vault%>%dplyr::select(all_of(common_names))
)



# DRAW PLOT ---------------------------------------------------------------
combined_df%>%
  annotate_driver_summary()%>%  
  filter(study == "VAULT")%>%
  filter_rs()%>%
  filter(actionable_summary != "Non_Driver")%>%
  group_by(clonality, actionable_summary)%>%
  filter(Consensus_annotated == "Biomarker")%>%
  group_by(Hugo_Symbol, clonality)%>%
  reframe(count = n())%>%
  filter(grepl("CLONAL", clonality))%>%
  ggplot(aes("x", count, fill = Hugo_Symbol))+
  geom_col()+
  facet_wrap(~clonality, scales = "free_y")+
  theme_classic(base_size = 20)+
  ggtitle("Biomarker Variants")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), legend.text = element_text(size = 10),
        legend.title = element_blank())+
  scale_fill_manual(values = colours)+
  ylab("Variant Count")


# WRITE PLOT --------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Figure3C_genes_drivers_clonality.png"))


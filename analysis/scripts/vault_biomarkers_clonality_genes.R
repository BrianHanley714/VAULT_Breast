# For VAULT biomarkers - identify which genes they are in

rm(list = ls(all = TRUE))
# LIBRARIES ---------------------------------------------------------------
library(tidyverse)

# PATHS -------------------------------------------------------------------
BASE = "/Users/hanleyb/Dropbox (The Francis Crick)/HoLSTF_Breast/Github_Repo"
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
matched_patients = read.delim(MATCHED_PT_MBTCGA)
matched_patients_char = read.delim(MATCHED_PT_CHAR_MBTCGA)
tumour_IDs = read.delim(IDMAP)
matched_patients = c(matched_patients$PATIENT_ID, unique(tumour_IDs$sample[match(matched_patients$PATIENT_ID, tumour_IDs$metabricId)]))
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


vault$study = "VAULT"
common_names = Reduce(intersect,list(names(mbtcga), names(vault)))


combined_df = bind_rows(mbtcga%>%dplyr::select(common_names),
                        vault%>%dplyr::select(common_names)
)

colours = c("#B2DF8A93", "#FDBF6Fdb", "#FB9A99db", "#1F78B4db", "#B15928ff",
            "#6A3D9Aff", "#FB9A9993", "#FF7F00db", "#CAB2D6ff", "#33A02Cdb", "#FB9A99ff",
            "#E31A1Cff", "#A6CEE3db", "#FF7F00ff", "#B2DF8Ab7", "#B15928b7", "#CAB2D6db",
            "#33A02Cb7", "#E31A1Cdb", "#A6CEE3ff", "#33A02Cff", "#CAB2D6b7", "#B2DF8Aff",
            "#1F78B4ff", "#6A3D9Adb", "#FFFF99db", "#B15928db", "#E31A1Cb7", "#FF7F00b7",
            "#FFFF99b7", "#A6CEE393", "#B2DF8A6f", "#FFFF99ff", "#B2DF8Adb")


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


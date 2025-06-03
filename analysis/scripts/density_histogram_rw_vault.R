# Density histogram of mutational cancer cell fraction

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


rs_col = "#54278F"
sr_col = "#33A02Cb7"
sr_col_1 = "#33A02Cb7"
sr_col_2 = "#B2DF8Aff"
gt_col = "grey"
clonal_col = "#33A02Cb7"
subclonal_col= "#54278F"
vault$study = "VAULT"


# COMBINED DATAFRAME ------------------------------------------------------
common_names = Reduce(intersect,list(names(mbtcga), names(vault)))

combined_df = bind_rows(mbtcga%>%dplyr::select(common_names),
                        vault%>%dplyr::select(common_names)
)

combined_df%>%
  filter_matched()%>%
  filter_rs()%>%
  mutate(study = factor(study, levels = c("VAULT", "TCGA", "METABRIC")))%>%
  ggplot(aes(ccf_expected_copies, fill = study))+
  # geom_density()+
  stat_density(geom = "area", position = "identity", alpha = 0.7, aes(y = ..scaled..))+
  scale_fill_manual(values= c("VAULT" = rs_col, "METABRIC" = sr_col_1, TCGA = sr_col_2))+
  theme_classic(base_size = 20)+
  xlab("Cancer Cell Fraction")+
  theme(legend.position = "top",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))+
  #  ggtitle(paste("Sampling for", run))+
  ylab("Mutational density (scaled)")



# WRITE PLOT --------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Figure2D_CCF_mutational_density.png"))



# Create a figure showing the selected genes by CCF

rm(list = ls(all = TRUE))
# LIBRARIES ---------------------------------------------------------------
library(tidyverse)

# PATHS -------------------------------------------------------------------

BASE = here::here()
BASE = "/Users/hanleyb/Dropbox (The Francis Crick)/HoLSTF_Breast/Github_Repo"
OUT_DIR = file.path(BASE, "analysis", "figures")
VARIANTS_VAULT = file.path(BASE, "data","variants", "variant_calls_VAULT.txt")
VARIANTS_MBTCGA = file.path(BASE, "data","variants", "variant_calls_TCGA_MB.txt")
MATCHED_PT_MBTCGA = file.path(BASE, "data","metadata", "matched_patients.txt")
MATCHED_PT_CHAR_MBTCGA = file.path(BASE, "data","metadata", "matched_patients_characteristics.txt")
IDMAP = file.path(BASE, "data","metadata", "tumouridmap_MB.txt")
INCLUDED_PATIENTS = file.path(BASE, "data","metadata", "cases_included.xlsx")
CLINDATA = file.path(BASE, "data", "metadata", "clinical_data.txt")
MUTSIGPATH = file.path(BASE, "data", "variants", "MutSig2CV_sig_genes.txt")
DNDSPATH = file.path(BASE, "data", "variants", "dndscv_selected_genes.txt")

# LOAD DATA ---------------------------------------------------------------
rs_patients = read.delim(INCLUDED_PATIENTS)[,1]
vault = read.delim(VARIANTS_VAULT)
mbtcga = read.delim(VARIANTS_MBTCGA)
matched_patients = read.delim(MATCHED_PT_MBTCGA)
matched_patients_char = read.delim(MATCHED_PT_CHAR_MBTCGA)
tumour_IDs = read.delim(IDMAP)
matched_patients = c(matched_patients$PATIENT_ID, unique(tumour_IDs$sample[match(matched_patients$PATIENT_ID, tumour_IDs$metabricId)]))
clinical_data = read.delim(CLINDATA)
selected = read.delim(DNDSPATH)[,1]
mutsig2cv = read.delim(MUTSIGPATH)

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


# SIGNIFICANTLY MUTATED GENES------------------------------------------------

mutsig2cv = mutsig2cv$gene[mutsig2cv$q<0.25]
drivers = unique(c(mutsig2cv, selected))

# get order for plotting
order = 
  combined_df%>%
  filter(!is.na(ccf_expected_copies))%>%
  filter_rs()%>%
  filter(study == "VAULT")%>%
  filter(Hugo_Symbol %in% drivers)%>%
  mutate(id = paste(Hugo_Symbol, Start_Position, Tumor_Sample_Barcode))%>%
  group_by(Hugo_Symbol)%>%
  mutate(mean = mean(ccf_expected_copies, na.rm = T))%>%
  arrange(mean)%>%pull(Hugo_Symbol)%>%unique()


# DRAW PLOT ---------------------------------------------------------------
combined_df%>%
  filter(!is.na(ccf_expected_copies))%>%
  filter_rs()%>%
  filter(study == "VAULT")%>%
  filter(Hugo_Symbol %in% drivers)%>%
  mutate(id = paste(Hugo_Symbol, Start_Position, Tumor_Sample_Barcode))%>%
  group_by(Hugo_Symbol)%>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, levels = rev(order)))%>%
  mutate(text_colour = if_else(Hugo_Symbol %in% c("NELFE", "PTH2", "MYEOV"), "new", "known"))%>%
  ggplot(aes(reorder(id, -ccf_expected_copies_lower), ccf_expected_copies, col = clonality))+
  geom_point()+
  geom_errorbar(aes(ymax = ccf_expected_copies_upper,
                    ymin = ccf_expected_copies_lower), 
                width = 0.2)+
  facet_wrap(~Hugo_Symbol+text_colour, labeller = 
               labeller(
                 Hugo_Symbol = ~ paste(.),
                 text_colour = ~ paste(.),
                 .multi_line = T), scales = "free_x", nrow=1)+
  theme_classic(base_size = 15)+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        strip.text = element_text(face = "bold")
  )+
  ylab("Cancer Cell Fraction")+
  scale_color_manual(values = c("CLONAL" = clonal_col, "SUBCLONAL" = subclonal_col))




# WRITE PLOT --------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Figure3J_clonality_selected_drivers.png"))

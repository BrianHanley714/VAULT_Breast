# Create a plot looking at the proportion of subclonal variants per driver type

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
ENDPOINT = file.path(BASE, "data", "metadata", "VAULT_endpoint_data.txt")
VAULT_UNFILTERED = file.path(BASE, "data","variants", "exposure_signatures_unfiltered_data.txt")
VAULT_UNFILTERED_CSI = file.path(BASE, "data","variants", "cosinse_sim_signatures_unfiltered_data.txt")

# LOAD DATA ---------------------------------------------------------------
rs_patients = read.delim(INCLUDED_PATIENTS)[,1]
vault = read.delim(VARIANTS_VAULT)
mbtcga = read.delim(VARIANTS_MBTCGA)
matched_patients_char = read.delim(MATCHED_PT_CHAR_MBTCGA)
tumour_IDs = read.delim(IDMAP)
matched_patients = c(matched_patients_char$PATIENT_ID, unique(tumour_IDs$sample[match(matched_patients_char$PATIENT_ID, tumour_IDs$metabricId)]))
matched_patients = matched_patients[!is.na(matched_patients)]
clinical_data = read.delim(CLINDATA)
ffpe_exposure = read.delim(VAULT_UNFILTERED)
ffpe_csi = read.delim(VAULT_UNFILTERED_CSI)
endpoint_data = read.delim(ENDPOINT)
# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "custom_filters.R"))
source(file.path(BASE, "src", "rot.lab.R"))
source(file.path(BASE, "src", "annotate_driver_summary.R"))

rs_col = "#54278F"
sr_col = "#33A02Cb7"



# COMBINED DATAFRAME ------------------------------------------------------
# vault$study = "VAULT"
# common_names = Reduce(intersect,list(names(mbtcga), names(vault)))
# 
# 
# combined_df = bind_rows(mbtcga%>%dplyr::select(common_names),
#                         vault%>%dplyr::select(common_names)
# )


base_size = 20




# DRAW PLOTS --------------------------------------------------------------
plot_grid(
  ffpe_exposure%>%
    mutate(tumour_sample = rownames(.))%>%
    mutate(tumour_sample = toupper(sub(".unfiltered","", sub("^[^_]*_", "", tumour_sample))))%>%
    filter(tumour_sample %in% rs_patients)%>%
    mutate(Trial.ID = substr(tumour_sample, 1, 5))%>%
    left_join(endpoint_data, by = "Trial.ID")%>%
    ggplot(aes(Time..in.formalin..days., ffpe_sig_unrepaired))+
    geom_point(size = 4)+
    stat_cor(size = 8)+
    theme_classic(base_size = base_size)+
    geom_smooth(method = "lm", 
                col = "black")+
    ylab("FFPEsig_exposure")+
    rot.lab()+
    theme(axis.text.x = element_text(size = 15))+
    xlab("Time in formalin (days)"),
ffpe_csi%>%
  mutate(type = paste(study, filtered_variants, sep = "_"),
         type = factor(type, levels = c("VAULT_unfiltered", "VAULT_filtered", "TCGA_filtered")))%>%
  ggplot(aes(type, simil))+
  geom_violin(aes(fill = type), width= 1.1)+
  geom_boxplot(width = 0.5, outlier.alpha = 0)+
  geom_jitter(aes(type, simil), width = 0.11)+
  theme_classic(base_size = base_size)+
  stat_compare_means(label = "p.signif", comparisons = list(c("VAULT_unfiltered", "VAULT_filtered"), c("TCGA_filtered", "VAULT_unfiltered")), size = 10)+
  scale_fill_manual(values = c("VAULT_unfiltered" = rs_col, "VAULT_filtered" = rs_col, "TCGA_filtered" = sr_col))+
  ylab("Cosine Similarity to FFPESig")+
  rot.lab()+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15)),

labels = LETTERS[1:2],
label_size = 25
)

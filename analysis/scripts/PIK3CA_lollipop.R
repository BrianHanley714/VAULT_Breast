# Create lollipop plot of PICK3CA variants in VAULT, MB and TCGA

rm(list = ls(all = TRUE))
# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(gggenes)
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
PIK3CA_DOM = file.path(BASE, "data", "metadata", "PIK3CAdomains_P42336_EBI_10022025.tsv")

# LOAD DATA ---------------------------------------------------------------
rs_patients = read.delim(INCLUDED_PATIENTS)[,1]
vault = read.delim(VARIANTS_VAULT)
mbtcga = read.delim(VARIANTS_MBTCGA)
matched_patients = read.delim(MATCHED_PT_MBTCGA)
matched_patients_char = read.delim(MATCHED_PT_CHAR_MBTCGA)
tumour_IDs = read.delim(IDMAP)
matched_patients = c(matched_patients$PATIENT_ID, unique(tumour_IDs$sample[match(matched_patients$PATIENT_ID, tumour_IDs$metabricId)]))
clinical_data = read.delim(CLINDATA)
pik3ca_domains = read.delim(PIK3CA_DOM)

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



# COMBINED DATAFRAME ------------------------------------------------------
vault$study = "VAULT"
common_names = Reduce(intersect,list(names(mbtcga), names(vault)))


combined_df = bind_rows(mbtcga%>%dplyr::select(common_names),
                        vault%>%dplyr::select(common_names)
)
# PIK3CA

pik3ca_domains = pik3ca_domains%>%
  filter(Source.Database == "pfam")%>%
  mutate(start = as.numeric(sub("\\.\\..*", "", Matches)),
         end = as.numeric(sub(".*\\.\\.", "", Matches)),
         gene = "PIK3CA")%>%
  mutate(Name = case_when(Name == "Phosphatidylinositol 3- and 4-kinase" ~"Kinase",
                          Name == "Phosphoinositide 3-kinase family, accessory domain (PIK domain)"~"Accessory",
                          Name == "Phosphoinositide 3-kinase C2"~ "C2",
                          Name == "PI3-kinase family, ras-binding domain"~ "ras-binding",
                          Name == "PI3-kinase family, p85-binding domain"~ "p85-binding")
  )
mutations = combined_df%>%
  filter(!is.na(clonality)& clonality != "INDETERMINATE")%>%
  # mutate(ONCOGENIC = if_else(grepl("Onco", ONCOGENIC), ONCOGENIC, "Non_Driver"))%>%
  # mutate(Sensitivity_biomarker = if_else(HIGHEST_LEVEL!= "" &!is.na(HIGHEST_LEVEL), "Therapeutic_Biomarkers", ONCOGENIC),
  #        Sensitivity_biomarker = if_else(grepl("Non_Driver|Therapeutic_Biomarkers", Sensitivity_biomarker), Sensitivity_biomarker, "Driver_Other"),
  # )%>%
  filter_rs()%>%
  filter_matched()%>%
  mutate(study = if_else(study %in% c("METABRIC", "TCGA"), "SingReg", "RepSamp"))%>%
  filter(Hugo_Symbol == "PIK3CA")%>%
  dplyr::select(AMINO_ACID_START, study, clonality, HGVSp_short) 

# PIK3CA lollipop ---------------------------------------------------------


mutations = mutations%>%
  group_by(study, clonality, AMINO_ACID_START)%>%
  reframe(Frequency = n())%>%
  ungroup()%>%
  mutate(AMINO_ACID_START = as.numeric(AMINO_ACID_START))%>%
  mutate(Frequency = if_else(clonality == "SUBCLONAL", Frequency - (2*Frequency), Frequency),
         study = factor(study, levels = c("SingReg", "RepSamp")),
         clonality =if_else(study == "SingReg", "TCGA/METABRIC", clonality),
         clonality = factor(clonality, levels = c("TCGA/METABRIC", 
                                                  "CLONAL", 
                                                  "SUBCLONAL")),
         Frequency2 = if_else(clonality == "TCGA/METABRIC", NA, Frequency)
  )




# DRAW PLOT ---------------------------------------------------------------
ggplot() +
  geom_segment(data = mutations, aes(x = AMINO_ACID_START, xend = AMINO_ACID_START, y = 0, yend = Frequency, color = clonality), size = 1) +
  geom_point(data = mutations, aes(x = AMINO_ACID_START, y = Frequency, color = clonality, fill = clonality), 
             size = 1, shape = 21,  stroke = 1.2) +
  geom_segment(data = mutations, aes(x = AMINO_ACID_START, xend = AMINO_ACID_START, y = 0, yend = Frequency2, color = clonality), size = 1) +
  geom_point(data = mutations, aes(x = AMINO_ACID_START, y = Frequency2, color = clonality, fill = clonality), 
             size = 1, shape = 21,  stroke = 1.2) +
  geom_rect(aes(xmin = 0, xmax = unique(pik3ca_domains$Protein.Length), ymin = -0.1, ymax = 0.1), fill = "grey", col = "black")+
  
  geom_gene_arrow(data = pik3ca_domains, 
                  aes(xmin = start,xmax = end, y = 0, fill = Name),  arrowhead_width = grid::unit(x = 0, "mm"), arrowhead_height = grid::unit(0, "mm")) +
  scale_fill_manual(values = c("Kinase" = colours[1],
                               "Accessory" = colours[2],
                               "C2" = colours[3],
                               "ras-binding" = colours[4],
                               "p85-binding" = colours[5],
                               "CLONAL" = clonal_col, "SUBCLONAL" = subclonal_col,
                               "TCGA/METABRIC" = "lightgrey",
                               "RepSamp" = rs_col),
                    breaks = unique(pik3ca_domains$Name))+
  scale_y_continuous(trans = "pseudo_log", breaks = c(-1, -10, -100, 0, 1, 10, 100)) +
  scale_color_manual(values = c("TCGA/METABRIC" = "lightgrey", "CLONAL" = clonal_col, "SUBCLONAL" = subclonal_col,
                                "SingReg" = "lightgrey",
                                "RepSamp" = rs_col), )+
  theme_classic(base_size = 20)+
  xlab("Amino Acid Position")+
  ylab("Variant Count")+
  ggtitle("PIK3CA")+
  guides(col = "none")+
  labs(fill = NULL)+
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold"))+
  annotate("text", x = unique(pik3ca_domains$Protein.Length)*0.05, y = max(mutations$Frequency), label = "CLONAL", col = clonal_col, fontface = "bold", size = 6)+
  annotate("text", x = unique(pik3ca_domains$Protein.Length)*0.05, y = min(mutations$Frequency), label = "SUBCLONAL", col = subclonal_col, fontface = "bold", size = 6)


# WRITE PLOT --------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Figure3D_PIK3CA_lollipop.png"))


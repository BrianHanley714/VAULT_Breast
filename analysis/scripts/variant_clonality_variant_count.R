# Create a geom_col() plot summarizing variant clonality

rm(list = ls(all = TRUE))
# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(cowplot)
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

# LOAD DATA ---------------------------------------------------------------
rs_patients = read.delim(INCLUDED_PATIENTS)[,1]
vault = read.delim(VARIANTS_VAULT)
mbtcga = read.delim(VARIANTS_MBTCGA)
matched_patients = read.delim(MATCHED_PT_MBTCGA)
matched_patients_char = read.delim(MATCHED_PT_CHAR_MBTCGA)
tumour_IDs = read.delim(IDMAP)
matched_patients = c(matched_patients$PATIENT_ID, unique(tumour_IDs$sample[match(matched_patients$PATIENT_ID, tumour_IDs$metabricId)]))

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
common_names = Reduce(intersect,list(names(mbtcga), names(vault)))


# DRAW PLOT ---------------------------------------------------------------
combined_df = bind_rows(mbtcga%>%dplyr::select(common_names),
                        vault%>%dplyr::select(common_names)
)


# DRAW PLOT ---------------------------------------------------------------
plot_grid(combined_df%>%
            filter_matched()%>%
            filter_rs()%>%
            filter(!is.na(clonality))%>%
            group_by(study, Tumor_Sample_Barcode, clonality)%>%
            reframe(num_mut = n())%>%
            group_by(Tumor_Sample_Barcode)%>%
            mutate(sum = sum(num_mut),
                   prop = num_mut/sum,
                   order = if_else(clonality == "SUBCLONAL", prop, 0),
                   order = max(order))%>%
            filter(!duplicated(Tumor_Sample_Barcode))%>%
            ggplot(aes(reorder(Tumor_Sample_Barcode, order), sum))+
            geom_col(width=1)+
            scale_y_log10()+
            ylab("Count")+
            theme_classic(base_size = 19)+
            theme(axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank()),
          
          combined_df%>%
            filter_matched()%>%
            filter_rs()%>%
            filter(!is.na(clonality))%>%
            group_by(study, Tumor_Sample_Barcode, clonality)%>%
            reframe(num_mut = n())%>%
            group_by(Tumor_Sample_Barcode)%>%
            mutate(sum = sum(num_mut),
                   prop = num_mut/sum,
                   order = if_else(clonality == "SUBCLONAL", prop, 0),
                   order = max(order))%>%
            ggplot()+
            geom_col(aes(reorder(Tumor_Sample_Barcode, order), prop, fill = clonality), width = 1)+
            geom_point(aes(reorder(Tumor_Sample_Barcode, order), -0.05, col = study), shape = "|", size = 5)+
            theme_classic(base_size = 20)+
            theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank())+
            scale_color_manual(values= c("VAULT" = rs_col, "METABRIC" = sr_col_1, TCGA = sr_col_2))+
            scale_fill_brewer(palette = "Set3")+
            ylab("Variant Proportion")+
            xlab("Patient")+
            theme(legend.position = "bottom")+
            guides(color = guide_legend(override.aes = list(shape = 20))),
          ncol = 1,
          rel_heights = c(1,6)
)


# WRITE PLOT --------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Figure2E_variant_clonality.png"))

# create lollipop plots for the PTH2 and MYEOV variants comparing to positions in Cosmic

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
CLINDATA = file.path(BASE, "data", "metadata", "clinical_data.txt")
PTH_COSM = file.path(BASE, "data", "variants", "PTH2_variants_Cosmic_10022025.tsv")
PTH_DOM = file.path(BASE, "data", "metadata", "PTH2domains_Q96A98_EBI_10022025.tsv")
MYEOV_COSM = file.path(BASE, "data", "variants", "MYEOV_variants_Cosmic_10022025.tsv")

# LOAD DATA ---------------------------------------------------------------
rs_patients = read.delim(INCLUDED_PATIENTS)[,1]
vault = read.delim(VARIANTS_VAULT)
mbtcga = read.delim(VARIANTS_MBTCGA)
matched_patients = read.delim(MATCHED_PT_MBTCGA)
matched_patients_char = read.delim(MATCHED_PT_CHAR_MBTCGA)
tumour_IDs = read.delim(IDMAP)
matched_patients = c(matched_patients$PATIENT_ID, unique(tumour_IDs$sample[match(matched_patients$PATIENT_ID, tumour_IDs$metabricId)]))
clinical_data = read.delim(CLINDATA)
PTH2_variants_Cosmic = read.delim(PTH_COSM)
PTH2_domains = read.delim(PTH_DOM)
MYEOV_variants_Cosmic = read.delim(MYEOV_COSM)

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "custom_filters.R"))
source(file.path(BASE, "src", "rot.lab.R"))
source(file.path(BASE, "src", "plotting_features.R"))
#source(file.path(BASE, "src", "annotate_driver_summary.R"))

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


# PTH2 Lollipop -----------------------------------------------------------
PTH2_domains = PTH2_domains%>%
  filter(Source.Database == "pfam")%>%
  mutate(start = as.numeric(sub("\\.\\..*", "", Matches)),
         end = as.numeric(sub(".*\\.\\.", "", Matches)),
         gene = "PTH2")%>%
  mutate(Name = case_when(Name == "TIP39 peptide" ~"TIP39")
  )


mutations_PTH2 = bind_rows(PTH2_variants_Cosmic%>%
                             group_by(Position)%>%
                             reframe(Count = sum(Count))%>%
                             mutate(study = "COSMIC", 
                                    Position = as.integer(Position)),
                           
                           combined_df%>%
                             filter_rs()%>%
                             filter(study == "VAULT")%>%
                             filter(Hugo_Symbol == "PTH2")%>%
                             dplyr::select(Position = AMINO_ACID_START, study)%>%
                             group_by(Position)%>%
                             reframe(Count = n())%>%
                             mutate(study = "VAULT", 
                                    Position = as.integer(Position))
)%>%
  mutate(Count2 = as.numeric(if_else(study == "VAULT", Count, NA)))


# MYEOV Lollipop -----------------------------------------------------------
# not able to find pfam domains for MYEOV
PTH2_domains = PTH2_domains%>%
  filter(Source.Database == "pfam")%>%
  mutate(start = as.numeric(sub("\\.\\..*", "", Matches)),
         end = as.numeric(sub(".*\\.\\.", "", Matches)),
         gene = "PTH2")%>%
  mutate(Name = case_when(Name == "TIP39 peptide" ~"TIP39")
  )

mutations
mutations_MYEOV = bind_rows(MYEOV_variants_Cosmic%>%
                              group_by(Position)%>%
                              reframe(Count = sum(Count))%>%
                              mutate(study = "COSMIC", 
                                     Position = as.integer(Position)),
                            
                            combined_df%>%
                              filter_rs()%>%
                              filter(study == "VAULT")%>%
                              filter(Hugo_Symbol == "MYEOV")%>%
                              dplyr::select(Position = AMINO_ACID_START, study)%>%
                              group_by(Position)%>%
                              reframe(Count = n())%>%
                              mutate(study = "VAULT", 
                                     Position = as.integer(Position))
)%>%
  mutate(Count2 = as.numeric(if_else(study == "VAULT", Count, NA)))



# DRAW PLOT ---------------------------------------------------------------
plot_grid(
  ggplot() +
    geom_segment(data = mutations_PTH2, aes(x = Position, xend = Position, y = 0, yend = Count, color = study), size = 1) +
    geom_point(data = mutations_PTH2, aes(x = Position, y = Count, color = study, fill = study), 
               size = 1, shape = 21,  stroke = 1.2) +
    geom_segment(data = mutations_PTH2, aes(x = Position, xend = Position, y = 0, yend = Count2, color = study), size = 1)+ 
    geom_point(data = mutations_PTH2, aes(x = Position, y = Count2, color = study, fill = study), 
               size = 1, shape = 21,  stroke = 1.2) +
    geom_rect(aes(xmin = 0, xmax = 101, ymin = -5, ymax = 0), fill = "grey", col = "black")+
    
    # geom_gene_arrow(data = PTH2_domains, 
    #                 aes(xmin = start,xmax = end, y = 0.45, fill = Name),  arrow_body_height = grid::unit(x = 14, "mm"), arrowhead_width = grid::unit(x = 0, "mm"), arrowhead_height = grid::unit(0, "mm")) +
    scale_fill_manual(values = c("TIP39" = colours[1], "COSMIC" = "lightgrey", VAULT = rs_col),
                      breaks = unique(PTH2_domains$Name))+
    
    scale_color_manual(values = c("COSMIC" = "lightgrey", VAULT = rs_col))+
    theme_classic(base_size = 20)+
    #scale_y_log10()+
    xlab("Amino Acid Position")+
    ylab("Variant Count")+
    ggtitle("PTH2")+
    labs(fill = NULL)+
    theme(legend.position = c(0.9, 0.8),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.1, vjust = -10, face = "bold")),
  
  ggplot() +
    geom_segment(data = mutations_MYEOV, aes(x = Position, xend = Position, y = 0, yend = Count, color = study), size = 1) +
    geom_point(data = mutations_MYEOV, aes(x = Position, y = Count, color = study, fill = study), 
               size = 1, shape = 21,  stroke = 1.2) +
    geom_segment(data = mutations_MYEOV, aes(x = Position, xend = Position, y = 0, yend = Count2, color = study), size = 1)+ 
    geom_point(data = mutations_MYEOV, aes(x = Position, y = Count2, color = study, fill = study), 
               size = 1, shape = 21,  stroke = 1.2) +
    geom_rect(aes(xmin = 0, xmax = 311, ymin = -5, ymax = 0), fill = "grey", col = "black")+
    
    # geom_gene_arrow(data = MYEOV_domains, 
    #                 aes(xmin = start,xmax = end, y = 0.45, fill = Name),  arrow_body_height = grid::unit(x = 14, "mm"), arrowhead_width = grid::unit(x = 0, "mm"), arrowhead_height = grid::unit(0, "mm")) +
    scale_fill_manual(values = c( "COSMIC" = "lightgrey", VAULT = rs_col))+
    # 
    scale_color_manual(values = c("COSMIC" = "lightgrey", VAULT = rs_col))+
    theme_classic(base_size = 20)+
    # scale_y_log10()+
    xlab("Amino Acid Position")+
    ylab("Variant Count")+
    ggtitle("MYEOV")+
    
    labs(fill = NULL)+
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.1, vjust = -10, face = "bold"),
          plot.margin = margin(t = 0, r = 10, b = 10, l = 10))+
    guides(fill = "none"),
  ncol = 1,
  axis = "v"
)
  
  

# WRITE PLOT --------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Figure3I_MYEOV_PTH2_lollipop.png"))


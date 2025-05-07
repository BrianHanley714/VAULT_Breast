
rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(ggpubr)

# PATHS -------------------------------------------------------------------
BASE = here::here()
BASE = "/Users/hanleyb/Documents/GitHub/VAULT_Breast/"
OUT_DIR = file.path(BASE, "analysis", "figures")
VARIANTS_VAULT = file.path(BASE, "data","variants", "variant_calls_VAULT.txt")
INCLUDED_PATIENTS = file.path(BASE, "data","metadata", "cases_included.xlsx")
CLINDATA = file.path(BASE, "data", "metadata", "clinical_data.txt")


# LOAD DATA ---------------------------------------------------------------
rs_patients = read.delim(INCLUDED_PATIENTS)[,1]
vault = read.delim(VARIANTS_VAULT)
clinical_data = read.delim(CLINDATA)


# ORGANISE DATA -----------------------------------------------------------
TMB_cases = vault%>%
  group_by(Tumor_Sample_Barcode)%>%
  count()%>%
  mutate(TMB = n/36.7)%>%
  filter(TMB>10)%>%
  pull(Tumor_Sample_Barcode)

levels = vault%>%
  filter(Tumor_Sample_Barcode %in% TMB_cases)%>%
  filter(clonality == "CLONAL")%>%
  group_by(Tumor_Sample_Barcode)%>%
  count%>%
  mutate(Tumor_Sample_Barcode = substr(Tumor_Sample_Barcode, 1, 5))%>%
  arrange(-n)%>%
  
  pull(Tumor_Sample_Barcode)


df_plot = vault%>%
  filter(Tumor_Sample_Barcode %in% TMB_cases)
  

base_size = 25 # set base_size

# create bar plot
plot= df_plot%>%
  mutate(Tumor_Sample_Barcode = substr(Tumor_Sample_Barcode, 1, 5))%>%
  mutate(Tumor_Sample_Barcode = factor(Tumor_Sample_Barcode,levels = levels))%>%
  filter(!is.na(clonality))%>%
  group_by(clonality, Tumor_Sample_Barcode)%>%
  count%>%
  group_by(Tumor_Sample_Barcode)%>%
  mutate(sum = sum(n),
         prop = n/sum)%>%
  ggplot(aes(prop, Tumor_Sample_Barcode, fill = clonality))+
  geom_col()+
  scale_fill_brewer(palette = "Set2")+
  theme_classic(base_size = base_size)


# extract legend from plot
legend <- get_legend(
plot+
  theme(legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())+
  xlab("Proportion")
)




# DRAW PLOT ---------------------------------------------------------------
plot_grid(
  plot_grid(
  df_plot%>%
    mutate(Tumor_Sample_Barcode = substr(Tumor_Sample_Barcode, 1, 5))%>%
    mutate(Tumor_Sample_Barcode = factor(Tumor_Sample_Barcode,levels = levels))%>%
    ggplot(aes(ccf_expected_copies))+
    geom_density()+
    facet_wrap(~Tumor_Sample_Barcode, ncol = 1, scales = "free_y", strip.position = "left", as.table = F)+
    theme_classic(base_size = base_size) +
    theme(legend.position = "none",
                                                axis.text.y = element_blank(),
                                                axis.title.y = element_blank(),
                                                axis.ticks.y = element_blank())+
    xlab("Cancer Cell Fraction")
  ,
  
plot+
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())+
  xlab("Proportion")),
legend,
ncol = 1,
rel_heights = c(1,.05)
)



# WRITE PLOT --------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Extended_Data_Figure5_high_TMB_cases.png"))

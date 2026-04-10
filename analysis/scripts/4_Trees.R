rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(CONIPHER)
library(igraph)
library(qs)
library(openxlsx)

# LOAD DATA ---------------------------------------------------------------
BASE = here::here()
BASE = "/Users/hanleyb/Documents/GitHub/VAULT_Breast"
OUT_DIR = file.path(BASE, "analysis", "figures")
TREES_PATHS = file.path(BASE, "data", "Trees")
VAULT_CLIN_PATH = file.path(BASE, "data", "metadata", "clinical_data.txt")
VAULT_VARIANTS_PATH = file.path(BASE, "data", "variants", "variant_calls_VAULT.rdata")
CASES_INCL_PATH = file.path(BASE, "data", "metadata", "cases_included.xlsx")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_features.R"))
source(file.path(BASE, "src", "custom_tree_plotting.R"))
source(file.path(BASE, "src", "conipher_helper_functions.R"))


# LOAD DATA ---------------------------------------------------------------
trees_list = list.files(TREES_PATHS, pattern = "_ck_tree.RDS", full.names = T)
trial_ids = list.files(TREES_PATHS, pattern = "_ck_tree.RDS")%>%substr(1, 5)
trees_list = lapply(trees_list, function(hf){
  tree = readRDS(hf)
  return(tree)
})
names(trees_list) = trial_ids
clinical_data = read.delim(VAULT_CLIN_PATH)
maf_out = qread(VAULT_VARIANTS_PATH)
ck_included_cases_vault = read.xlsx(CASES_INCL_PATH, colNames = F)[,1]

drivers = maf_out%>%dplyr::filter(Consensus_annotated != "Non_Driver")

drivers = drivers%>%
  mutate(variant_name = paste(Trial_ID, Hugo_Symbol, Chromosome, Start_Position, ALT.x, sep = ";"))

drivers$rescued = "No"

drivers_called_rs = drivers%>%
  dplyr::filter(Tumor_Sample_Barcode %in% ck_included_cases_vault,
                rescued == "No"
  )%>%pull(variant_name)

drivers_called_ffpe = drivers%>%
  dplyr::filter(SingReg_Prim_v_Met %in% c("PRIMARY", "RECURRENCE_RS"),
                rescued == "No"
  )%>%pull(variant_name)

drivers = drivers%>%
  mutate(called_in_rs = variant_name %in% drivers_called_rs,
         called_in_ffpe = variant_name %in% drivers_called_ffpe)



# Output Figure -----------------------------------------------------------
png(file.path(OUT_DIR, "6_trees.png"), width = 1800, height = 600, res = 150)
layout(matrix(c(1:(length(trial_ids))), nrow = 3, ncol = 9, byrow = T))
for(i in trial_ids){
  par(mar = c(2, 2, 2, 1))
  plot_pyclone_tree(trees_list[[which(names(trees_list) == i)]], 1, -1.175, i, driver_line_width = 0.1)
}
dev.off()

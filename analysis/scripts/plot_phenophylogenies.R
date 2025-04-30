# Create a figure showing the selected genes by CCF

rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(CONIPHER)
library(igraph)

# PATHS -------------------------------------------------------------------

BASE = here::here()
BASE = "/Users/hanleyb/Dropbox (The Francis Crick)/HoLSTF_Breast/Github_Repo"
OUT_DIR = file.path(BASE, "analysis", "figures")
DRIVERS_LOC = file.path(BASE, "data", "variants", "drivers.txt")


# LOAD DATA ---------------------------------------------------------------
drivers = read.delim(DRIVERS_LOC)

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "conipher_helper_functions.R"))
source(file.path(BASE, "src", "custom_tree_plotting.R"))


# LOAD TREES --------------------------------------------------------------
tree_locations = list.files(file.path(BASE, "data", "trees"), full.names = T, pattern = "Ki67.RDS")
hf_num = substr(basename(tree_locations), 1,5)

trees_list_ki67 = lapply(hf_num, function(hf){
  loc = which(grepl(hf, tree_locations))
  location = tree_locations[loc]
  tree = readRDS(location)
  return(tree)
})


# GROUP TREES -------------------------------------------------------------
NETless30 = c("HF041", "HF283", "HF388")
recurrence = "HF057"
NETgreater30 = c("HF039", "HF299", "HF420", "HF105")
order_hf = c(NETless30, recurrence, NETgreater30)


# WRITE TREES -------------------------------------------------------------
png(file.path(OUT_DIR, "Figure5C_phenophylogenies_Ki67.png"), width = 1800, height = 600, res = 150)
layout(matrix(c(1:length(order_hf)), nrow = 2, ncol = 4, byrow = T))
for(i in 1:length(hf_num)){
  par(mar = c(2, 2, 2, 1))
  plot_pyclone_tree_Ki67(trees_list_ki67[[which(hf_num == order_hf[i])]], 1, -1.175, order_hf[i])
}
dev.off()

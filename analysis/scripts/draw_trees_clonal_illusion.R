# Draw the phylogenies of trees where matched SingReg and metastatic samples were available

 
rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(CONIPHER)
library(igraph)

# PATHS -------------------------------------------------------------------
BASE = here::here()
BASE = "/Users/hanleyb/Documents/GitHub/VAULT_Breast/"
OUT_DIR = file.path(BASE, "analysis", "figures")
DRIVERS_LOC = file.path(BASE, "data", "variants", "drivers.txt")


# LOAD DATA ---------------------------------------------------------------
drivers = read.delim(DRIVERS_LOC)

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "conipher_helper_functions.R"))
source(file.path(BASE, "src", "custom_tree_plotting.R"))

# GROUP TREES -------------------------------------------------------------
sr_trees = c("HF052", "HF168", "HF199", "HF255") # single region captures metastasising subclone
trees_repseq = c("HF044", "HF390", "HF414", "HF291", "HF268", "HF170", "HF299")
miss_trees = c("HF270", "HF132", "HF068", "HF032", "HF195", "HF033") #both sample types fail to capture metastasising subclone
trees_repseq_ck = c("HF044", "HF390", "HF414", "HF291", "HF268", "HF170") #RepSeq captures metastasising subclone
all = c(trees_repseq, sr_trees, miss_trees)

# LOAD TREES --------------------------------------------------------------
tree_locations = list.files(file.path(BASE, "data", "trees"), full.names = T)

tree_locations = grep("Ki67", tree_locations, invert = T, value = T) # get rid of the Ki67 trees
tree_locations = grep(paste(all, collapse = "|"), tree_locations, invert = T, value = T) # get rid of the Ki67 trees
all = substr(basename(tree_locations), 1, 5)

trees_list = lapply(all, function(hf){
  loc = which(grepl(hf, tree_locations))
  location = tree_locations[loc]
  tree = readRDS(location)
  return(tree)
})






# WRITE FIGURE 4A ---------------------------------------------------------
png(file.path(OUT_DIR, "Extended_Data_Figure8_pervasive_clonal_illusion.png"), width = 900, height = 1200, res = 150)
layout(matrix(c(c(1:9), 11, 10, 12), ncol = 3, nrow = 4, byrow = T))
for(i in 1:length(all)){
  par(mar = c(2, 2, 2, 1))
  plot_pyclone_tree(trees_list[[which(all == all[i])]], 1, -1.175, all[i])
}
dev.off()

plot.new()
legend("bottomleft", legend = c("CK_RepSeq", "Ki67_RepSamp","Unenriched_RepSeq", "Primary_SingReg", "Clone not identified"), 
       fill = c("#1F78B4db", "#6A3D9Aff", "#E31A1Cdb","#33A02Cff", "grey"), 
       title = expression(bold("Sample Type")))




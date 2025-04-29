# Create a figure showing the selected genes by CCF

rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(CONIPHER)
library(igraph)

# PATHS -------------------------------------------------------------------
BASE = "/Users/hanleyb/Dropbox (The Francis Crick)/HoLSTF_Breast/Github_Repo"
BASE = here::here()
OUT_DIR = file.path(BASE, "analysis", "figures")
DRIVERS_LOC = file.path(BASE, "data", "variants", "drivers.txt")


# LOAD DATA ---------------------------------------------------------------
drivers = read.delim(DRIVERS_LOC)

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "conipher_helper_functions.R"))
source(file.path(BASE, "src", "custom_tree_plotting.R"))


# LOAD TREES --------------------------------------------------------------
tree_locations = list.files(file.path(BASE, "data", "trees"), full.names = T)

tree_locations = grep("Ki67", tree_locations, invert = T, value = T) # get rid of the Ki67 trees
trees_list = lapply(all, function(hf){
  loc = which(grepl(hf, tree_locations))
  location = tree_locations[loc]
  tree = readRDS(location)
  return(tree)
})
# GROUP TREES -------------------------------------------------------------
sr_trees = c("HF052", "HF168", "HF199", "HF255") # single region captures metastasising subclone
trees_repseq = c("HF044", "HF390", "HF414", "HF291", "HF268", "HF170", "HF299")
miss_trees = c("HF270", "HF132", "HF068", "HF032", "HF195", "HF033") #both sample types fail to capture metastasising subclone
trees_repseq_ck = c("HF044", "HF390", "HF414", "HF291", "HF268", "HF170") #RepSeq captures metastasising subclone
all = c(trees_repseq, sr_trees, miss_trees)





# WRITE FIGURE 4A ---------------------------------------------------------
png(file.path(OUT_DIR, "Figure4A_RepSeq_captures_met_subclone.png"), width = 1800, height = 600, res = 150)
layout(matrix(c(1:length(trees_repseq_ck)), nrow = 1, ncol = 6, byrow = T))
for(i in 1:length(trees_repseq_ck)){
  par(mar = c(2, 2, 2, 1))
  plot_pyclone_tree(trees_list[[which(all == trees_repseq_ck[i])]], 1, -1.175, trees_repseq[i])
}
dev.off()


# WRITE FIGURE 4B ---------------------------------------------------------
png(file.path(OUT_DIR, "Figure4B_SingReg_captures_met_subclone.png"), width = 1800, height = 600, res = 150)
layout(matrix(c(1:length(sr_trees)), nrow = 1, ncol = 4, byrow = T))
for(i in 1:length(sr_trees)){
  par(mar = c(2, 2, 2, 1))
  plot_pyclone_tree(trees_list[[which(all == sr_trees[i])]], 1, -1.175, sr_trees[i])
}
dev.off()

# WRITE FIGURE 4C ---------------------------------------------------------
png(file.path(OUT_DIR, "Figure4C_Missed_met_subclone.png"), width = 1800, height = 600, res = 150)
layout(matrix(c(1:length(miss_trees)), nrow = 1, ncol = 6, byrow = T))
for(i in 1:length(miss_trees)){
  par(mar = c(2, 2, 2, 1))
  plot_pyclone_tree(trees_list[[which(all == miss_trees[i])]], 1, -1.175, miss_trees[i])
}
dev.off()


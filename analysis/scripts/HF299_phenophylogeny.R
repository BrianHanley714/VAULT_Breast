# Draw the pheno-phylogeny for HF299

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


# LOAD TREE ---------------------------------------------------------------
HF299_tree = readRDS(file.path(BASE, "data", "trees", "HF299.tree.RDS"))


# WRITE TREE --------------------------------------------------------------
dev.off()
png(file.path(OUT_DIR, "Figure5D_HF299_phenophylogeny.png"), width = 1000, height = 600, res = 150)
plot_pyclone_tree(HF299_tree,1, -1.175, "HF299")
dev.off()



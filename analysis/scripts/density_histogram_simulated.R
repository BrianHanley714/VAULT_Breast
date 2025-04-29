# Density histogram of mutational cancer cell fraction for simulated data

rm(list = ls(all = TRUE))
# LIBRARIES ---------------------------------------------------------------
library(tidyverse)

# PATHS -------------------------------------------------------------------
run = "run001"
BASE = "/Users/hanleyb/Dropbox (The Francis Crick)/HoLSTF_Breast/Github_Repo"
BASE = here::here()
OUT_DIR = file.path(BASE, "analysis", "figures")
input_dir = file.path(BASE, "data/simulations", run, "out/")


# LOAD DATA ---------------------------------------------------------------
if(file.exists(input_dir)){
  num_mut_out = read.delim(paste0(input_dir, "num_mut_out.txt"))
  clones_out = read.delim(paste0(input_dir, "clone_details.txt"))
  mat_repsampling = read.delim(paste0(input_dir, "simulation_data_repseq.txt"))
  ccf_comparisons_rep = read.delim(paste0(input_dir, "ccf_comparisons_data_repseq.txt"))
  mat_sampling = read.delim(paste0(input_dir, "simulation_data_ffpe.txt"))
  ccf_comparisons_sr = read.delim(paste0(input_dir, "ccf_comparisons_data_sr.txt"))
  mat_repsampling = as_tibble(mat_repsampling)
  mat_sampling = as_tibble(mat_sampling)
}


# GRAPHICAL PARAMETERS ----------------------------------------------------
rs_col = "#54278F"
sr_col = "#33A02Cb7"
gt_col = "grey"



# DRAW PLOT ---------------------------------------------------------------
as_tibble(ccf_comparisons_sr)%>%
  mutate(Sampling = "SingReg")%>%
  bind_rows(as_tibble(ccf_comparisons_rep)%>%mutate(Sampling = "RepSamp"))%>%
  bind_rows(as_tibble(num_mut_out)%>%mutate(ccfs = CCF,
                                            true_ccfs = CCF)%>%dplyr::select(ccfs,true_ccfs, -c(CCF, ID))%>%mutate(Sampling = "Ground Truth"))%>%
  filter(ccfs>0.01)%>%
  mutate(Sampling = factor(Sampling, levels = c("SingReg", "Ground Truth", "RepSamp")))%>%
  group_by(Sampling)%>%
  slice_sample(n=500)%>%
  ggplot(aes(x=ccfs, col = Sampling, fill = Sampling))+
  stat_density(geom = "area", position = "identity", alpha = 0.7, aes(y = ..scaled..))+
  theme_classic(base_size = 20)+
  scale_color_manual(values = c("RepSamp" = rs_col, "SingReg" = sr_col, "Ground Truth" = gt_col))+
  scale_fill_manual(values = c("RepSamp" = rs_col, "SingReg" = sr_col, "Ground Truth" = gt_col))+
  xlab("Cancer Cell Fraction")+
  theme(legend.position = "top",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))+
  #  ggtitle(paste("Sampling for", run))+
  ylab("Mutational density (scaled)")


# WRITE PLOT --------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Figure2C_mutational_density_simulated.png"))     

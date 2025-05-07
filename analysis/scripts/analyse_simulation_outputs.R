# demonstrate simulation output supplementary info

rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(cowplot)
library(scales)
library(ggpubr)

# PATHS -------------------------------------------------------------------
run = "run001"
BASE = here::here()
BASE = "/Users/hanleyb/Documents/GitHub/VAULT_Breast/"
OUT_DIR = file.path(BASE, "analysis", "figures")
input_dir = file.path(BASE, "data/simulations", run, "out/")


# LOAD DATA ---------------------------------------------------------------
if(file.exists(input_dir)){
  ccf_comparisons_rep = read.delim(paste0(input_dir, "ccf_comparisons_data_repseq.txt"))
  ccf_comparisons_sr = read.delim(paste0(input_dir, "ccf_comparisons_data_sr.txt"))
  num_mut_out = read.delim(paste0(input_dir, "num_mut_out.txt"))
  clones_out = read.delim(paste0(input_dir, "clone_details.txt"))
  mat_repsampling = read.delim(paste0(input_dir, "simulation_data_repseq.txt"))
  ccf_comparisons_rep = read.delim(paste0(input_dir, "ccf_comparisons_data_repseq.txt"))
  mat_sampling = read.delim(paste0(input_dir, "simulation_data_ffpe.txt"))
  ccf_comparisons_sr = read.delim(paste0(input_dir, "ccf_comparisons_data_sr.txt"))
  mat_repsampling = as_tibble(mat_repsampling)
  mat_sampling = as_tibble(mat_sampling)
}



# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "rot.lab.R"))

# GRAPHICAL PARAMETERS ----------------------------------------------------
rs_col = "#54278F"
sr_col = "#33A02Cb7"
base_size = 10
point_size = 1
line_col = "darkgrey"


# ORGANISE PLOTS INTO LIST ------------------------------------------------
{
plot_list = list()
i = 1

plot_list[[i]] = as_tibble(mat_repsampling)%>%
  mutate(Sampling = "RepSamp")%>%
  bind_rows(as_tibble(mat_sampling)%>%mutate(Sampling = "SingReg"))%>%
  ggplot(aes(Sampling, num_unique_subclones, col = Sampling))+
  geom_jitter(height=0, size = point_size)+
  theme_classic(base_size = base_size)+
  ylab("Subclones (n)")+
  scale_color_manual(values = c("RepSamp" = rs_col, "SingReg" = sr_col))+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank())
  



i=i+1
plot_list[[i]] = as_tibble(mat_repsampling)%>%
  mutate(Sampling = "RepSamp")%>%
  bind_rows(as_tibble(mat_sampling)%>%mutate(Sampling = "SingReg"))%>%
  mutate(index  = row_number())%>%
  pivot_longer(cols = c(ccf_diff_higher_true, ccf_diff_lower_true), names_to = "diff.ccf", values_to = "ccf_difference")%>%
  ggplot(aes(Sampling, ccf_difference))+
  geom_jitter(height = 0, size = point_size, aes(col = diff.ccf))+
  ylab("CCF error")+
  scale_color_manual(values = c("ccf_diff_higher_true" = "maroon", "ccf_diff_lower_true" = "navy"))+
  theme_classic(base_size = base_size)+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank())


i=i+1
plot_list[[i]] = as_tibble(mat_repsampling)%>%
  mutate(Sampling = "RepSamp")%>%
  bind_rows(as_tibble(mat_sampling)%>%mutate(Sampling = "SingReg"))%>%
  ggplot(aes(Sampling, ks_stat, col = Sampling))+
  geom_boxplot(outlier.alpha = 0)+
  geom_jitter(size = point_size)+
  scale_y_log10()+
  theme_classic(base_size = base_size)+
  ylab("Distance (KS)")+
  scale_color_manual(values = c("RepSamp" = rs_col, "SingReg" = sr_col))+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank())
  


i=i+1
plot_list[[i]] = as_tibble(mat_repsampling)%>%
  mutate(Sampling = "RepSamp")%>%
  bind_rows(as_tibble(mat_sampling)%>%mutate(Sampling = "SingReg"))%>%
  ggplot(aes(Sampling, sp_rho, col = Sampling))+
  geom_boxplot(outlier.alpha = 0)+
  geom_jitter(height = 0, size = point_size)+
  theme_classic(base_size = base_size)+
  theme(axis.title.x = element_blank())+
  scale_color_manual(values = c("RepSamp" = rs_col, "SingReg" = sr_col))+
  theme(legend.position = "none")+
  ylab("Correlation (SR)")+
  theme(legend.position = "none")
  

i=i+1
plot_list[[i]] = 
  as_tibble(ccf_comparisons_rep)%>%
    mutate(Clonal_Illusion = if_else(true_ccfs<0.8 & ccfs >0.8, "Illusion of Clonality", "Not IOC"))%>%
    mutate(Clonal_Illusion = factor(Clonal_Illusion, levels = c("Not IOC", "Illusion of Clonality")))%>%
    ggplot(aes(true_ccfs, ccfs))+
  geom_smooth(method = "lm", col = line_col)+
    geom_point(aes(col = Clonal_Illusion), size = point_size)+
  scale_color_manual(values = c("Not IOC" = rs_col, "Illusion of Clonality" = sr_col))+
    theme_classic(base_size = base_size)+
  theme(legend.position = "none")+
  ylab("CCF (obs)")+
  xlab("CCF (true)")

i=i+1
plot_list[[i]] =   
  as_tibble(ccf_comparisons_sr)%>%
    mutate(Clonal_Illusion = if_else(true_ccfs<0.8 & ccfs >0.8, "Illusion of Clonality", "Not IOC"))%>%
    mutate(Clonal_Illusion = factor(Clonal_Illusion, levels = c("Not IOC", "Illusion of Clonality")))%>%
    ggplot(aes(true_ccfs,ccfs))+
  geom_smooth(method = "lm", col =line_col)+
    geom_point(aes(col = Clonal_Illusion), size = point_size)+
  scale_color_manual(values = c("Not IOC" = rs_col, "Illusion of Clonality" = sr_col))+
    theme_classic(base_size = base_size)+
  theme(legend.position = "none")+
  ylab("CCF (obs)")+
  xlab("CCF (true)")
 

i=i+1
plot_list[[i]] = clones_out%>%
  mutate(clone_order = factor(clone_order, levels = c("second", "third", "fourth", "fifth", "sixth", "seventh", "terminal")))%>%
  filter(!is.na(clone_order))%>%
  ggplot(aes(clone_order, Fitness))+
  geom_boxplot(outlier.alpha=0)+
  geom_jitter(height = 0, size = point_size)+
  theme_classic(base_size = base_size)+
  rot.lab()+
  theme(axis.title.x = element_blank())

i=i+1
plot_list[[i]] = clones_out%>%
  ggplot(aes(ID, Fitness))+
  geom_point(size = point_size)+
  theme_classic(base_size = base_size)+
  xlab("Clone order of emergence")+
  rot.lab()

i=i+1
plot_list[[i]] = clones_out%>%
  ggplot(aes(CCF, Fitness, col = ID))+
  geom_point(size = point_size)+
  theme_classic(base_size = base_size)+
  scale_x_log10()+
  theme(legend.position = "none")+
  rot.lab()




i=i+1
plot_list[[i]] = mat_repsampling%>%
  mutate(Sampling = "RepSamp")%>%
  bind_rows(mat_sampling%>%mutate(Sampling = "SingReg"))%>%
  ggplot(aes(prop_total_subclones, max_fit, col = Sampling))+
  geom_point(size = point_size)+
  theme_classic(base_size = base_size)+
  ylab("Max Fitness")+
  xlab("Subclones %")+
  ylim(0, max(c(mat_repsampling$max_fit+2, mat_sampling$max_fit)))+
  annotate(x = 1, y = max(clones_out$Fitness), geom = "point", col = line_col, size = point_size)+
  scale_color_manual(values = c("RepSamp" = rs_col, "SingReg" = sr_col))+
  theme(legend.position = "none")+
  scale_x_log10()
  
i=i+1
plot_list[[i]] =  mat_repsampling%>%
  mutate(Sampling = "RepSamp")%>%
  bind_rows(mat_sampling%>%mutate(Sampling = "SingReg"))%>%
  ggplot(aes(prop_total_subclones, max_fit_over_0.01, col = Sampling))+
  geom_point(size = point_size)+
  theme_classic(base_size = base_size)+
  ylab("Max Fitness")+
  xlab("Subclones %")+
  ylim(0, max(c(mat_repsampling$max_fit+2, mat_sampling$max_fit)))+
  annotate(x = 1, y = max(clones_out$Fitness[clones_out$CCF >0.01]),col = line_col, geom = "point", size = point_size)+
  scale_color_manual(values = c("RepSamp" = rs_col, "SingReg" = sr_col))+
  theme(legend.position = "none")+
  scale_x_log10()
  


i=i+1
plot_list[[i]]   = mat_repsampling%>%
  mutate(Sampling = "RepSamp")%>%
  bind_rows(mat_sampling%>%mutate(Sampling = "SingReg"))%>%
  ggplot(aes(prop_total_subclones, mean_fit, col = Sampling))+
  geom_point(size = point_size)+
  theme_classic(base_size = base_size)+
  ylab("Mean Fitness")+
  xlab("Subclones %")+
  scale_color_manual(values = c("RepSamp" = rs_col, "SingReg" = sr_col))+
  ylim(0, max(c(mat_repsampling$mean_fit, mat_sampling$mean_fit)))+
  annotate(x = 1, y = mean(clones_out$Fitness), col = line_col, geom = "point", size = point_size)+
  theme(legend.position = "none")+
  scale_x_log10()
  
i=i+1    
plot_list[[i]]   =  mat_repsampling%>%
  mutate(Sampling = "RepSamp")%>%
  bind_rows(mat_sampling%>%mutate(Sampling = "SingReg"))%>%
  ggplot(aes(prop_total_subclones, mean_fit_over_0.01, col = Sampling))+
  geom_point(size = point_size)+
  theme_classic(base_size = base_size)+
  ylab("Mean Fitness")+
  xlab("Subclones %")+
  ylim(0, max(c(mat_repsampling$mean_fit, mat_sampling$mean_fit)))+
  scale_color_manual(values = c("RepSamp" = rs_col, "SingReg" = sr_col))+
  annotate(x = 1, y = mean(clones_out$Fitness[clones_out$CCF >0.01]), geom = "point", col = line_col, size = point_size)+
  #annotate(x = 0.9, y = mean(clones_out$Fitness[clones_out$CCF >0.01])+4, geom = "text", label = "Mean\nWhole\nTumour\n(CCF>0.01)")+
  theme(legend.position = "none")+
  scale_x_log10()


true_ITH  = length(num_mut_out$CCF[num_mut_out$CCF<=0.8])/length(num_mut_out$CCF[num_mut_out$CCF>0.8])
true_ITH_over_0.01 = length(num_mut_out$CCF[num_mut_out$CCF<=0.8 & num_mut_out$CCF>0.01])/length(num_mut_out$CCF[num_mut_out$CCF>0.8])


i=i+1
plot_list[[i]] = mat_repsampling%>%
    mutate(Sampling = "RepSamp")%>%
    bind_rows(mat_sampling%>%mutate(Sampling = "SingReg"))%>%
    add_row(tibble(Sampling = "Whole", ITH = true_ITH))%>%
    mutate(Sampling = factor(Sampling, levels = c("Whole", "RepSamp", "SingReg")))%>%
    ggplot(aes(Sampling, ITH, col = Sampling))+
    geom_jitter(size = point_size, height = 0)+
    theme_classic(base_size = base_size)+
  scale_color_manual(values = c("Whole" = line_col, "RepSamp" = rs_col, "SingReg" = sr_col))+
    #stat_compare_means(label.x = 1.5)+
    ylab("ITH_index")+
    rot.lab()+
  theme(axis.title.x = element_blank(),
        legend.position = "none")
  
i=i+1
plot_list[[i]] =   mat_repsampling%>%
    mutate(Sampling = "RepSamp")%>%
    bind_rows(mat_sampling%>%mutate(Sampling = "SingReg"))%>%
    dplyr::select(Sampling, ITH_over_0.01)%>%
    add_row(tibble(Sampling = "Whole", ITH_over_0.01 = true_ITH_over_0.01))%>%
    mutate(Sampling = factor(Sampling, levels = c("Whole", "RepSamp", "SingReg")))%>%
    ggplot(aes(Sampling, ITH_over_0.01, col = Sampling))+
    geom_jitter(size = point_size, height = 0)+
    theme_classic(base_size = base_size)+
  scale_color_manual(values = c("Whole" = line_col, "RepSamp" = rs_col, "SingReg" = sr_col))+
  theme(axis.title.x = element_blank(),
        legend.position = "none")+
    ylab("ITH_index")+
    rot.lab()

}


# DRAW PLOTS --------------------------------------------------------------
plot_grid(plotlist = plot_list, labels = LETTERS[2:(length(plot_list)+1)], ncol = 3)

# WRITE PLOT --------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Extended_Data_Figure3_simulation_outputs.png"))     



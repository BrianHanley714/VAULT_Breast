rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(qs)
library(openxlsx)
library(ggpubr)
# LOAD DATA ---------------------------------------------------------------
BASE = here::here()
BASE = "/Users/hanleyb/Documents/GitHub/VAULT_Breast"
OUT_DIR = file.path(BASE, "analysis", "figures")
cluster_info_path = file.path(BASE, "data", "Trees", "tree_clusters_info.RDS")
VAULT_CLIN_PATH = file.path(BASE, "data", "metadata", "clinical_data.txt")
VAULT_VARIANTS_PATH = file.path(BASE, "data", "variants", "variant_calls_VAULT.rdata")
CASES_INCL_PATH = file.path(BASE, "data", "metadata", "cases_included.xlsx")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_features.R"))
source(file.path(BASE, "src", "custom_tree_plotting.R"))
source(file.path(BASE, "src", "conipher_helper_functions.R"))

get_ccf_per_sample = function(check_sample_type, ffpe_sample_index){
  
  
  sapply(1:length(cluster_locations_list), function(fn){
    
    cluster_info = cluster_locations_list[[fn]]
    unique_samples = unique(cluster_info$SAMPLE)
    if(check_sample_type== "FFPE"){
      FFPE_samples = unique_samples[unique_samples %in% c(ffpe_sampl_info%>%dplyr::filter(Sample_Type %in% c("PRIMARY", "RECURRENCE_RS"))%>%pull(Sample_ID))]
      RS_samples = grep("ET", unique_samples, value = T)
      
      
      clusters_in_ffpe = cluster_info%>%
        dplyr::filter(SAMPLE %in% FFPE_samples[ffpe_sample_index])%>%dplyr::filter(clonality!= "absent")%>%pull(clusterID)
      
      
      clusters_absent = cluster_info%>%
        dplyr::filter(!SAMPLE %in% FFPE_samples[ffpe_sample_index])%>%
        dplyr::filter(clusterID %in% clusters_in_ffpe,
                      clonality == "absent")%>%
        pull(clusterID)%>%table()
      
      clusters_limited_to_sample = names(clusters_absent)[clusters_absent == length(unique_samples)-1]
      
      
      mean_ccf = cluster_info%>%
        dplyr::filter(clusterID%in% clusters_limited_to_sample,
                      SAMPLE == FFPE_samples[ffpe_sample_index])%>%pull(meanCCF)
    }
    
    if(check_sample_type== "RepSeq"){
      FFPE_samples = unique_samples[unique_samples %in% c(ffpe_sampl_info%>%dplyr::filter(Sample_Type %in% c("PRIMARY", "RECURRENCE_RS"))%>%pull(Sample_ID))]
      RS_samples = grep("ET", unique_samples, value = T)
      
      
      clusters_in_rs = cluster_info%>%
        dplyr::filter(SAMPLE %in% RS_samples[ffpe_sample_index])%>%dplyr::filter(clonality!= "absent")%>%pull(clusterID)
      
      
      clusters_absent = cluster_info%>%
        dplyr::filter(!SAMPLE %in% RS_samples[ffpe_sample_index])%>%
        dplyr::filter(clusterID %in% clusters_in_rs,
                      clonality == "absent")%>%
        pull(clusterID)%>%table()
      
      clusters_limited_to_sample = names(clusters_absent)[clusters_absent == length(unique_samples)-1]
      
      
      mean_ccf = cluster_info%>%
        dplyr::filter(clusterID%in% clusters_limited_to_sample,
                      SAMPLE == RS_samples[ffpe_sample_index])%>%pull(meanCCF)
    }
    max = if_else(max(mean_ccf) == -Inf, NA, max(mean_ccf))
    return(max)
  }
  )
}

# LOAD DATA ---------------------------------------------------------------
cluster_locations_list = readRDS(cluster_info_path)



# ORGANISE DATA -----------------------------------------------------------



table(is.na(c(get_ccf_per_sample(check_sample_type = "FFPE", ffpe_sample_index = 1))))
large_ffpe_subclones = c(get_ccf_per_sample(check_sample_type = "FFPE", ffpe_sample_index = 1),
                         get_ccf_per_sample(check_sample_type = "FFPE", ffpe_sample_index = 2)
)

large_ffpe_subclones = large_ffpe_subclones[!is.na(large_ffpe_subclones)]
table(is.na(get_ccf_per_sample(check_sample_type = "RepSeq", ffpe_sample_index = 1)))
large_rs_subclones = get_ccf_per_sample(check_sample_type = "RepSeq", ffpe_sample_index = 1)
large_rs_subclones = large_rs_subclones[!is.na(large_rs_subclones)]
t.test(large_rs_subclones, large_ffpe_subclones)

all_clusters = do.call("rbind", cluster_locations_list)

primary_FFPEs = all_clusters%>%
  dplyr::filter(truncal == FALSE)%>%
  mutate(group = paste(substr(SAMPLE, 1, 5), clusterID, sep = "_"))%>%
  dplyr::mutate(PRIMARY = SAMPLE %in% c(ffpe_sampl_info%>%dplyr::filter(Sample_Type %in% c("PRIMARY", "RECURRENCE_RS"))%>%pull(Sample_ID)))%>%
  dplyr::filter(PRIMARY)%>%
  dplyr::select(nMuts, meanCCF, CCF_CI_high, CCF_CI_low, group)



RS_samples = all_clusters%>%
  dplyr::filter(truncal == FALSE)%>%
  mutate(group = paste(substr(SAMPLE, 1, 5), clusterID, sep = "_"))%>%
  dplyr::filter(grepl("_ET", SAMPLE))%>%
  dplyr::select(nMuts, meanCCF, CCF_CI_high, CCF_CI_low, group)

out_df = left_join(primary_FFPEs, RS_samples, by = "group", suffix = c(".FFPE", ".RS"))  

out_df%>%
  dplyr::filter(meanCCF.FFPE >50 | meanCCF.RS >50,
                !(meanCCF.FFPE >50 & meanCCF.RS >50),
                (abs(meanCCF.FFPE -meanCCF.RS )>30)
  )%>%
  pivot_longer(cols = c(meanCCF.FFPE, meanCCF.RS))%>%
  mutate(name = if_else(name == "meanCCF.FFPE", "SOC", "RepSeq"),
         name = factor(name, levels = c("SOC", "RepSeq"))
         )%>%
  mutate(value = if_else(value >100, 100, value))%>%
  ggplot(aes(name, value))+
  
  geom_violin(aes(fill = name))+
  geom_line(aes(group = group))+
  geom_boxplot(outlier.alpha = 0, width = 0.1)+
  stat_compare_means(paired = T, 
                     label = "p.signif", size = basesize/2,
                     label.x.npc = 0.5,
                     label.y.npc = 0.9
  )+
  #geom_text(aes(label = substr(group, 1, 5)), position = position_jitter(), size = 2)+
  theme_classic(base_size = basesize)+
  scale_fill_manual(values = c(sr_col, rs_col))+
  theme(legend.position = "none",
        axis.title.x = element_blank())+
  ylab("Cancer Cell Fraction")


# Write plot --------------------------------------------------------------


ggsave(file.path(OUT_DIR, "5j_assess_for_IOC.png"))

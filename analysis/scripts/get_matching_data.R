
rm(list = ls(all = TRUE))
# LIBRARIES ---------------------------------------------------------------
library(MatchIt)
library(tidyverse)
library(cowplot)

# PATHS -------------------------------------------------------------------
BASE = here::here()
BASE = "/Users/hanleyb/Documents/GitHub/VAULT_Breast/"
OUT_DIR = file.path(BASE, "analysis", "figures")
MB_PATH = file.path(BASE, "data","metadata", "whole_METABRIC_metadata.txt")
TCGA_PATH = file.path(BASE, "data","metadata", "whole_TCGA_Breast_metadata.txt")
VAULT_PATH = file.path(BASE, "data","metadata", "clinical_data.txt")

# LOAD DATA ---------------------------------------------------------------
metabric = read.delim(MB_PATH)
tcga = read.delim(TCGA_PATH)
vault = read.delim(VAULT_PATH)



# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_features.R"))

# plotting function for figure
plot_change = function(test){
  if(class(comparator_df%>%pull(test))=="character"){
    plot = comparator_df%>%
      group_by(study, .data[[test]])%>%
      count()%>%
      group_by(study)%>%
      mutate(sum = sum(n),
             prop = n/sum)%>%
      mutate(label = if_else(duplicated(paste(sum, study)), NA, sum))%>%
      ggplot(aes(study, prop, fill = .data[[test]]))+
      geom_col()+
      theme_classic()+
      theme(legend.position = c(0.5,0.90),
            legend.margin = margin(t = -5),  
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            # plot.margin = margin(t =5, r = 5, b = 5, l = 5),      
            # legend.box.margin = margin(0, 0, 0, 0),
            axis.title.x = element_blank())+
      scale_fill_manual(values = colours)+
      guides(fill = guide_legend(nrow = 1))+
      geom_text(aes(x = study, label = label, y = 0.5))+
      ylab(paste(test, "(%)"))+
      labs(fill = NULL)+
      
      scale_y_continuous(breaks = c(0, 0.5, 1,1.2))+
      ylim(0,1.2)
    return(plot)
  }
  
  if(class(comparator_df%>%pull(test))=="numeric"){
    plot = comparator_df%>%
      ggplot(aes(study, .data[[test]]))+
      geom_boxplot()+
      theme_classic()
    return(plot)
  }
  else{print(paste("cannot render for data of class", class(comparator_df%>%pull(test))))}
  
}


# function to run propensity based matching
run_matching = function(df_chosen, covariates, covariates_HF, ratio_max){
  df =df_chosen
  a = df%>%
    filter(!(stage == 0))%>%
    dplyr::select(covariates)%>%
    mutate(Trial = 0)
  
  
  b = vault%>%
    dplyr::select(covariates_HF)%>%
    mutate(Trial = 1)
  
  names(b) = names(a) 
  
  
  combined_df = as.data.frame(rbind(a, b))
  
  
  formula <- as.formula(paste("Trial", "~", paste(covariates[2:(length(covariates))], collapse = "+")))
  
  get_distance = function(ratio_t_c, covariates){
    formula <- as.formula(paste("Trial", "~", paste(covariates[2:(length(covariates))], collapse = "+")))
    m.out1 = matchit(formula, data = combined_df, method = "optimal", distance = "glm", replace = F, ratio = ratio_t_c)
    a=summary(m.out1)
    a=as.data.frame(a$sum.matched)
    return(mean(abs(a$`Std. Mean Diff.`[2:nrow(a)])))
  }
  
  vector = c()
  for(i in  1:round(sum(combined_df$Trial == 0)/sum(combined_df$Trial == 1))){
    vector[i] = get_distance(i, covariates)
    names(vector)[i] = i
  }
  
  
  ratio_chosen =which.min(vector)
  
  m.out1 = matchit(formula, data = combined_df, method = "optimal", distance = "glm", replace = F, ratio = as.numeric(ratio_chosen))
  
  
  matched_data = match.data(m.out1)
  matched_data = matched_data%>% filter(!(grepl("HF", PATIENT_ID)))%>%mutate(match = "after")
  return(matched_data)
}


# get ICD10 codes for VAULT
vault$ICD10 = case_when(vault$Histological.Subtype == "IDC_NST" ~ "8500/3",
                        vault$Histological.Subtype == "ILC"~ "8520/3" ,
                        vault$Histological.Subtype == "MixedILC/IDC"~"8522/3",
                        vault$Histological.Subtype == "MucinousCarcinoma"~"8480/3")

# ORGANISE DATA -----------------------------------------------------------

vault$subtype = case_when(vault$Receptor.Subtype == "ER+/HER2-" ~"Luminal",
          vault$Receptor.Subtype == "HER2+" ~"HER2-overexpressing",
          vault$Receptor.Subtype == "TN" ~"TN")

tcga$AGE

metabric$Gender = if_else(metabric$Gender == FALSE, "F", NA)


tcga_before = tcga%>%
  filter(!is.na(stage) & !is.na(molsubtype) & !is.na(Gender) & !is.na(ICD0))%>%
  filter(ICD0 %in% unique(vault$ICD10))

metabric_before = metabric%>%  
  filter(!is.na(stage) & !is.na(molsubtype) & !is.na(Gender) & !is.na(ICD0), !(is.na(Histological.grade)))%>%
  filter(ICD0 %in% unique(vault$ICD10))%>%
  filter(Histological.grade %in% unique(vault$Histological.grade))






covariates_HF = c("Trial.ID", "subtype", "Gender", "stage", "Age", "ICD10")
covariates = c("PATIENT_ID", "molsubtype", "Gender", "stage", 'AGE', "ICD0")


TCGA_after = run_matching(tcga_before, covariates, covariates_HF, "high")
TCGA_after = as.data.frame(TCGA_after)
covariates = c("PATIENT_ID", "molsubtype", "Gender", "stage", 'AGE', "ICD0", "Histological.grade")
covariates_HF = c("Trial.ID", "subtype", "Gender", "stage", "Age", "ICD10", "Histological.grade")
Metabric_after = run_matching(metabric_before, covariates, covariates_HF, "high")

names(covariates_HF) = covariates
TCGA_after$Histological.grade = NA

# CREATE COMPARATOR DATAFRAME --------------------------------------------- 

comparator_df = bind_rows(
vault%>%dplyr::select(covariates_HF)%>%mutate(study = "VAULT"),
tcga%>%dplyr::select(covariates)%>%mutate(study = "TCGA_before"),
metabric%>%dplyr::select(covariates)%>%mutate(study = "MB_before"),
TCGA_after%>%dplyr::select(covariates)%>%mutate(study = "TCGA_after"),
Metabric_after%>%dplyr::select(covariates)%>%mutate(study = "MB_after")
)
comparator_df = comparator_df%>%
  select(everything(), "ICD10" = ICD0)
comparator_df$study = factor(comparator_df$study, levels = c("MB_before", "TCGA_before", "VAULT", "MB_after", "TCGA_after"))
comparator_df$stage = as.character(comparator_df$stage)
comparator_df$ICD10 = if_else(comparator_df$ICD10 %in% unique(vault$ICD10), comparator_df$ICD10, "other")




# DRAW PLOT ---------------------------------------------------------------


plot_grid(
plot_change("Histological.grade"),
plot_change("stage"),
plot_change("ICD10"),
plot_change("molsubtype"),
plot_change("Gender"),
plot_change("AGE"),
ncol = 1
)


# WRITE PLOT --------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Extended_Data_Figure4_matching_TCGA_MB_VAULT.png"))




# WRITE MATCHED PATIENTS --------------------------------------------------
write.table(bind_rows(Metabric_after, TCGA_after), file.path(BASE, "data", "metadata", "matched_patients_characteristics.txt"), sep = "\t", quote = F, row.names = F, col.names = T)



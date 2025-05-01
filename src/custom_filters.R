# filter for repseq CK enriched cases
filter_rs = function(dataframe){dataframe%>%filter(if_else(study == "VAULT", Tumor_Sample_Barcode %in% rs_patients, TRUE))}


# filter for matched cohort in TCGA and METABRIC
filter_matched = function(dataframe){dataframe%>%filter(if_else(study != "VAULT", Tumor_Sample_Barcode %in% matched_patients, TRUE))}

#filter for VAULT Rx naive
filter_VAULT_RXnaive = function(dataframe){dataframe%>%filter(if_else(study == "VAULT", substr(Tumor_Sample_Barcode, 1, 5) %in% c(clinical_data$Trial.ID[is.na(clinical_data$neo_adj_treatment_class) & clinical_data$recurrence_current_sample_is_recurence == "No"]), TRUE))}


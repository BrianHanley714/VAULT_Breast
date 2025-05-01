annotate_driver_summary = function(dataframe){dataframe%>%
    mutate(ONCOGENIC = if_else(grepl("Onco", ONCOGENIC_ONCOKB)|grepl("annotated", CGI.Oncogenic.Summary), ONCOGENIC_ONCOKB, "Non_Driver"))%>%
    mutate(Biomarker_oncoKB = if_else(HIGHEST_LEVEL_ONCOKB!= "" &!is.na(HIGHEST_LEVEL_ONCOKB), "Biomarker", ONCOGENIC),
           Biomarker_oncoKB = if_else(grepl("Non_Driver|Biomarker", Biomarker_oncoKB), Biomarker_oncoKB, "Driver_Other"),
    )%>%
    mutate(ONCOGENICITY = if_else(grepl("Onco", ONCOGENICITY), ONCOGENICITY, "Non_Driver"),
           Biomarker_ACMG = if_else(grepl("1|2", ACTIONABILITY_TIER), "Biomarker", ONCOGENICITY),
           Biomarker_ACMG = if_else(grepl("Non_Driver|Biomarker", Biomarker_ACMG), Biomarker_ACMG, "Driver_Other")
    )%>%
    mutate(Consensus_annotated = if_else(Biomarker_ACMG == "Non_Driver" & Biomarker_oncoKB == "Non_Driver", "Non_Driver",
                                         if_else(Biomarker_ACMG== "Biomarker" & Biomarker_oncoKB != "Biomarker", Biomarker_ACMG,
                                                 if_else(Biomarker_oncoKB== "Biomarker" & Biomarker_ACMG != "Biomarker", Biomarker_oncoKB,
                                                         if_else(Biomarker_oncoKB== "Biomarker" & Biomarker_ACMG == "Biomarker", "Biomarker", "Driver_Other"))))

                                                         )%>%
    #mutate(Consensus_annotated = if_else(Biomarker_ACMG == "Non_Driver" & Biomarker_oncoKB != "Non_Driver", Biomarker_oncoKB, Biomarker_ACMG))%>%
    mutate(rec_type = clinical_data$Receptor.Subtype[match(substr(Tumor_Sample_Barcode, 1, 5), clinical_data$Trial.ID)],
           rec_type = if_else(study == "TCGA", matched_patients_char$molsubtype[match(Tumor_Sample_Barcode, matched_patients_char$PATIENT_ID)], rec_type))%>%
    mutate(trial_specific = if_else(!is.na(MCG_n_trials_all_bc)&MCG_n_trials_all_bc>0|
                                      !is.na(MCG_n_trial_TN)&MCG_n_trial_TN>0 &rec_type == "TN"|
                                      !is.na(MCG_n_trial_HER2pos)&MCG_n_trial_HER2pos>0 &rec_type == "HER2-overexpressing"|
                                      !is.na(MCG_n_trials_ERpos)&MCG_n_trials_ERpos>0 &rec_type == "Luminal", 
                                    "trial_inclusion_criterion",
                                    "non_trial_inclusion"))%>%
    mutate(#actiionable_summary = if_else(Consensus_annotated == "Non_Driver" & trial_specific == "trial_inclusion_criterion", trial_specific, Consensus_annotated),
      actionable_summary = if_else(grepl("Biomarker|Driver_Other", Consensus_annotated) & trial_specific == "trial_inclusion_criterion", paste(Consensus_annotated, trial_specific, sep = ";"), Consensus_annotated))}


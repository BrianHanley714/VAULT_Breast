#devtools::install_github("im3sanger/dndscv")
library(mutSignatures)

library(tidyverse)
library(qs)
library(ggnewscale) 
library(BSgenome.Hsapiens.UCSC.hg19)
RS_col = "#6A3D9Adb" 
SR_col = "#33A02Cb7"



# PATHS -------------------------------------------------------------------
BASE = "/Users/hanleyb/Documents/GitHub/VAULT_Breast"
OUT_DIR = file.path(BASE, "analysis", "figures")
VAULT_VARIANTS_PATH = file.path(BASE, "data", "variants", "variant_calls_VAULT.rdata")

NIKZAINAL_SIGS_PATH = file.path(BASE, "data", "metadata", "breast_signatures_signal_nikzainal.csv")



nikzainal_sigs = read.csv(NIKZAINAL_SIGS_PATH,check.names = F)
maf_out = qread(VAULT_VARIANTS_PATH)
hg19 <- BSgenome.Hsapiens.UCSC.hg19

maf_out$clonality[is.na(maf_out$clonality)] = "INDETERMINATE"
x = maf_out%>%mutate(Chromosome = paste0("chr", Chromosome))%>%dplyr::filter(Tumor_Sample_Barcode %in% ck_included_cases_vault)

x <- attachContext(mutData = x,
                   chr_colName = "Chromosome",
                   start_colName = "Start_Position",
                   end_colName = "End_Position",
                   nucl_contextN = 3,
                   BSGenomeDb = hg19)

x <- removeMismatchMut(mutData = x,                  # input data.frame
                       refMut_colName = "Reference_Allele",       # column name for ref base
                       context_colName = "context",  # column name for context
                       refMut_format = "N")   


x = x%>%dplyr::select(-mutType)%>%as.data.frame%>%
  dplyr::filter(Variant_Type == "SNP")


x <- attachMutType(mutData = x,                      # as above
                   ref_colName = "Reference_Allele",              # column name for ref base
                   var_colName = "Tumor_Seq_Allele2",              # column name for mut base
                   context_colName = "context") 








rownames(nikzainal_sigs) = nikzainal_sigs$Signature
nikzainal_sigs = nikzainal_sigs[,3:length(colnames(nikzainal_sigs))]
nikzainal_sigs = t(nikzainal_sigs)
nikzainal_sigs = nikzainal_sigs[,colnames(nikzainal_sigs) %in% paste0("SBS", c(1, 5, 2, 18, 8, 13, 3, 127))] # these are the sbs signatures seen in >10% of breast cancers according to https://signal.mutationalsignatures.com/
cosmic_signatures = mutSignatures::getCosmicSignatures()

cosmic_signatures_df = mutSignatures::as.data.frame(cosmic_signatures)


nikzainal_sigs = nikzainal_sigs[match(cosmic_signatures@mutTypes$mutTypes, rownames(nikzainal_sigs)),]
rownames(nikzainal_sigs) ==rownames( cosmic_signatures_df )


#cosmic_signatures_df = mutSignatures::as.data.frame(cosmic_signatures)
#cosmic_signatures_df = cosmic_signatures_df[,!colnames(cosmic_signatures_df) %in% c("COSMIC.1", "COSMIC.6")]
#cosmic_signatures_df_ffpesig = cbind(cosmic_signatures_df, ffpe_sig_rep = ffpe_sig_match[,2])

plot_output = function(test_includes_SBS6){
if(test_includes_SBS6){
nikzainal_sigs_ffpe = as.data.frame(cbind(nikzainal_sigs, SBS6 = cosmic_signatures_df[,colnames(cosmic_signatures_df) == "COSMIC.6"]))
mut.counts = countMutTypes(mutTable = x,
                           mutType_colName = "mutType",
                           sample_colName = "SampleID")
} else{
  nikzainal_sigs_ffpe = as.data.frame(nikzainal_sigs)
  filtered_x = x%>%dplyr::filter(Tumor_Sample_Barcode != "HF057_ET001_D01")
  mut.counts = countMutTypes(mutTable = filtered_x,
                             mutType_colName = "mutType",
                             sample_colName = "SampleID")
}
  out = resolveMutSignatures(mut.counts, mutSignatures::as.mutation.signatures(nikzainal_sigs_ffpe))
out$Results
msigPlot(out$Results$count.result)

out_df = mutSignatures::as.data.frame(out$Results$count.result)



levels =maf_out%>%
  group_by(Tumor_Sample_Barcode)%>%
  reframe(count = n())%>%
  dplyr::filter(Tumor_Sample_Barcode %in% ck_included_cases_vault,
                count>50)%>%
  arrange(count)%>%
  mutate(Tumor_Sample_Barcode = substr(Tumor_Sample_Barcode, 1, 5))%>%
  pull(Tumor_Sample_Barcode)



print(plot_grid(
  maf_out%>%
    group_by(Tumor_Sample_Barcode)%>%
    reframe(count = n())%>%
    dplyr::filter(Tumor_Sample_Barcode %in% ck_included_cases_vault,
                  count>50)%>%
    mutate(Tumor_Sample_Barcode = factor(substr(Tumor_Sample_Barcode, 1, 5), levels = levels)
    )%>%
    ggplot(aes(Tumor_Sample_Barcode, count))+
    geom_col()+theme_classic()+
    scale_y_log10()+
    geom_hline(yintercept = 380, linetype = "dashed")+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank()
    ), 
  maf_out%>%
    #dplyr::filter(Variant_Type == "SNP",grep)%>%
    group_by(Tumor_Sample_Barcode, clonality)%>%
    reframe(count = n())%>%
    dplyr::filter(Tumor_Sample_Barcode %in% ck_included_cases_vault,
    )%>%
    
    group_by(Tumor_Sample_Barcode)%>%
    mutate(sum = sum(count),
           Tumor_Sample_Barcode = substr(Tumor_Sample_Barcode, 1, 5),
           Prop = count/sum)%>%
    dplyr::filter(Tumor_Sample_Barcode %in% levels)%>%
    mutate(Tumor_Sample_Barcode = factor(Tumor_Sample_Barcode, levels = levels))%>%
    
    ggplot()+
    geom_col(aes(Tumor_Sample_Barcode, Prop, fill = clonality))+
    geom_point(aes(Tumor_Sample_Barcode,1.17), alpha = 0)+
    theme_classic()+
    theme(legend.position = c(0.5, 0.96),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    guides(fill = guide_legend(nrow = 1, title = element_blank()))+
    scale_fill_brewer(palette = "Set2", na.value = "grey80")
  
  , 
  
  out_df%>%
    mutate(signatures = rownames(.))%>%
    pivot_longer(c(everything(), -signatures))%>%
    group_by(name)%>%
    mutate(order = sum(value))%>%
    dplyr::filter(name %in% ck_included_cases_vault)%>%
    mutate(name = factor(substr(name, 1, 5), levels = levels))%>%
    dplyr::filter(value !=0)%>%
    group_by(name)%>%
    mutate(prop = value/order)%>%
    dplyr::filter(!is.na(name))%>%
    ggplot()+
    geom_point(aes(name,1.05), alpha = 0)+
    geom_col(aes(name,prop, fill = signatures))+
    #geom_hline(yintercept = 380, linetype = "dashed")+
    theme_classic()+
    theme(axis.title.x = element_blank(),
          legend.position = c(0.5, 0.97))+
    guides(fill = guide_legend(nrow = 1, title = element_blank()))+
    rot.lab()+
    scale_fill_brewer(palette = "Set1")+
    ylab("Prop"),
  ncol = 1,
  rel_heights = c(1, 1.2, 3)
))
return(out_df)
}

out_df = plot_output(test_includes_SBS6 = T)
ggsave(file.path(OUT_DIR, "SupplementaryFigure_Signatures.png"))

out_df = plot_output(test_includes_SBS6 = F)
ggsave(file.path(OUT_DIR, "3G_signatures.png"))


out_df%>%
  mutate(signatures = rownames(.))%>%
  pivot_longer(c(everything(), -signatures))%>%
  group_by(name)%>%
  mutate(order = sum(value))%>%
  dplyr::filter(name %in% ck_included_cases_vault)%>%
  mutate(name = factor(substr(name, 1, 5), levels = levels))%>%
  dplyr::filter(value !=0)%>%
  group_by(name)%>%
  mutate(prop = value/order)%>%
  dplyr::filter(!is.na(name))


stat_df = out_df  
stat_df = sweep((stat_df), 2, colSums(stat_df, na.rm = TRUE), FUN = "/")
stat_df = t(stat_df)%>%as.data.frame()
rowSums(stat_df)
stat_df$TMB = rowSums(t(out_df))

stat_df = stat_df%>%
  mutate(clocklike_sig = SBS1+SBS5,
         ABOBEC = SBS13+SBS2
  )

summary(lm(SBS13 ~TMB, stat_df[!rownames(stat_df) %in% c("HF079_ET001_D01", "HF057_ET001_D01"),]))
summary(lm(SBS2 ~TMB, stat_df[!rownames(stat_df) %in% c("HF079_ET001_D01", "HF057_ET001_D01"),]))
summary(lm(SBS3 ~TMB, stat_df[!rownames(stat_df) %in% c("HF079_ET001_D01", "HF057_ET001_D01"),]))
summary(lm(SBS18 ~TMB, stat_df[!rownames(stat_df) %in% c("HF079_ET001_D01", "HF057_ET001_D01"),]))
summary(lm(SBS8 ~TMB, stat_df[!rownames(stat_df) %in% c("HF079_ET001_D01", "HF057_ET001_D01"),]))
summary(lm(SBS5 ~TMB, stat_df[!rownames(stat_df) %in% c("HF079_ET001_D01", "HF057_ET001_D01"),]))
summary(lm(SBS1 ~TMB, stat_df[!rownames(stat_df) %in% c("HF079_ET001_D01", "HF057_ET001_D01"),]))
summary(lm(clocklike_sig ~TMB, stat_df[!rownames(stat_df) %in% c("HF079_ET001_D01", "HF057_ET001_D01"),]))
summary(lm(ABOBEC ~TMB, stat_df[!rownames(stat_df) %in% c("HF079_ET001_D01", "HF057_ET001_D01"),]))


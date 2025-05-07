# Consort Diagram

rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------
library(DiagrammeR)
library(glue)
library(DiagrammeRsvg)
library(rsvg)

# LOAD DATA ---------------------------------------------------------------
BASE = here::here()
BASE = "/Users/hanleyb/Documents/Github/VAULT_Breast"
OUT_DIR = file.path(BASE, "analysis", "figures")


# ASSIGN VARIABLES --------------------------------------------------------
total_breast_surgeries = 5508
total_cancer_surgeries = 3317
non_cancer_surgeries = total_breast_surgeries - total_cancer_surgeries
total_invcancer_surgeries = 3000
tum_remains_n = 410
tum_remains_prop = round((tum_remains_n / total_invcancer_surgeries)*100, digits = 1)
ci = prop.test(tum_remains_n, total_invcancer_surgeries)
ci_lower = round(ci$conf.int[1]*100, digits = 1)
ci_upper = round(ci$conf.int[2]*100,digits = 1)
vault_recruited = 85
dilute_sample = 5
bacterial_contamination = 1
non_carcinoma = 2
small_tumour = 1
no_tumour_remains = 1
low_tumour_purity = 3
formalin_over_fixation = 2
excluded_col = "mistyrose"
vault_sorted = 75
vault_reported = 70
snv_indel_only = 1
snv_indel_cn = 69


# DRAW PLOT ---------------------------------------------------------------

out = grViz(glue("
digraph consort {{
  graph [rankdir = TB, layout = dot]

  node [shape = box, style = filled, fillcolor = lavender, fontname = Helvetica]

  Start [label=<Total Breast Surgeries<br/> 10/09/2018 – 09/03/2022 <br/> (n=<b><u>{total_breast_surgeries}</u></b> )>]
  NON_CANCER [label=<<u>Excluded</u> <br/>- Non-Cancer Surgeries <br/>(n=<b><u>{non_cancer_surgeries}</u></b>)>, fillcolor = {excluded_col}]
  CANCER [label=<Cancer Surgeries<br/> (n=<b><u>{total_cancer_surgeries}</u></b>)>]
  
  NON_INVASIVE [label=<<u>Excluded</u><br/>-Carcinoma In-Situ (n=<b><u>500</u></b>)<br/>-Borderline Tumours (n=<b><u>500</u></b>)>, fillcolor = {excluded_col}]
  INVASIVE [label=<Invasive Cancer Surgeries<br/>(n=<b><u>{total_invcancer_surgeries}</u></b>)>]  
  TUMOUR_REMAINS [label=<Tumour-Remains
    <br/>(n=<b><u>{tum_remains_n}</u></b>)>]
  VAULT_RECRUITED [label=<Recruited to VAULT<br/>(n=<b><u>{vault_recruited}</u></b>)>]
  TUMOUR_REMAINS_PROP [label=<Tumour-Remains (%)
    <br/><b><u>{tum_remains_prop}%</u></b>  (95% CI {ci_lower}-{ci_upper}%)>]
  NOT_SORTED [label=<<u>Excluded </u>
    <br/>-Dilute Sample (n=<b><u>{dilute_sample}</u></b>) 
    <br/>-Bacterial Contamination (n=<b><u>{bacterial_contamination}</u></b>)
    <br/>-Non-Carcinoma Histology (n=<b><u>{non_carcinoma}</u></b>)
    <br/>-Small invasive tumour (7mm) (n=<b><u>{small_tumour}</u></b>)
    <br/>-No tumour remains (n=<b><u>{no_tumour_remains = 1}</u></b>)>, fillcolor = {excluded_col}]
  SORTED [label=<Fixed FACS and WES<br/>(n=<b><u>{vault_sorted}</u></b>)>]
  GENOMICS_FAILED [label=<<u>Excluded</u><br/>
    -Low Tumour Purity (n=<b><u>{low_tumour_purity}</u></b>)<br/>
    -Over 175 days in formalin (n=<b><u>{formalin_over_fixation}</u></b>)>, fillcolor = {excluded_col}]
  IN_STUDY [label=<Genomics Reported<br/>(n=<b><u>{vault_reported}</u></b>)>]
  SNV_INDEL_ONLY [label=<SNVs and INDELS only <br/>(n=<b><u>{snv_indel_only}</u></b>)>]
  SNV_INDEL_CN [label=<SNVs,INDELS and COPY NUMBER <br/>(n=<b><u>{snv_indel_cn}</u></b>)>]

  Start -> NON_CANCER
  Start -> CANCER
  CANCER -> NON_INVASIVE  
  CANCER -> INVASIVE
  INVASIVE -> TUMOUR_REMAINS
  TUMOUR_REMAINS -> VAULT_RECRUITED
  TUMOUR_REMAINS -> TUMOUR_REMAINS_PROP
  VAULT_RECRUITED -> NOT_SORTED
  VAULT_RECRUITED -> SORTED
  SORTED -> GENOMICS_FAILED
  SORTED -> IN_STUDY
  IN_STUDY -> SNV_INDEL_ONLY
  IN_STUDY -> SNV_INDEL_CN
}}
")
)


# WRITE PLOT --------------------------------------------------------------
svg <- export_svg(out)
rsvg_png(charToRaw(svg), file = file.path(OUT_DIR, "Extended_Data_Figure1_consort.png"), width = 800, height = 600)



###############################################################################
######################### Prenatal Exposome Analysis ##########################
###############################################################################
setwd("./TFM")
exposome <- read.csv2("./db/exposome/data_HELIXsubcohort_18_09_23_FILTERED_F.csv")
getwd()
vars <- c(
  # Air pollution   # Traffic load
  "h_no2_ratio_preg_Log",
  "h_pm10_ratio_preg_None",
  "h_pm25_ratio_preg_None",
  "h_trafload_preg_pow1over3",
  "h_trafnear_preg_pow1over3",
  
  # Built environment   # Natural spaces
  "h_builtdens300_preg_Sqrt",
  "h_fdensity300_preg_Log",
  "h_frichness300_preg_None",
  "h_landuseshan300_preg_None",
  "h_popdens_preg_Sqrt",
  "h_walkability_mean_preg_None",
  "h_ndvi100_preg_None",

  # Sociodemographic
  "h_age_None",
  "e3_gac_None",
  "e3_asmokyn_p_None",
  "e3_alcpreg_yn_None",
  "h_mbmi_None",
  "hs_wgtgain_None",
  "h_edumc_None",
  
  # Metals
  "hs_cd_m_Log2",
  "hs_hg_m_Log2",
  "hs_pb_m_Log2",
  
  # OCs
  "hs_dde_madj_Log2",
  "hs_hcb_madj_Log2",
  "hs_pcb138_madj_Log2",
  "hs_pcb153_madj_Log2",
  "hs_pcb180_madj_Log2",
  
  # PFAS
  "hs_pfhxs_m_Log2",
  "hs_pfna_m_Log2",
  "hs_pfoa_m_Log2",
  "hs_pfos_m_Log2",
  
  # Diet
  "h_fish_preg_Ter",
  "h_fruit_preg_Ter",
  "h_legume_preg_Ter",
  "h_veg_preg_Ter",
  "h_dairy_preg_Ter",
  "h_meat_preg_Ter"
)

exposome <- exposome[c(vars,"X.7")]
exposome$e3_asmokyn_p_None <- ifelse(exposome$e3_asmokyn_p_None == "yes", 1, 0)

rownames(exposome) <- exposome$X.7
exposome$X.7 <- NULL
write.csv2(exposome, "./results/ExWAS/exposome_filtered.csv", row.names = T)


#
codebook <- read.csv2("./db/exposome/CODEBOOK_ANALYSIS_AUGUSTO.csv")


codebook_fil<-codebook[codebook$Variable_name_TRANS %in% vars,] %>%arrange(Group)
write_xlsx(codebook_fil[6:12], "./results/ExWAS/codebool_fil.xlsx")



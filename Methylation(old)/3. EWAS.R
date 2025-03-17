###############################################################################
################################ EWAS Analysis ################################ 
###############################################################################

# //////// R environment 

rm(list=ls())
setwd("./TFM")

library(bumphunter)
library(limma)
library(minfi)
library(RColorBrewer)
library(matrixStats)
library(minfiData)
library(stringr)
library(plyr)
library(ggplot2)
library(reshape2)
library(GOplot)
library(snpStats)
library(glmnet)
library(fastDummies)
library(dplyr)
library(ggrepel)
library(EnhancedVolcano)
library(tidyverse)
library(RColorBrewer)
library(generics)
library(parsnip)
library(yardstick)
library(readr)
library(readxl)
library(fastDummies)
library(data.table)

# Performances in test set of model refitted in training set
#require(sharp)
# detach("package:sharp",unload= TRUE) # Detach the sharp function (CRAN)
# #Clean version
# # install.packages("./package_versions/clean/sharp-main",
# #                  repos= NULL,
# #                  type = "source", force=TRUE) # Load Sharp function modified to apply weights
# #Modified version
# install.packages("./packages_R/sharp-main",
#                  repos= NULL,
#                  type = "source", force=TRUE) # Load Sharp function modified to apply weights
# # install.packages("sharp",force=TRUE)
# install.packages("./packages_R/sharp-main", repos= NULL, type = "source", force=TRUE) 
library(sharp)

# //////// LOADING OMICS DATASET

# /// Load HelixID by group (N/S)
helixid_n <- c(read.csv("./helixid_n.csv", row.names = 1)$x)
helixid_s <- c(read.csv("./helixid_s.csv", row.names = 1)$x)

#Loading cpgs list
ewas_cat_df<-read_tsv("ewas_catalog_cpgs.tsv")

#Selecting cpgs associated with BP
ewas_cat_df$trait<-tolower(ewas_cat_df$trait)
ewas_cat_df<-ewas_cat_df[ewas_cat_df$trait %in% c("systolic blood pressure", "diastolic blood pressure"),]
ewas_cat_cpgs <- unique(ewas_cat_df$cpg)

#PUBMED: (EWAS AND blood pressure) OR (dna methylation AND blood pressure) OR (EWAS AND hyperten*) OR (dna methylation hyperten*)
### (Wang et al., 2013)
wang_cpgs <- c("cg09772827346", "cg14958635", "cg11719157630", "cg14371590222", "cg27431859691", 
               "cg079623151375", "cg06701500565", "cg15118204191", "cg04463638148")

### (Boström et al., 2016)
bostrom_cpgs <- c("cg00161968", "cg00875989", "cg06251539", "cg08706258", "cg09134341", "cg10146710", 
                  "cg10596925", "cg10640093", "cg12360759", "cg15612682", "cg16076930", "cg16118212", 
                  "cg16500810", "cg18643762", "cg20841073", "cg21344124", "cg21996137", "cg22011370", 
                  "cg22295383", "cg23945265", "cg25521086", "cg25544164")

### (Domínguez-Barragán et.al., 2023)
dguez_cpgs <- c("cg05399785", "cg12728588", "cg03084350", "cg02107842", "cg05575921", "cg27395200", 
                "cg08774868", "cg13518625", "cg23900905", "cg23761815", "cg00574958", "cg02079413", 
                "cg01678580", "cg03819286", "cg02650017", "cg18181703", "cg20761853", "cg00711496")

### (Irvin et al., 2021)
irvin_cpgs <- c("cg14476101", "cg06690548", "cg23999170", "cg00574958", "cg00574958", "cg19693031",
                "cg06690548", "cg07598370", "cg19693031", "cg18120259")

### (Kazmi et al., 2020)
kazmi_cpgs <- c("cg15935121", "cg26332488", "cg07797660", "cg25203007", "cg04427651", "cg06754224",
                "cg00039326", "cg14663208", "cg20146909", "cg03725309", "cg18933331", "cg23999170", 
                "cg16246545", "cg14476101", "cg19693031", "cg19266329", "cg24955196", "cg12593793",
                "cg17453456", "cg00936728", "cg04275362", "cg21777154", "cg21534578", "cg12417775",
                "cg08035323", "cg01243072", "cg11938080", "cg21990144", "cg18119407", "cg22007809", 
                "cg15616915", "cg09639152", "cg00959259", "cg22213445", "cg24960291", "cg22959409", 
                "cg06690548", "cg07990556", "cg05473987", "cg20131596", "cg27054084", "cg22885332",
                "cg04104695", "cg21066063", "cg21618521", "cg18120259", "cg03125341", "cg21429551", 
                "cg03068497", "cg19390658", "cg07621224", "cg22510074", "cg06688763", "cg15261712", 
                "cg22103219", "cg12816198", "cg05632420", "cg00008629", "cg04665046", "cg17443080", 
                "cg03383434", "cg07856667", "cg15995714", "cg02116864", "cg00805360", "cg17061862", 
                "cg03393444", "cg14099685", "cg07160014", "cg11376147", "cg06178669", "cg15920975", 
                "cg20379593", "cg00574958", "cg17058475", "cg20374917", "cg10601624", "cg09680149", 
                "cg06826457", "cg14741228", "cg05242244", "cg06679990", "cg22608507", "cg09001549", 
                "cg00716257", "cg22143352", "cg02946885", "cg00944304", "cg10941749", "cg26916780", 
                "cg07175797", "cg05288253", "cg20805367", "cg13917614", "cg22361181", "cg08857797", 
                "cg14179401", "cg11133963", "cg18824549", "cg02976539", "cg14818621", "cg06599169", 
                "cg11719283", "cg22304262", "cg02711608", "cg21766592", "cg15114651", "cg07626482", 
                "cg00711496", "cg17417856", "cg22052056", "cg24890964", "cg20239391", "cg01385679", 
                "cg10589813", "cg14090647", "cg15333769", "cg07141002", "cg17863679")

bibliography_cpgs <- unique(c(ewas_cat_cpgs, wang_cpgs, bostrom_cpgs, dguez_cpgs, irvin_cpgs, kazmi_cpgs))


# /// Loading data
# ···· Metadata .RData filepath
metadataFile <- "./script/denoising_2024/metadata/HELIX_SVA_common_OmicsMetadata_20231026.RData" 

# ···· Omic data (Rdata file with an ExpressionSet or GenomicRatioSet)
# ··················· INDICAR RDATA NUEVO CON LA WINSORIZACIÓN
omicFile <- "./results/denoising/denoising_methyl/methylome_subcohort_winsor_ICC_filtered.RData"

# ···· Phenotype data 
phenotype <- readRDS("./db/pheno/final/bp_wide_validN5332023-10-16.rds")

rownames(phenotype) <- phenotype$HelixID

names(phenotype)[grep("cohort.x", names(phenotype))] <- "cohort"
names(phenotype)[grep("h_ethnicity_c.x", names(phenotype))] <- "h_ethnicity_c"
names(phenotype)[grep("e3_sex_Time1", names(phenotype))] <- "e3_sex"


# ···· Modify phenotype, cohort and sex, with dummies groups (1 vs. all columns)
phenotype$h_ethnicity_c <- ifelse(phenotype$h_ethnicity_c %in% c("Asian","Pakistani"), "Asian_pakistani", phenotype$h_ethnicity_c)

phenotype <- dummy_cols(phenotype,
                        select_columns = "h_ethnicity_c",
                        remove_selected_columns = TRUE)

phenotype <- dummy_cols(phenotype,
                        select_columns = "cohort",
                        remove_selected_columns = TRUE)

phenotype <- dummy_cols(phenotype,
                        select_columns = "e3_sex",
                        remove_selected_columns = TRUE)

phenotype <- phenotype[,colnames(phenotype)!=c("h_ethnicity_c_NA")]

#Convert var to factor
eth_dummies <- names(phenotype)[grep("h_ethnicity_", names(phenotype))]
phenotype[eth_dummies] <- lapply(phenotype[eth_dummies], as.factor)

coh_dummies <- names(phenotype)[grep("cohort", names(phenotype))]
phenotype[coh_dummies] <- lapply(phenotype[coh_dummies], as.factor)


# /// Main variables definition
# 1. Omic used: Define omics used in analysis omicLayer variable.
#   Possible values are: 
#       Methyl   : Methylation data
#       Gene     : Genomic data
#       Prot     : Proteomic data
#       miRNAs   : micro-RNAs
#       Serum    : Serum methabolome
#       Urine    : Urine metabolome
omicsLayer <- "Methyl"

# 2. ID column name in the phenotype dataset 
#   It relates the phenotype to the Omic data. (by default SampleID)
phenotypeID <- "SampleID"


# 3. Variables to be used from phenotype dataset (by default ALL)
#  Examples:
#       phenoVariables <- c('bmi', 'sex', 'age', 'smoke')  # 
#       phenoVariables <- 'ALL'         # If we need ALL the variables present in phenotype dataset
#       phenoVariables <- ''            # If no value is defined then ALL the variables present in the phenotype will be used


# LAS VARIABLES ASIGNADAS EN LOS DATOS FILTRADOS
# c("hs2_zdia_bp.v3_2017_Time1","hs2_zsys_bp.v3_2017_Time1","hs2_zbmi_who_Time1","hs2_height_c_zscore_Time1","e3_sex_Time1","hs2_visit_age_years_Time1")

phenotypeVariables <- c("ALL")

# 4. Ethnicity - Population
#   Ethnic parameter used to make a subset of the data and also to take in to account
#   the PCs related to the population. 
#   For European (ONLY if we are using EUR population and no other) we use the PCs calculated from the EUR population, 
#   otherwise we use the PCs computed from ALL populations 
#   Define ethnicity to filter
#   Possible values are: 
#       AFR     :  African
#       AMR     :  Ad Mixed American
#       EAS     :  East Asian
#       EUR     :  European
#       SAS     :  South Asian
#       MIXED   :  Mixed population (not classified in the other categories)
#       ALL     : No filter
#       
#  Example:
#       ethnic <- 'ALL'                     # All populations - no filter
#       ethnic <- ''                        # All populations - no filter
#       ethnic <- c('AFR', 'AMR', 'SAS')    # Use African, Mixed American and South Asian populations
# 

ethnic <- c('ALL')

#### STEP 3. Data preparation
source("./script/denoising_2024/functions/generic_functions_denoising_v2.R")
source("./script/denoising_2024/functions/mlmer_local.R")

omics <- getfullOmicsPhenotype(omicFile, metadataFile, phenotype, phenotypeID, phenotypeVariables, omicsLayer, ethnic )    
ids1<-colnames(omics$data)

# Hacemos dummies con Methyl.Array
pData(omics$data) <- DataFrame(dummy_cols(pData(omics$data), 
                                          select_columns = "Methyl.Array",
                                          remove_selected_columns = FALSE), row.names = pData(omics$data)$SampleID)

pData(omics$data) <- DataFrame(dummy_cols(pData(omics$data), 
                                          select_columns = "Methyl.BATCH_Conversion",
                                          remove_selected_columns = FALSE), row.names = pData(omics$data)$SampleID)

columns_to_modify <- c("Neu", "Mono", "Bcell", "CD4T", "CD8T", "NK", "Eos")
pData(omics$data)[columns_to_modify] <- lapply(pData(omics$data)[columns_to_modify], round, 4)

#Scale continuos vars
pData(omics$data)[c(columns_to_modify,"hs2_visit_age_years_Time1")] <- scale(pData(omics$data)[c(columns_to_modify,"hs2_visit_age_years_Time1")],center = TRUE,scale = TRUE)

#Comprovació
colnames(omics$data)<-omics$data$SampleID
all(ids1==colnames(omics$data))

#Filter by group (N/S)
omics_N <- omics$data[,which(pData(omics$data)$HelixID.x %in% helixid_n)]
omics_S <- omics$data[,which(pData(omics$data)$HelixID.x %in% helixid_s)]
omics_all <- omics$data[,which(pData(omics$data)$HelixID.x %in% c(helixid_n,helixid_s))]

imp_data_n <- data.frame(t(getBeta(omics_N)))
imp_data_s <- data.frame(t(getBeta(omics_S)))
imp_data <- data.frame(t(getBeta(omics_all)))

rownames(imp_data_n) <- pData(omics_N)$HelixID.x
rownames(imp_data_s) <- pData(omics_S)$HelixID.x
rownames(imp_data) <- pData(omics_all)$HelixID.x

# /// Load covars
array_dummies <- names(pData(omics_all))[grep("Methyl.Array_", names(pData(omics_all)))]
array_dummies <- array_dummies[-length(array_dummies)]
batch_dummies <- names(pData(omics_all))[grep("Methyl.BATCH_Conversion_", names(pData(omics_all)))]
batch_dummies <- batch_dummies[-length(batch_dummies)]

#delete 1 var for dummies!!!!
lcovars_n <- c(array_dummies, batch_dummies, "CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu", "NK", "Eos",
               "e3_bw", "hs2_visit_age_years_Time1", "e3_sex_male","h_ethnicity_c_Asian_pakistani", 
               "h_ethnicity_c_Caucasian", 'cohort_BIB', "cohort_EDEN", "cohort_KANC")

lcovars_s <- c(array_dummies, batch_dummies, "CD8T", "CD4T" , "NK", "Bcell", "Mono","NK", "Eos",
               "Neu", "e3_bw", "hs2_visit_age_years_Time1", "e3_sex_male", "cohort_SAB")

lcovars_all <- c(array_dummies, batch_dummies, "CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu", "NK", "Eos",
                 "e3_bw", "hs2_visit_age_years_Time1", "e3_sex_male","h_ethnicity_c_Asian_pakistani",
                 "h_ethnicity_c_Caucasian", 'cohort_BIB', "cohort_EDEN", "cohort_KANC", "cohort_MOBA",
                 "cohort_SAB")

outcomes<- c("hs2_zdia_bp.v3_2017_Time1","hs2_zsys_bp.v3_2017_Time1",
             "hs2_zdia_bp.v3_2017_Time2","hs2_zsys_bp.v3_2017_Time2")

covars_n<-data.frame(pData(omics_N)[lcovars_n])
outcome_n<-data.frame(pData(omics_N)[outcomes])
rownames(covars_n)<-pData(omics_N)$HelixID.x
rownames(outcome_n)<-pData(omics_N)$HelixID.x

covars_s<- data.frame(pData(omics_S)[lcovars_s])
outcome_s<-data.frame(pData(omics_S)[outcomes])
rownames(covars_s)<-pData(omics_S)$HelixID.x
rownames(outcome_s)<-pData(omics_S)$HelixID.x

covars_all<-data.frame(pData(omics_all)[lcovars_all])
outcome_all<-data.frame(pData(omics_all)[outcomes])
rownames(covars_all)<-pData(omics_all)$HelixID.x
rownames(outcome_all)<-pData(omics_all)$HelixID.x

#Delete NAs
nas_id_n <- rownames(covars_n[!complete.cases(cbind(covars_n,outcome_n)),])
imp_data_n <- imp_data_n[!rownames(imp_data_n) %in% nas_id_n,]
covars_n <- covars_n[!rownames(covars_n) %in% nas_id_n,]
outcome_n <- outcome_n[!rownames(outcome_n) %in% nas_id_n,]

nas_id <- rownames(covars_all[!complete.cases(cbind(covars_all,outcome_all)),])
imp_data <- imp_data[!rownames(imp_data) %in% nas_id,]
covars_all <- covars_all[!rownames(covars_all) %in% nas_id,]
outcome_all <- outcome_all[!rownames(outcome_all) %in% nas_id,]

# ipws
ipw <- read.csv2("./db/ipws_bp_1feb24.csv")
rownames(ipw) <- ipw$HelixID
ipw <- ipw[ipw$HelixID %in% rownames(imp_data),]
ipw_n<- ipw[ipw$HelixID %in% rownames(imp_data_n),]
ipw_s<- ipw[ipw$HelixID %in% rownames(imp_data_s),]

#Comprovacions
sort_by_rownames <- function(df) {
  df[order(rownames(df)),]}

imp_data<-sort_by_rownames(imp_data)
imp_data_n<-sort_by_rownames(imp_data_n)
imp_data_s<-sort_by_rownames(imp_data_s)
covars_all<-sort_by_rownames(covars_all)
covars_n<-sort_by_rownames(covars_n)
covars_s<-sort_by_rownames(covars_s)
outcome_all<-sort_by_rownames(outcome_all)
outcome_n<-sort_by_rownames(outcome_n)
outcome_s<-sort_by_rownames(outcome_s)
ipw<-sort_by_rownames(ipw)
ipw_n<-sort_by_rownames(ipw_n)
ipw_s<-sort_by_rownames(ipw_s)

all(rownames(covars_all)==rownames(imp_data))
all(rownames(covars_all)==rownames(outcome_all))
all(rownames(covars_n)==rownames(imp_data_n))
all(rownames(covars_n)==rownames(outcome_n))
all(rownames(covars_s)==rownames(imp_data_s))
all(rownames(covars_s)==rownames(outcome_s))
all(rownames(ipw_s)==rownames(outcome_s))
all(rownames(ipw_n)==rownames(outcome_n))
all(rownames(ipw)==rownames(outcome_all))

#Resampling function
BalancedResampling_gaus <- function(data, tau, Z, ...) {
  s <- NULL
  formule <- list()
  for (i in 1:(length(names(Z)))){
    formule[[i]] <- paste("c(s, sample(which(",paste(paste0("Z$",names(Z)[i]),"== 1)"),", size = tau * sum(",paste(paste0("Z$",names(Z)[i]),"== 1)"),"))",sep="")
  }
  for (i in 1:(length(names(Z)))) {
    s <- eval(parse(text=formule[[i]]))
  }
  return(s)
}

#factors to numeric
covars_all[] <- lapply(covars_all, function(col) {
  if (is.factor(col)) as.numeric(as.character(col)) else col})
covars_s[] <- lapply(covars_s, function(col) {
  if (is.factor(col)) as.numeric(as.character(col)) else col})
covars_n[] <- lapply(covars_n, function(col) {
  if (is.factor(col)) as.numeric(as.character(col)) else col})

#EWAS bibliography filter-------------------------------------------------------
imp_data_fil <- imp_data[colnames(imp_data) %in% bibliography_cpgs]
imp_data_n_fil <- imp_data_n[colnames(imp_data_n) %in% bibliography_cpgs]
imp_data_s_fil <- imp_data_s[colnames(imp_data_s) %in% bibliography_cpgs]

#EWAS_North
penalty_n <- c(rep(1, ncol(imp_data_n_fil)),rep(0, ncol(covars_n))) 

n_dia1 <- VariableSelection(xdata = cbind(imp_data_n_fil, covars_n), 
                            ydata = outcome_n[1],
                            family = "gaussian", 
                            penalty.factor = penalty_n,
                            Lambda = LambdaSequence(lmax = 20, lmin = 1e-10, cardinal = 50),
                            seed = c(12345), 
                            weights = ipw_n$ipw,  
                            resampling = BalancedResampling_gaus, 
                            Z = covars_n[c('cohort_BIB', "cohort_EDEN", "cohort_KANC")], 
                            tau=0.8, 
                            K = 100, 
                            n_cat=3 , 
                            pi_list=seq(0.51, 0.99, by = 0.01), 
                            beep=NULL, 
                            verbose=TRUE, 
                            n_cores=16)

n_sys1 <- VariableSelection(xdata = cbind(imp_data_n_fil, covars_n), 
                            ydata = outcome_n[2],
                            family = "gaussian", 
                            penalty.factor = penalty_n,
                            Lambda = LambdaSequence(lmax = 20, lmin = 1e-10, cardinal = 50),
                            seed = c(12345), 
                            weights = ipw_n$ipw,  
                            resampling = BalancedResampling_gaus, 
                            Z = covars_n[c('cohort_BIB', "cohort_EDEN", "cohort_KANC")], 
                            tau=0.8, 
                            K = 100, 
                            n_cat=3 , 
                            pi_list=seq(0.51, 0.99, by = 0.01), 
                            beep=NULL, 
                            verbose=TRUE, 
                            n_cores=16)

n_dia2 <- VariableSelection(xdata = cbind(imp_data_n_fil, covars_n,outcome_n[2]), 
                            ydata = outcome_n[3],
                            family = "gaussian", 
                            penalty.factor = c(penalty_n,0),
                            Lambda = LambdaSequence(lmax = 20, lmin = 1e-10, cardinal = 50),
                            seed = c(12345), 
                            weights = ipw_n$ipw,  
                            resampling = BalancedResampling_gaus, 
                            Z = covars_n[c('cohort_BIB', "cohort_EDEN", "cohort_KANC")], 
                            tau=0.8, 
                            K = 100, 
                            n_cat=3 , 
                            pi_list=seq(0.51, 0.99, by = 0.01), 
                            beep=NULL, 
                            verbose=TRUE, 
                            n_cores=16)

n_sys2 <- VariableSelection(xdata = cbind(imp_data_n_fil, covars_n,outcome_n[3]), 
                            ydata = outcome_n[4],
                            family = "gaussian", 
                            penalty.factor = c(penalty_n,0),
                            Lambda = LambdaSequence(lmax = 20, lmin = 1e-10, cardinal = 50),
                            seed = c(12345), 
                            weights = ipw_n$ipw,  
                            resampling = BalancedResampling_gaus, 
                            Z = covars_n[c('cohort_BIB', "cohort_EDEN", "cohort_KANC")], 
                            tau=0.8, 
                            K = 100, 
                            n_cat=3 , 
                            pi_list=seq(0.51, 0.99, by = 0.01), 
                            beep=NULL, 
                            verbose=TRUE, 
                            n_cores=16)


res_list_n <- list("Results zdia t1"= n_dia1,
                   "Results zsys t1"= n_sys1,
                   "Results zdia t2"= n_dia2,
                   "Results zsys t2"= n_sys2)

write_rds(res_list_n,"./results/EWAS/results_n_ewascat.RDS")

#EWAS_South
penalty_s <- c(rep(1, ncol(imp_data_s_fil)),rep(0, ncol(covars_s))) 

s_dia1 <- VariableSelection(xdata = cbind(imp_data_s_fil, covars_s), 
                            ydata = outcome_s[1],
                            family = "gaussian", 
                            penalty.factor = penalty_s,
                            Lambda = LambdaSequence(lmax = 20, lmin = 1e-10, cardinal = 50),
                            seed = c(12345), 
                            weights = ipw_s$ipw,  
                            resampling = BalancedResampling_gaus, 
                            Z = covars_s[c('cohort_SAB')], 
                            tau=0.8, 
                            K = 100, 
                            n_cat=3 , 
                            pi_list=seq(0.51, 0.99, by = 0.01), 
                            beep=NULL, 
                            verbose=TRUE, 
                            n_cores=16)

s_sys1 <- VariableSelection(xdata = cbind(imp_data_s_fil, covars_s), 
                            ydata = outcome_s[2],
                            family = "gaussian", 
                            penalty.factor = penalty_s,
                            Lambda = LambdaSequence(lmax = 20, lmin = 1e-10, cardinal = 50),
                            seed = c(12345), 
                            weights = ipw_s$ipw,  
                            resampling = BalancedResampling_gaus, 
                            Z = covars_s[c('cohort_SAB')], 
                            tau=0.8, 
                            K = 100, 
                            n_cat=3 , 
                            pi_list=seq(0.51, 0.99, by = 0.01), 
                            beep=NULL, 
                            verbose=TRUE, 
                            n_cores=16)

s_dia2 <- VariableSelection(xdata = cbind(imp_data_s_fil, covars_s,outcome_s[2]), 
                            ydata = outcome_s[3],
                            family = "gaussian", 
                            penalty.factor = c(penalty_s,0),
                            Lambda = LambdaSequence(lmax = 20, lmin = 1e-10, cardinal = 50),
                            seed = c(12345), 
                            weights = ipw_s$ipw,  
                            resampling = BalancedResampling_gaus, 
                            Z = covars_s[c('cohort_SAB')], 
                            tau=0.8, 
                            K = 100, 
                            n_cat=3 , 
                            pi_list=seq(0.51, 0.99, by = 0.01), 
                            beep=NULL, 
                            verbose=TRUE, 
                            n_cores=16)

s_sys2 <- VariableSelection(xdata = cbind(imp_data_s_fil, covars_s,outcome_s[3]), 
                            ydata = outcome_s[4],
                            family = "gaussian", 
                            penalty.factor = c(penalty_s,0),
                            Lambda = LambdaSequence(lmax = 20, lmin = 1e-10, cardinal = 50),
                            seed = c(12345), 
                            weights = ipw_s$ipw,  
                            resampling = BalancedResampling_gaus, 
                            Z = covars_s[c('cohort_SAB')], 
                            tau=0.8, 
                            K = 100, 
                            n_cat=3 , 
                            pi_list=seq(0.51, 0.99, by = 0.01), 
                            beep=NULL, 
                            verbose=TRUE, 
                            n_cores=16)


res_list_s <- list("Results zdia t1"= s_dia1,
                   "Results zsys t1"= s_sys1,
                   "Results zdia t2"= s_dia2,
                   "Results zsys t2"= s_sys2)

write_rds(res_list_s,"./results/EWAS/results_s_ewascat.RDS")

#EWAS_all
penalty <- c(rep(1, ncol(imp_data_fil)),rep(0, ncol(covars_all))) 

all_dia1 <- VariableSelection(xdata = cbind(imp_data_fil, covars_all), 
                  ydata = outcome_all[1],
                  family = "gaussian", 
                  penalty.factor = penalty,
                  Lambda = LambdaSequence(lmax = 20, lmin = 1e-10, cardinal = 50),
                  seed = c(12345), 
                  weights = ipw$ipw,  
                  resampling = BalancedResampling_gaus, 
                  Z = covars_all[c('cohort_BIB', "cohort_EDEN", "cohort_KANC", "cohort_MOBA","cohort_SAB")], 
                  tau=0.8, 
                  K = 100, 
                  n_cat=3 , 
                  pi_list=seq(0.51, 0.99, by = 0.01), 
                  beep=NULL, 
                  verbose=TRUE, 
                  n_cores=16)

all_sys1 <- VariableSelection(xdata = cbind(imp_data_fil, covars_all), 
                  ydata = outcome_all[2],
                  family = "gaussian", 
                  penalty.factor = penalty,
                  Lambda = LambdaSequence(lmax = 20, lmin = 1e-10, cardinal = 50),
                  seed = c(12345), 
                  weights = ipw$ipw,  
                  resampling = BalancedResampling_gaus, 
                  Z = covars_all[c('cohort_BIB', "cohort_EDEN", "cohort_KANC", "cohort_MOBA","cohort_SAB")], 
                  tau=0.8, 
                  K = 100, 
                  n_cat=3 , 
                  pi_list=seq(0.51, 0.99, by = 0.01), 
                  beep=NULL, 
                  verbose=TRUE, 
                  n_cores=16)


all_dia2 <- VariableSelection(xdata = cbind(imp_data_fil, covars_all,outcome_all[2]), 
                  ydata = outcome_all[3],
                  family = "gaussian", 
                  penalty.factor = c(penalty,0),
                  Lambda = LambdaSequence(lmax = 20, lmin = 1e-10, cardinal = 50),
                  seed = c(12345), 
                  weights = ipw$ipw,  
                  resampling = BalancedResampling_gaus, 
                  Z = covars_all[c('cohort_BIB', "cohort_EDEN", "cohort_KANC", "cohort_MOBA","cohort_SAB")], 
                  tau=0.8, 
                  K = 100, 
                  n_cat=3 , 
                  pi_list=seq(0.51, 0.99, by = 0.01), 
                  beep=NULL, 
                  verbose=TRUE, 
                  n_cores=16)

all_sys2 <- VariableSelection(xdata = cbind(imp_data_fil, covars_all,outcome_all[3]), 
                  ydata = outcome_all[4],
                  family = "gaussian", 
                  penalty.factor = c(penalty,0),
                  Lambda = LambdaSequence(lmax = 20, lmin = 1e-10, cardinal = 50),
                  seed = c(12345), 
                  weights = ipw$ipw,  
                  resampling = BalancedResampling_gaus, 
                  Z = covars_all[c('cohort_BIB', "cohort_EDEN", "cohort_KANC", "cohort_MOBA","cohort_SAB")], 
                  tau=0.8, 
                  K = 100, 
                  n_cat=3 , 
                  pi_list=seq(0.51, 0.99, by = 0.01), 
                  beep=NULL, 
                  verbose=TRUE, 
                  n_cores=16)


res_list_all <- list("Results zdia t1"= all_dia1,
                     "Results zsys t1"= all_sys1,
                     "Results zdia t2"= all_dia2,
                     "Results zsys t2"= all_sys2)

write_rds(res_list_all,"./results/EWAS/results_all_ewascat.RDS")



# Results N
source("./codi/Functions/function_volcano_plot.R")

for (m in 1:length(res_list_n)){
  name<- names(res_list_n[m])
  stab_m <- res_list_n[[m]]
  class(stab_m) <- "variable_selection"
  
  selected_m <- SelectedVariables(stab_m)
  selected_names <- names(selected_m[selected_m==1])
  selprop_m <- SelectionProportions(stab_m)
  selected_selprop <- selprop_m[selected_m==1]
  pi_list=seq(0.51, 0.99, by = 0.01)
  argmax_id <- ArgmaxId(stab_m)[2]
  
  tmp <- t(stab_m$Beta[argmax_id, colnames(stab_m$selprop), ])
  beta_m <- apply(tmp, 2, FUN = function(x) {mean(x[x != 0])})
  selected_beta <- beta_m[selected_m==1]
  
  write.csv(
    cbind("Names"=selected_names,
          "1_Selprop"=1-selected_selprop,
          "Beta"=selected_beta), row.names = FALSE,
    file = paste0("./results/EWAS/Enrichment/north_",name, ".csv")
  )
  
  cat(name, ":",  kableExtra::kable(selected_selprop), sep ="\n" )
  volcano_plot(filename = paste0("N_",name,".svg"))  
  
  svg(paste0("./results/EWAS/Calibration_Plots/Cal_Plot_N_", name, ".svg"))
  CalibrationPlot(stab_m)
  dev.off()

}

# Results S

for (m in 1:length(res_list_s)){
  name<- names(res_list_s[m])
  stab_m <- res_list_s[[m]]
  class(stab_m) <- "variable_selection"
  
  selected_m <- SelectedVariables(stab_m)
  selected_names <- names(selected_m[selected_m==1])
  selprop_m <- SelectionProportions(stab_m)
  selected_selprop <- selprop_m[selected_m==1]
  pi_list=seq(0.51, 0.99, by = 0.01)
  argmax_id <- ArgmaxId(stab_m)[2]
  
  tmp <- t(stab_m$Beta[argmax_id, colnames(stab_m$selprop), ])
  beta_m <- apply(tmp, 2, FUN = function(x) {mean(x[x != 0])})
  selected_beta <- beta_m[selected_m==1]
  
  write.csv(
    cbind("Names"=selected_names,
          "1_Selprop"=1-selected_selprop,
          "Beta"=selected_beta), row.names = FALSE,
    file = paste0("./results/EWAS/Enrichment/south_",name, ".csv")
  )
  
  cat(name, ":",  kableExtra::kable(selected_selprop), sep ="\n" )
  volcano_plot(filename = paste0("S_",name,".svg"))
  
  svg(paste0("./results/EWAS/Calibration_Plots/Cal_Plot_S_", name, ".svg"))
  CalibrationPlot(stab_m)
  xlim
  dev.off()
  
}

# Results all
for (m in 1:length(res_list_all)){
  name<- names(res_list_all[m])
  stab_m <- res_list_all[[m]]
  class(stab_m) <- "variable_selection"
  
  selected_m <- SelectedVariables(stab_m)
  selected_names <- names(selected_m[selected_m==1])
  selprop_m <- SelectionProportions(stab_m)
  selected_selprop <- selprop_m[selected_m==1]
  pi_list=seq(0.51, 0.99, by = 0.01)
  argmax_id <- ArgmaxId(stab_m)[2]
  
  tmp <- t(stab_m$Beta[argmax_id, colnames(stab_m$selprop), ])
  beta_m <- apply(tmp, 2, FUN = function(x) {mean(x[x != 0])})
  selected_beta <- beta_m[selected_m==1]
  
  write.csv(
    cbind("Names"=selected_names,
          "1_Selprop"=1-selected_selprop,
          "Beta"=selected_beta), row.names = FALSE,
    file = paste0("./results/EWAS/Enrichment/all_",name, ".csv")
  )

  cat(name, ":",  kableExtra::kable(selected_selprop), sep ="\n" )
  volcano_plot(filename = paste0("all_",name,".svg"))
  
  svg(paste0("./results/EWAS/Calibration_Plots/Cal_Plot_all_", name, ".svg"))
  CalibrationPlot(stab_m)
  dev.off()
  
}

#Incremental plots
#North
inc_n_dia1 <-Incremental(xdata = cbind(imp_data_n_fil, covars_n), 
                  ydata = outcome_n[1],
                  family = "gaussian",
                  stability = n_dia1,
                  seed = c(12345), 
                  resampling = BalancedResampling_gaus,
                  Z = covars_n[c('cohort_BIB', "cohort_EDEN", "cohort_KANC")], 
                  tau=0.8, 
                  K = 100)

IncrementalPlot(inc_n_dia1)

inc_n_sys1 <-Incremental(xdata = cbind(imp_data_n_fil, covars_n), 
                         ydata = outcome_n[2],
                         family = "gaussian",
                         stability = n_sys1,
                         seed = c(12345), 
                         resampling = BalancedResampling_gaus,
                         Z = covars_n[c('cohort_BIB', "cohort_EDEN", "cohort_KANC")], 
                         tau=0.8, 
                         K = 100)
IncrementalPlot(inc_n_sys1)


inc_n_dia2 <-Incremental(xdata = cbind(imp_data_n_fil, covars_n, outcome_n[1]), 
                         ydata = outcome_n[3],
                         family = "gaussian",
                         stability = n_dia2,
                         seed = c(12345), 
                         resampling = BalancedResampling_gaus,
                         Z = covars_n[c('cohort_BIB', "cohort_EDEN", "cohort_KANC")], 
                         tau=0.8, 
                         K = 100)
IncrementalPlot(inc_n_dia2)

inc_n_sys2 <-Incremental(xdata = cbind(imp_data_n_fil, covars_n, outcome_n[2]), 
                         ydata = outcome_n[4],
                         family = "gaussian",
                         stability = n_sys2,
                         seed = c(12345), 
                         resampling = BalancedResampling_gaus,
                         Z = covars_n[c('cohort_BIB', "cohort_EDEN", "cohort_KANC")], 
                         tau=0.8, 
                         K = 100)
IncrementalPlot(inc_n_sys2)

#South
inc_s_dia1 <-Incremental(xdata = cbind(imp_data_s_fil, covars_s), 
                         ydata = outcome_s[1],
                         family = "gaussian",
                         stability = s_dia1,
                         seed = c(12345), 
                         resampling = BalancedResampling_gaus,
                         Z = covars_s[c('cohort_SAB')], 
                         tau=0.8, 
                         K = 100)

IncrementalPlot(inc_s_dia1)

inc_s_sys1 <-Incremental(xdata = cbind(imp_data_s_fil, covars_s), 
                         ydata = outcome_s[2],
                         family = "gaussian",
                         stability = s_sys1,
                         seed = c(12345), 
                         resampling = BalancedResampling_gaus,
                         Z = covars_s[c('cohort_SAB')], 
                         tau=0.8, 
                         K = 100)
IncrementalPlot(inc_s_sys1)


inc_s_dia2 <-Incremental(xdata = cbind(imp_data_s_fil, covars_s, outcome_s[1]), 
                         ydata = outcome_s[3],
                         family = "gaussian",
                         stability = s_dia2,
                         seed = c(12345), 
                         resampling = BalancedResampling_gaus,
                         Z = covars_s[c('cohort_SAB')], 
                         tau=0.8, 
                         K = 100)
IncrementalPlot(inc_s_dia2)

inc_s_sys2 <-Incremental(xdata = cbind(imp_data_s_fil, covars_s, outcome_s[2]), 
                         ydata = outcome_s[4],
                         family = "gaussian",
                         stability = s_sys2,
                         seed = c(12345), 
                         resampling = BalancedResampling_gaus,
                         Z = covars_s[c('cohort_SAB')], 
                         tau=0.8, 
                         K = 100)
IncrementalPlot(inc_s_sys2)

#All

inc_all_dia1 <-Incremental(xdata = cbind(imp_data_fil, covars_all), 
                         ydata = outcome_all[1],
                         family = "gaussian",
                         stability = all_dia1,
                         seed = c(12345), 
                         resampling = BalancedResampling_gaus,
                         Z = covars_all[c('cohort_BIB', "cohort_EDEN", "cohort_KANC","cohort_SAB")], 
                         tau=0.8, 
                         K = 100)

IncrementalPlot(inc_all_dia1)

inc_all_sys1 <-Incremental(xdata = cbind(imp_data_fil, covars_all), 
                         ydata = outcome_all[2],
                         family = "gaussian",
                         stability = all_sys1,
                         seed = c(12345), 
                         resampling = BalancedResampling_gaus,
                         Z = covars_all[c('cohort_BIB', "cohort_EDEN", "cohort_KANC")], 
                         tau=0.8, 
                         K = 100)
IncrementalPlot(inc_all_sys1)


inc_all_dia2 <-Incremental(xdata = cbind(imp_data_fil, covars_all, outcome_all[1]), 
                         ydata = outcome_all[3],
                         family = "gaussian",
                         stability = all_dia2,
                         seed = c(12345), 
                         resampling = BalancedResampling_gaus,
                         Z = covars_all[c('cohort_BIB', "cohort_EDEN", "cohort_KANC")], 
                         tau=0.8, 
                         K = 100)
IncrementalPlot(inc_all_dia2)

inc_all_sys2 <-Incremental(xdata = cbind(imp_data_fil, covars_all, outcome_all[2]), 
                         ydata = outcome_all[4],
                         family = "gaussian",
                         stability = all_sys2,
                         seed = c(12345), 
                         resampling = BalancedResampling_gaus,
                         Z = covars_all[c('cohort_BIB', "cohort_EDEN", "cohort_KANC")], 
                         tau=0.8, 
                         K = 100)
IncrementalPlot(inc_all_sys2)

#Save plots
ggsave(filename = "./results/EWAS/Incremental_plots/inc_n_dia1.png", plot = IncrementalPlot(inc_n_dia1), width = 8, height = 6, dpi = 300)
ggsave(filename = "./results/EWAS/Incremental_plots/inc_n_sys1.png", plot = IncrementalPlot(inc_n_sys1), width = 8, height = 6, dpi = 300)
ggsave(filename = "./results/EWAS/Incremental_plots/inc_n_dia2.png", plot = IncrementalPlot(inc_n_dia2), width = 8, height = 6, dpi = 300)
ggsave(filename = "./results/EWAS/Incremental_plots/inc_n_sys2.png", plot = IncrementalPlot(inc_n_sys2), width = 8, height = 6, dpi = 300)

ggsave(filename = "./results/EWAS/Incremental_plots/inc_s_dia1.png", plot = IncrementalPlot(inc_s_dia1), width = 8, height = 6, dpi = 300)
ggsave(filename = "./results/EWAS/Incremental_plots/inc_s_sys1.png", plot = IncrementalPlot(inc_s_sys1), width = 8, height = 6, dpi = 300)
ggsave(filename = "./results/EWAS/Incremental_plots/inc_s_dia2.png", plot = IncrementalPlot(inc_s_dia2), width = 8, height = 6, dpi = 300)
ggsave(filename = "./results/EWAS/Incremental_plots/inc_s_sys2.png", plot = IncrementalPlot(inc_s_sys2), width = 8, height = 6, dpi = 300)

ggsave(filename = "./results/EWAS/Incremental_plots/inc_all_dia1.png", plot = IncrementalPlot(inc_all_dia1), width = 8, height = 6, dpi = 300)
ggsave(filename = "./results/EWAS/Incremental_plots/inc_all_sys1.png", plot = IncrementalPlot(inc_all_sys1), width = 8, height = 6, dpi = 300)
ggsave(filename = "./results/EWAS/Incremental_plots/inc_all_dia2.png", plot = IncrementalPlot(inc_all_dia2), width = 8, height = 6, dpi = 300)
ggsave(filename = "./results/EWAS/Incremental_plots/inc_all_sys2.png", plot = IncrementalPlot(inc_all_sys2), width = 8, height = 6, dpi = 300)

#filter cpgs spls_N/S-----------------------------------------------------------


#filter cpgs spls_all-----------------------------------------------------------






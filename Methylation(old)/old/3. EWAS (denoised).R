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
library(wateRmelon)
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
# Performances in test set of model refitted in training set
#require(sharp)
# detach("package:sharp",unload= TRUE) # Detach the sharp function (CRAN)
# #Clean version
# # install.packages("./package_versions/clean/sharp-main",
# #                  repos= NULL,
# #                  type = "source", force=TRUE) # Load Sharp function modified to apply weights
# #Modified version
# install.packages("./package_versions/modified_offset/sharp-main",
#                  repos= NULL,
#                  type = "source", force=TRUE) # Load Sharp function modified to apply weights
# # install.packages("sharp",force=TRUE)
library(sharp)

# //////// LOADING OMICS DATASET

# /// Load HelixID by group (N/S)
helixid_n <- c(read.csv("./helixid_n.csv", row.names = 1)$x)
helixid_s <- c(read.csv("./helixid_s.csv", row.names = 1)$x)

#Loading cpgs list
ewas_cat_df<-read_tsv("ewas_catalog_cpgs.tsv")
sPLS_n<-read_xlsx("./results/sPLS_results/Methyl/methyl_sPLS_n.xlsx")
sPLS_s<-read_xlsx("./results/sPLS_results/Methyl/methyl_sPLS_s.xlsx")
sPLS<-read_xlsx("./results/sPLS_results/Methyl/methyl_sPLS_all.xlsx")

#Selecting cpgs associated with BP
ewas_cat_df$trait<-tolower(ewas_cat_df$trait)
ewas_cat_df<-ewas_cat_df[ewas_cat_df$trait %in% c("systolic blood pressure", "diastolic blood pressure"),]
ewas_cat_cpgs <- unique(ewas_cat_df$cpg)

spls_n_cpgs <-sPLS_n$results[which(sPLS_n$results$SPLS_MIN.TP=="1"),]
spls_s_cpgs <-sPLS_s$results[which(sPLS_s$results$SPLS_MIN.TP=="1"),]
spls_all_cpgs <-sPLS$results[which(sPLS$results$SPLS_MIN.TP=="1"),]



# Loading methylation data:

omics_N<- readRDS("./results/denoising/denoising_methyl/methyl_denoised_n.RDS")
omics_S <- readRDS("./results/denoising/denoising_methyl/methyl_denoised_s.RDS")
omics_all <- readRDS("./results/denoising/denoising_methyl/methyl_denoised.RDS")

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

#filter cpgs EWAS catalog-------------------------------------------------------
imp_data_fil <- imp_data[colnames(imp_data) %in% ewas_cat_cpgs]
imp_data_n_fil <- imp_data_n[colnames(imp_data_n) %in% ewas_cat_cpgs]
imp_data_s_fil <- imp_data_s[colnames(imp_data_s) %in% ewas_cat_cpgs]

pred_all <- cbind(imp_data, )

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


#EWAS_all
penalty <- c(rep(1, ncol(imp_data_fil)),rep(0, ncol(covars_all))) 
names(outcome_all)
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
#filter cpgs spls_N/S-----------------------------------------------------------


#filter cpgs spls_all-----------------------------------------------------------






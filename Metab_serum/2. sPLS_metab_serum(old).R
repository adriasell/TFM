################################# Options

options(max.print = 1000000)
options(dplyr.print_max = 1000000)
rm(list=ls())
date <- Sys.Date()

################################# Set R environment 

list.of.packages <- c(
  "haven","foreach","Matrix","MASS","parallel","glmnet","spls","MXM","dsa","dlnm","splines","mgcv","doParallel",
  "ranger","palmerpenguins","tidyverse","kableExtra","haven", "BiocManager", "Biobase", "sharp", "writexl"
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}
for(package.i in list.of.packages){ #loading packages
  suppressPackageStartupMessages(
    library(
      package.i, 
      character.only = TRUE
    )
  )
}

################################# Defining working directory.

setwd("./TFM")
getwd()

################################# Loading data.

# /// Load HelixID by group (N/S)

helixid_n <- c(read.csv("./helixid_n.csv", row.names = 1)$x)
helixid_s <- c(read.csv("./helixid_s.csv", row.names = 1)$x)

# /// Load omic data (X)
omics_denoised <- readRDS("./results/denoising/denoising_metab/serum/metab_serum_denoised.RDS")
omics_n_denoised <- omics_denoised[,which(pData(omics_denoised)$HelixID.x %in% helixid_n)]
omics_s_denoised <- omics_denoised[,which(pData(omics_denoised)$HelixID.x %in% helixid_s)]


imp_data_n <- t(exprs(omics_n_denoised))
imp_data_s <- t(exprs(omics_s_denoised))
imp_data <- t(exprs(omics_denoised))


rownames(imp_data_n) <- pData(omics_n_denoised)$HelixID.x
rownames(imp_data_s) <- pData(omics_s_denoised)$HelixID.x
rownames(imp_data) <- pData(omics_denoised)$HelixID.x


# /// Load outcomes (Y)

phenotype <- readRDS("./db/pheno/final/bp_wide_validN5332023-10-16.rds")
outcomes <- c("hs2_zdia_bp.v3_2017_Time1", "hs2_zsys_bp.v3_2017_Time1",
              "hs2_zdia_bp.v3_2017_Time2", "hs2_zsys_bp.v3_2017_Time2")
rownames(phenotype) <- phenotype$HelixID

outcome_n <- phenotype[phenotype$HelixID %in% helixid_n, outcomes]
outcome_n <- outcome_n[rownames(outcome_n) %in% rownames(imp_data_n),]
outcome_s <- phenotype[phenotype$HelixID %in% helixid_s, outcomes]
outcome_all <- phenotype[, outcomes]


#Delete NAs
any(is.na(outcome_s))
any(is.na(outcome_n)) 
any(is.na(outcome_all))

nas_id_n <- rownames(outcome_n[!complete.cases(outcome_n),])
outcome_n <- outcome_n[!rownames(outcome_n) %in% nas_id_n,]
imp_data_n <- imp_data_n[!rownames(imp_data_n) %in% nas_id_n,]

nas_id <- rownames(outcome_all[!complete.cases(outcome_all),])
imp_data <- imp_data[!rownames(imp_data) %in% nas_id,]
outcome_all <- outcome_all[rownames(outcome_all) %in% rownames(imp_data),]

outcome_s<- outcome_s[rownames(outcome_s) %in% rownames(imp_data_s),]

#Comprovacions

sort_by_rownames <- function(df) {
  df[order(rownames(df)),]
}

outcome_all<-sort_by_rownames(outcome_all)
outcome_n<-sort_by_rownames(outcome_n)
outcome_s<-sort_by_rownames(outcome_s)
imp_data<-sort_by_rownames(imp_data)
imp_data_n<-sort_by_rownames(imp_data_n)
imp_data_s<-sort_by_rownames(imp_data_s)


all(rownames(outcome_all)==rownames(imp_data))
all(rownames(outcome_n)==rownames(imp_data_n))
all(rownames(outcome_s)==rownames(imp_data_s))

################################ RUN MODELS

# # # # # # # # #
# # #  sPLS # # #
# # # # # # # # #
set.seed(1899)

source(file="./codi/Functions/sPLS_function_mod.R")

#*##*#*##*#*##*#*##*#*##*#*##*#*##*#*##*#*##*#*##*#*##
# PROBLEM: DIFFERENT SUBSETS OF VARIABLES ARE SELECTED EACH TIME THE MODEL IS RUN.
#*##*#*##*#*##*#*##*#*##*#*##*#*##*#*##*#*##*#*##*#*##


sPLS_n = applySPLS(data.Y=as.matrix(outcome_n),data.X=as.matrix(scale(imp_data_n)))
sPLS_s = applySPLS(data.Y=as.matrix(outcome_s),data.X=as.matrix(scale(imp_data_s)))
sPLS = applySPLS(data.Y=as.matrix(outcome_all),data.X=as.matrix(scale(imp_data)))

sPLS_n$results[which(sPLS_n$results$SPLS_MIN.TP=="1"),]
sPLS_s$results[which(sPLS_s$results$SPLS_MIN.TP=="1"),]
sPLS$results[which(sPLS$results$SPLS_MIN.TP=="1"),]

write_xlsx(sPLS_n, "./results/sPLS_results/serum/serum_sPLS_n.xlsx")
write_xlsx(sPLS_s, "./results/sPLS_results/serum/serum_sPLS_s.xlsx")
write_xlsx(sPLS, "./results/sPLS_results/serum/serum_sPLS_all.xlsx")

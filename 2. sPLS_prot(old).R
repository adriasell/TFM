################################# Options

options(max.print = 1000000)
options(dplyr.print_max = 1000000)
rm(list=ls())
date <- Sys.Date()

################################# Set R environment 

list.of.packages <- c(
  "haven","foreach","Matrix","MASS","parallel","glmnet","spls","MXM","dsa","dlnm","splines","mgcv","doParallel",
  "ranger","palmerpenguins","tidyverse","kableExtra","haven", "BiocManager", "Biobase", "sharp", "mixOmics", "writexl"
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
omics_denoised <- readRDS("./results/denoising/denoising_prot/prot_denoised.RDS")
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
outcomes <- c("hs2_zdia_bp.v3_2017_Time1", "hs2_zsys_bp.v3_2017_Time1", "hs2_zdia_bp.v3_2017_Time2", "hs2_zsys_bp.v3_2017_Time2")
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

source(file="./codi/Functions/sPLS_function_mod.R")

#*##*#*##*#*##*#*##*#*##*#*##*#*##*#*##*#*##*#*##*#*##
# PROBLEM: DIFFERENT SUBSETS OF VARIABLES ARE SELECTED EACH TIME THE MODEL IS RUN.
#*##*#*##*#*##*#*##*#*##*#*##*#*##*#*##*#*##*#*##*#*##

set.seed(1899)

sPLS_n = applySPLS(data.Y=as.matrix(outcome_n),data.X=as.matrix(scale(imp_data_n)))
sPLS_s = applySPLS(data.Y=as.matrix(outcome_s),data.X=as.matrix(scale(imp_data_s)))
sPLS = applySPLS(data.Y=as.matrix(outcome_all),data.X=as.matrix(scale(imp_data)))

sPLS_n$results[which(sPLS_n$results$SPLS_MIN.TP=="1"),]
sPLS_s$results[which(sPLS_s$results$SPLS_MIN.TP=="1"),]
sPLS$results[which(sPLS$results$SPLS_MIN.TP=="1"),]

write_xlsx(sPLS_n, "./results/sPLS_results/prot/prot_sPLS_n.xlsx")
write_xlsx(sPLS_s, "./results/sPLS_results/prot/prot_sPLS_s.xlsx")
write_xlsx(sPLS, "./results/sPLS_results/prot/prot_sPLS_all.xlsx")












#-------------------------------------------------------------------------------
# # # # # # # # # # # # #
# # SparsePLS (sharp) # #
# # # # # # # # # # # # #


# // Load weights 

ipw <- read.csv2("./db/ipws_bp_1feb24.csv")
rownames(ipw) <-ipw$HelixID
ipw_n <-ipw[ipw$HelixID %in% rownames(outcome_n),"ipw"]
ipw_s <-ipw[ipw$HelixID %in% rownames(outcome_s),"ipw"]

#sPLS North
sPLS_n<-VariableSelection(xdata = scale(imp_data_n), 
                          ydata =  outcome_n,
                          Lambda = LambdaSequence(lmax = 35, lmin = 0.1, cardinal = 100),
                          pi_list=seq(0.51, 0.99, by = 0.01), 
                          family = "gaussian",
                          implementation = SparsePLS,
                          weights = ipw_n, 
                          seed = c(12345), scale = F,
                          resampling = "subsampling", 
                          tau=0.8, K = 100,  n_cat=3, 
                          beep=NULL, verbose=T)

CalibrationPlot(sPLS_n)
Stable(sPLS_n)
argmax_id_n <- ArgmaxId(sPLS_n)[1]
beta_spls_n <- apply(t(sPLS_n$Beta[argmax_id_n,colnames(sPLS_n$Beta),]), 2, FUN = function(x) {
  mean(x[x != 0])})

#sPLS south
sPLS_s<-VariableSelection(xdata = scale(imp_data_s), 
                          ydata =  outcome_s,
                          Lambda = LambdaSequence(lmax = 35, lmin = 0.1, cardinal = 100),
                          pi_list=seq(0.51, 0.99, by = 0.01), 
                          family = "gaussian",
                          implementation = SparsePLS,
                          weights = ipw_s, 
                          seed = c(12345), scale = F,
                          resampling = "subsampling", 
                          tau=0.8, K = 100, n_cat=3, 
                          beep=NULL, verbose=TRUE)


CalibrationPlot(sPLS_s)
Stable(sPLS_s)
argmax_id_s <- ArgmaxId(sPLS_s)[1]
beta_spls_s <- apply(t(sPLS_s$Beta[argmax_id_s,colnames(sPLS_s$Beta),]), 2, FUN = function(x) {
  mean(x[x != 0])})


#sPLS all
sPLS_all<-VariableSelection(xdata = scale(imp_data), 
                            ydata =  outcome_all,
                            Lambda = LambdaSequence(lmax = 35, lmin = 0.1, cardinal = 100),
                            pi_list=seq(0.51, 0.99, by = 0.01), 
                            family = "gaussian",
                            implementation = SparsePLS,
                            weights = ipw$ipw, 
                            seed = c(12345), scale = F,
                            resampling = "subsampling", 
                            tau=0.8, K = 100, n_cat=3, 
                            beep=NULL, verbose=TRUE)

dim(imp_data)  # For xdata before scaling
dim(outcome)  # For ydata

CalibrationPlot(sPLS_all)
Stable(sPLS_all)
argmax_id_all <- ArgmaxId(sPLS_all)[1]
beta_spls_all <- apply(t(sPLS_all$Beta[argmax_id,colnames(sPLS_all$Beta),]), 2, FUN = function(x) {
  mean(x[x != 0])})



#FUNCIONS______________________________________________________________________

##### Extract beta
getbeta<- function(y){
  argmax_id <- ArgmaxId(y)[1]
  
  beta <- apply(t(y$Beta[argmax_id,colnames(y$Beta),]), 2, 
                FUN = function(x) {
                  mean(x[x != 0])})
  
  df <- data.frame(beta[1:36], Stable(y)) #arreglar 1:36
  colnames(df) <- c("Beta", "Selected")
  
  df1 <- data.frame(beta[37:40], 0) #arreglar 37:40
  colnames(df1) <- c("Beta", "Selected")
  
  df2 <- rbind(df, df1)
  
  return(df2)
}



# We write a function for the resampling strategy (keep the same proportion of case/controls and cohorts as in the original data, in all partitions)
BalancedResampling_bin <- function(data, tau, Z, ...) {
  s <- NULL
  formule <- list()
  for (i in 1:(length(names(Z)))){
    formule[[i]] <- c(paste("c(s, sample(which((data == 0) ",paste("& (",paste0("Z$",names(Z)[i]),"== 1)"),")",", size = tau * sum((data == 0) ",paste("& (",paste0("Z$",names(Z)[i]),"== 1)"),")))",sep=""),paste("c(s, sample(which((data == 1) ",paste("& (",paste0("Z$",names(Z)[i]),"== 1)"),")",", size = tau * sum((data == 1) ",paste("& (",paste0("Z$",names(Z)[i]),"== 1)"),")))",sep=""))
  }
  for (i in 1:(length(names(Z)))) {
    for (z in 1:2) {
      s <- eval(parse(text=formule[[i]][z]))
    }
  }
  return(s)
}

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

#-------------------------------------------------------------------------------


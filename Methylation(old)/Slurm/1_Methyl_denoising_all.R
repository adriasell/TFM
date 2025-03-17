################################################################
#################### Methylation  denoising ####################
################################################################

library(lumi)
library(doParallel)
library(minfi)
library(data.table)
library(readxl)
library(xlsx)
library(omics)
library(ggplot2)
library(corrplot)
library(PCAtools)
library(fastDummies)
library(polycor)
library(dplyr)
################################# Defining working directory.
getwd()
setwd("./TFM_25")
output_path <- "./results/denoising/denoising_methyl"

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

source("./script/denoising_2024/functions/lm_local.R")
source("./script/denoising_2024/functions/mlmer_local.R")
source("./script/denoising_2024/functions/generic_functions_denoising_v2.R")

# We add Prot.plate dummies afterwards to denoise_vars
load("/home/isglobal.lan/aseto/methyl_data.RData")
omics <- methyl_data$omics
denoise_vars <- methyl_data$denoise_vars
print("Data loaded")

print("Starting denoising")

omics_denoised <- denoising_mlmer(omics$data, omicsLayer, denoise_vars, ncores = detectCores(), blocksize = 100)

# ///// Save denoised data

saveRDS(omics_denoised, file = paste0(output_path, "/methyl_denoised.RDS"))
print("File saved")


#-------------------------------------------------------------
#               Check noised VS. denoised data
#-------------------------------------------------------------

# Correlation plot variables:
#   Define variables to be used in correlation plot. Please set which variables should be
#   used as numerical and which as a factor
#   
#       num_vars: Numerical variables or variables to be used as numerical
#       factor_vars: Factor variables or variables to be used as factor
#       
#           Example: 
#               num_vars <- c("e3_sex", "blood_sam4", "Prot.plate")
#               factor_vars <- c("cohort", "h_ethnicity_c")
#               factor_vars <- NA                                       # If we don't have factor variables
# Define the number of first principal components to be used in check, 
# If nComponents is NULL or NA, all components will be used 
# 
#       Example: 
#           nComponents <- 10 # Checks first 10 compoents
#           nComponents <- NA # Checks all the components
#           
nComponents <- 10
print("Starting corr_plots")


check_noised_denoised(original = omics$data, denoised = omics_denoised, 
                      Vars = denoise_vars, prefix = omicsLayer, topPCA = nComponents,
                      output_path=paste0(output_path,"/corr_plot"))
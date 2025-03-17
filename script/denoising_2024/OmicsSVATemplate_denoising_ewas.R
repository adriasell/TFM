# #######################################################################################################
#   TEMPLATE TO PERFORM SVA IN OMICS
# 
# 
# ---------------------
#       BIG OMICS
# ---------------------
# 
# Methylation:
#   /PROJECTES/HELIX_OMICS/data_final/methyl/child/8y/blood_450K_QChelix_20170401/methylome_subcohort_v4.Rdata            # (filtered)
#   /PROJECTES/HELIX_OMICS/data_final/methyl/child/8y/blood_450K_QChelix_20170401/methylome_subcohort_notfitr_v4.Rdata    # (not filtered)
# 
# Gene Expression:
#   /PROJECTES/HELIX_OMICS/data_final/trans/child/8y/blood_HTAv2_QChelix_20180701/transcriptome_subcohort_notfitr_v3.RData
#   
# miRNAs:
#   /PROJECTES/HELIX_OMICS/data_final/mirna/child/8y/blood_Agilent_QChelix_20180701/mirna_subcohort_notfiltr_inclsex_v4.RData
#   /PROJECTES/HELIX_OMICS/data_final/mirna/child/8y/blood_Agilent_QChelixR12_20200201/mirna_r12_subcohort_notfiltr_inclsex_v1.RData # ** DE MOMENT UTILITZAR AQUEST !!! ** #
# 
#  
# -----------------------
#       SMALL OMICS
# -----------------------
# 
#  Proteome:
#   /PROJECTES/HELIX_OMICS/data_final/prot/child/8y/plasma_Luminex_QChelix_20170301/proteome_subcohort_v5.Rdata
#  
#  Serum metabolome:
#   /PROJECTES/HELIX_OMICS/data_final/metab/child/8y/serum.urine_Biocrates.NMR_QChelix_20170101/metab_serum_subcohort_v3.RData
#  
#  Urine metabolome:
#   /PROJECTES/HELIX_OMICS/data_final/metab/child/8y/serum.urine_Biocrates.NMR_QChelix_20170101/metab_urine_subcohort_v3.RData
#   
# #######################################################################################################

if (!require(c("omics"), quietly = TRUE))
    install.packages("omics")
if (!require(c("keycovar"), quietly = TRUE))
    install.packages("keycovar")
if (!require(c("corrplot "), quietly = TRUE))
    install.packages("corrplot ")

library(minfi)
library(data.table)
library(readxl)
library(xlsx)
library(dplyr)
library(omics)

library(ggplot2)
library(corrplot )

# SVA packages
library(isva)
library(SmartSVA)

# PCA
library(PCAtools)

setwd("/PROJECTES/HELIX_OMICS/software/templates/Omics_Multivariate/")

# Load external functions
source("scripts/functions/extract_data_by_xlsx_styles.R")
source("scripts/functions/generic_functions_denoising_v2.R")
source("scripts/functions/sva_association.R")
source("scripts/functions/mlmer_local.R")


# ----------------------------------------------
#       Source data definition and loading
# ----------------------------------------------

# Omics Metadata
metadataFile <- "./db/metadata_aug/HELIX_SVA_common_OmicsMetadata_20230612_aug.RData" 
#"/home/aanguita@isglobal.lan/Escritorio/DOWNLOAD_MULTIOMICS_12_DIC/ATH_BPmultiomics_AA/script/denoising_2024/metadata/HELIX_SVA_common_OmicsMetadata_20231026.RData"


# Omic Data
# Define filepath and filename where omics data is stored, Rdata file with an ExpressionSet or GenomicRatioSet
#   Ex: 
#    omicFile <- "/PROJECTES/HELIX_OMICS/data_final/prot/child/8y/plasma_Luminex_QChelix_20170301/proteome_subcohort_v5.Rdata"

omicFile <- "./db/ewas/methylome_subcohort_notfitr_N533_sep23_nosex.Rds"


# Phenotype Data
#   Load phenotype data in phenotype variable
#       Examples:
#           - If phenotype data is stored in txt file with tabs
#               phenotype <- read.table(<filename>, sep="\t")
#           - If phenotype data is stored in escel format
#               phenotype <- read_excel( <filename>, na = c("NA", "."))
phenotype <- readRDS("./db/pheno/final/bp_wide_validN5332023-10-16.rds")
rownames(phenotype) <- phenotype$HelixID
ids
names(phenotype)
dim(phenotype)
phenotype$SampleID <- ids[phenotype$HelixID,"SampleId"]
phenotype$HelixID == ids[phenotype$HelixID,"HelixID"]
table(phenotype$FINAL_ancestry)
names(phenotype)


    
# ----------------------------------------------
#           Main variables definition
# ----------------------------------------------

workpath <- "/home/aanguita@isglobal.lan/Escritorio/DOWNLOAD_MULTIOMICS_12_DIC/ATH_BPmultiomics_AA"


# Set working path
setwd(workpath)

# Omic used
#   Define omics used in analysis omicLayer variable.
#   Possible values are: 
#       Methyl   : Methylation data
#       Gene     : Genomic data
#       Prot     : Proteomic data
#       miRNAs   : micro-RNAs
#       Serum    : Serum methabolome
#       Urine    : Urine metabolome
omicsLayer <- "EWAS"


# Set ID column name in the phenotype dataset to be used to relate 
# the phenotype to the Omic data. (by default SampleID)
phenotypeID <- "SampleID"


# Set the variables to be used from phenotype dataset (by default ALL)
# 
#  Examples:
#       phenoVariables <- c('bmi', 'sex', 'age', 'smoke')  # 
#       phenoVariables <- 'ALL'         # If we need ALL the variables present in phenotype dataset
#       phenoVariables <- ''            # If no value is defined then ALL the variables present in the phenotype will be used
#       

phenotypeVariables <- c("hs2_zdia_bp.v3_2017_Time1","hs2_zsys_bp.v3_2017_Time1","hs2_zbmi_who_Time1","hs2_height_c_zscore_Time1","e3_sex_Time1","hs2_visit_age_years_Time1")


# Ethnicity - Population
# 
#   ethnic parameter is used to make a subset of the data and also to take in to account
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
#       
#       ALL     : No filter
#       
#  Example:
#       ethnic <- 'ALL'                     # All populations - no filter
#       ethnic <- ''                        # All populations - no filter
#       ethnic <- c('AFR', 'AMR', 'SAS')    # Use African, Mixed American and South Asian populations
# 

ethnic <- c('ALL')



# ----------------------------------------------
#               SVA Variables
# ----------------------------------------------


#   SVA execution mode:
#       Method 1: With all variables - only protect the phenotype of interest
#       Method 2: Protect all covariables - protect the phenotype of interest and all the covariates - apply SVA only in technical variables
#       Method 3: Protect some covariables - protect the phenotype of interest and some covariates - apply SVA in technical variables and some covariates
#   
#   Possible values:
#       'Interest'  : Apply Method 1.
#       'Full'      : Apply Method 2.
#       'Partial'   : Apply Method 3.
#       
#        Example:
#           svaMode <- c('Interest')                    # Applly method 1
#           svaMode <- c('Interest', 'Partial')         # Applly method 1 and 3
#           svaMode <- c('Interest', 'Full', 'Partial') # Applly method 1, 2 and 3

svaMode <- c('Partial') #c( 'Interest', 'Partial', 'Full')


# Variable of interest
#   Set the variable of interest - used to protect in SVA
#       Example:
#           interest_var <- 'hs_zbmi_who'

interest_var <- 'hs2_zsys_bp.v3_2017_Time1'


# Covariates
#   Set the covariates to be used in the model
#   
#       Example:
#           covariates <- c('sex', 'bmi', 'age') # We want to use sex, bmi and age in our model as covariates
#           covariates <- c('e3_sex_None', 'h_cohort', 'h_ethnicity_c')

covariates <- c("NK","Bcell","CD4T","CD8T","Eos","Mono","Neu",'cohort', 'h_ethnicity_c',"age_sample_months","hs2_zbmi_who_Time1","hs2_height_c_zscore_Time1","e3_sex_Time1","hs2_visit_age_years_Time1")


# Covariates to be protected in SVA
# 
#   Define only if you are using the 'Partial' method and you want to protect these covariates
#   to be used to compute the Surrogate variables (SVA)
#   
#   NOTE: this setting is only used in case you choose to use the method 'Partial' (Method 2)
#   
#       Example:
#           to_protect_covariates <- c('sex')  # The sex variable will be protected so that it is not used to compute the surrogate variables
#           to_protect_covariates <- c('sex', 'age')  # The variables sex and age  will be protected, these variables will not be used to compute the surrogate variables
to_protect_covariates <- c("NK","Bcell","CD4T","CD8T","Eos","Mono","Neu","cohort", "h_ethnicity_c","age_sample_months","hs2_zbmi_who_Time1","hs2_height_c_zscore_Time1","e3_sex_Time1","hs2_visit_age_years_Time1","hs2_zdia_bp.v3_2017_Time1")

# Plot Variables
# Define variables to use to display in the plot of the surrogate variables SV1 vs SV2 and SV2 vs SV3
# 
#   Example:
#       plot_vars <-  c("h_cohort", "e3_sex_None")  # Plot SV1 vs SV2 and SV2 vs SV3 for cohort and sex
#       plot_vars <- NA                             # No plots

plot_vars <-  c("cohort","e3_sex_Time1")



# ----------------------------------------------
#               Prepare data
# ----------------------------------------------

setwd(workpath) 

if(phenotypeVariables == 'ALL' || 'ALL' %in% phenotypeVariables) {
    omicSVA <- getfullOmicsPhenotype(omicFile, metadataFile, phenotype, phenotypeID, phenotypeVariables, omicsLayer, ethnic )    
} else {
    omicSVA <- getfullOmicsPhenotype(omicFile, metadataFile, phenotype, phenotypeID, 
                                     c( phenotypeID, phenotypeVariables, interest_var, covariates, to_protect_covariates), 
                                     omicsLayer, ethnic )
}



#-------------------------------------------------------------
#                        DENOISING 
#-------------------------------------------------------------
# Denoise data using linear mixed effects model
#-------------------------------------------------------------


# Denoising variables:
#   Define variables to be used to denoise data. Suggested vars to be used for denoising 
#   each omic are: 
#   
#       * Metylation
#           - nano_conc_ng_ul
#           - extr_batch
#           - Sample_Plate
#           - BATCH_Conversion - IMPORTANT
#           - BATCH_Infinium
#           - Array - MOST IMPORTANT
#           - Slide - IMPORTANT
#           - SampleType
# 
#       * Gene expression:
#           - round_gexp
#           - r_b IMPORTANT
#           - round_rna_extr
#           - extr_batch
#           - bio_conc_denat_ngul IMPORTANT
#           - bio_rin_denat IMPORTANT
#           - extraction_successful_class
#           - sample_milky_blood_class
#           - SampleType
#           
#       * miRNA
#           - initial_exclusion
#           - batch - IMPORTANT
#           - slide - IMPORTANT 
#           - final_selection
#           - bio_conc_denat_ngul - IMPORTANT
#           - bio_rin_denat - IMPORTANT
#           - extraction_successful_class
#           - round_rna_extr
#           - sample_milky_blood_class
#           - SampleType
#           
#       * Proteome
#           - plate - IMPORTANT
#           - SampleType
#           
#       * Urina
#           - Sample.ID
#           - Centre
#           - Child_ID
#           - SamplingType - IMPORTANT
#           
#       * Serum
#           - run_order1 - IMPORTANT
#           - run_order2 - IMPORTANT
#           - plate_batch - IMPORTANT
#           - SampleType
#   
#   
#   NOTE: 
#   
#   Cohort, sex and ethnicity are also important but perhaps they can be used in the final 
#   model during the data analysis process

#   SUGGESTED CONFIGURATION:
#       Please add or remove variables as needed
#       
#       * Methylation:
#           denoise_vars <- c("Methyl.nano_conc_ng_ul", "Methyl.extr_batch", "Methyl.Sample_Plate", 
#                             "Methyl.BATCH_Conversion", "Methyl.BATCH_Infinium", "Methyl.Array", 
#                             "Methyl.Slide", "Methyl.SampleType")
# 
#       * Gene Expression:
#           denoise_vars <- c("Gene.round_gexp", "Gene.r_b", "Gene.round_rna_extr", 
#                             "Gene.extr_batch", "Gene.bio_conc_denat_ngul", "Gene.bio_rin_denat", 
#                             "Gene.extraction_successful_class", "Gene.sample_milky_blood_class", "Gene.SampleType")
#   
#       * miRNAs:
#           denoise_vars <- c("miRNA.initial_exclusion", "miRNA.batch", "miRNA.slide", 
#                             "miRNA.final_selection", "miRNA.bio_conc_denat_ngul", "miRNA.bio_rin_denat", 
#                             "miRNA.extraction_successful_class", "miRNA.round_rna_extr", 
#                             "miRNA.sample_milky_blood_class", "miRNA.SampleType")
# 
#       * Proteome:
#           denoise_vars <- c("Prot.plate", "blood_sam4")
#  
#       * Serum metabolome:
#           denoise_vars <- c("Serum.run_order1", "Serum.run_order2", "Serum.plate_batch", "Serum.SampleType")
#  
#       * Urine metabolome:
#           denoise_vars <- c("Urine.Sample.ID", "Urine.Centre", "Urine.Child_ID", "Urine.SamplingType")
#       
#   
 
denoise_vars <- c("Methyl.nano_conc_ng_ul", "Methyl.extr_batch", "Methyl.Sample_Plate","Methyl.BATCH_Conversion", "Methyl.BATCH_Infinium", "Methyl.Array","Methyl.Slide", "Methyl.SampleType")


omicSVA$denoised <- denoising( omicSVA$data, omicsLayer, denoise_vars, ncores = 24)
    


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

num_vars <- c("Methyl.nano_conc_ng_ul","NK","Bcell","CD4T","CD8T","Eos","Mono","Neu","age_sample_months","hs2_zbmi_who_Time1","hs2_height_c_zscore_Time1","e3_sex_Time1","hs2_visit_age_years_Time1","hs2_zsys_bp.v3_2017_Time1","hs2_zdia_bp.v3_2017_Time1")
factor_vars <- c("h_ethnicity_c","cohort","Methyl.extr_batch", "Methyl.Sample_Plate","Methyl.BATCH_Conversion", "Methyl.BATCH_Infinium", "Methyl.Array","Methyl.Slide", "Methyl.SampleType")

# Define the number of first principal components to be used in check, 
# If nComponents is NULL or NA, all components will be used 
# 
#       Example: 
#           nComponents <- 10 # Checks first 10 compoents
#           nComponents <- NA # Checks all the components
#           
nComponents <- 10

check_noised_denoised( original = omicSVA$data, denoised = omicSVA$denoised, 
                       numVars=num_vars, factorVars=factor_vars, 
                       prefix = omicsLayer, topPCA = nComponents)




#-------------------------------------------------------------
#                           SVA 
#-------------------------------------------------------------
# Get SVA variables using different approaches
#-------------------------------------------------------------


# Get formulas for each method
formulas <- get_SVAFormulas(svaMode, interest_var, covariates, to_protect_covariates, class(omicSVA$data))

# Compute SVA for each method
resSVA <- get_SVA( svaMode, formulas, omicSVA, omicsLayer, plot_vars = plot_vars, interest_var = interest_var, prefix="")
    



#-------------------------------------------------------------
#                       Residualization
#-------------------------------------------------------------
# Resudialize data
#-------------------------------------------------------------

# TO BE DEFINED:

# Number of surrotate variables to be used
#   Set the number of surrogate variables (SVs) to be used. 
#       - If n_SVs = NA  then we use all the SVs
#   
#       Example:
#           n_SVs <- NA    # We use all variables
#           n_SVs <- 7     # We use first 7 SVs

n_SVs <- NA

# List of surrogate variables to be residualized. 
#   Set de Surrogate Variables position to be residualized. 
#   IMPORTANT: only the variables set in this variable are residualized. n_SVs is ignored
#       - If SVs_list = NA  or  SVs_list = NULL  ->  this variable is not used
# 
#       Example:
#           SVs_list <- c( 1:3, 5, 7) # Only variables 1,2,3,5 and 7 are residualized
#           SVs_list <- NA
#           
SVs_list <- NA


# Get residualized data
residualized_SVA <- get_SVsResiduals( resSVA = resSVA, n_SVs = n_SVs, dd = omicSVA$data, SVs_list = SVs_list)



#-------------------------------------------------------------
#                       Check PCA
#-------------------------------------------------------------


residualized_PCA <- sapply(names(residualized_SVA[!is.na(residualized_SVA)]), function(method, fulldata, interest_var ){
    
    p <- pca( exprs(fulldata[[method]]), metadata = pData(fulldata[[method]]), removeVar = 0.1)
    
    pdf(paste0(interest_var, "/PCA_", method,".pdf"))
        print(  biplot(p, pointSize = 1) )
    dev.off()
    
    pdf(paste0(interest_var, "/PCA_Pairs_", method,".pdf"))
        print(  pairsplot(p, pointSize = 1) )
    dev.off()
    
    return(p)
    
}, fulldata = residualized_SVA, interest_var = interest_var,  simplify = F)



#-------------------------------------------------------------
#                       Check Association PCA
#-------------------------------------------------------------


sapply(names(residualized_PCA), function(mode) {
    
    sapply( omicSVA$phenoVars, association, model = mode, df = pData(omicSVA$data), svs=residualized_PCA[[mode]]$rotated, curvar = interest_var, maxsvs = n_SVs, PCs = TRUE, simplify = F)
    
}, simplify = F ) 


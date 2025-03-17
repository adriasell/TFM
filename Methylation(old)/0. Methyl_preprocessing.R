################################################################
#################### Methylation  preproc ######################
################################################################

library(lumi)
library(doParallel)
library(minfi)
library(data.table)
library(readxl)
library(xlsx)
library(dplyr)
library(omics)
library(ggplot2)
library(corrplot)
library(isva)
library(SmartSVA)
library(PCAtools)
library(smplot2)
library(DescTools)
library(fastDummies)
library(polycor)

getwd()
setwd("./TFM")
output_path <- "./results/denoising/denoising_methyl"

##············· DATA PREPROCESSING

# ///  Data winsorization 
# Loading omics data
load("./db/ewas/methylome_subcohort_notfitr_v4.Rdata")

# Winsorization functions
getBeta2 <- function( ms, offset  = 0.001) {
  bs <- getBeta( ms )
  bs[ bs == 0 ] <- offset
  bs[ bs == 1 ] <- 1 - offset
  return(bs)
}

process_outliers <- function(bs){
  bs_without_outliers <- bs |> t() |>
    as.data.frame()
  bs_without_outliers[,] <- mclapply(bs_without_outliers,Winsorize) |>
    t()
  return(bs_without_outliers)
}

# Obtaining winsorization results
winsor_m <- methylome_subcohort_notfitr|>
  getBeta2(offset  = 0.001) |>
  process_outliers() |>
  t()

# Inserting winsorizated methylation values into the array
assay(methylome_subcohort_notfitr) <- winsor_m

# /// Filtering with ICC vector (remaining CpGs)
#··· Large vector with low ICCs CpGs 
lowICC <- readRDS("./db/cpgs_low_ICC.RDS")

# Filtering low ICCs
methylome_subcohort_filtered <- methylome_subcohort_notfitr[rownames(methylome_subcohort_notfitr) %in% lowICC,]

# Saving filtered winsorizated methylome
save(methylome_subcohort_filtered, file = paste0(output_path,"/methylome_subcohort_winsor_ICC_filtered.RData"))

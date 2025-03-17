################################# Set R environment 
library(writexl)
library(tidyverse)
library(doParallel)
library(parallel)
library(ggsankey)
library(networkD3)
library(sharp)
library(readxl)
library(dplyr)
library(tidyr)
library(fastDummies)
library(htmlwidgets)
library(webshot)
#install.packages("./packages_R/ggsankey-main", repos = NULL, type = "source")

################################# Defining working directory.
getwd()
setwd("./TFM")
output_path <-"./results/ExWAS"
set.seed(1899)

source("./codi/Functions/functions_ExWAS.R")
source("./codi/Functions/function_volcano_plot.R")

###Prepare the dataset----------------------------------------------------------
#Load data
list_X <- readRDS("./results/RGCCA/list_X.rds")
X_north <- list_X$X_north
X_south <- list_X$X_south
X_combined <- list_X$X_combined

#Load projections
projections <- readRDS("./results/RGCCA/projections.rds")
proj_north <- projections$`Projections N`
proj_south <- projections$`Projections S`
proj_all <-  projections$`Projections all`

#phenotype data
phenotype <- readRDS("./db/pheno/final/bp_wide_validN5332023-10-16.rds")

phenotype <- dummy_cols(phenotype,
                        select_columns = "h_cohort",
                        remove_selected_columns = T,
                        remove_first_dummy = F)

rownames(phenotype) <- phenotype$HelixID
phenotype$e3_sex_Time1 <- as.numeric(phenotype$e3_sex_Time1)
phenotype_all<-phenotype[rownames(phenotype) %in% rownames(proj_all$prot), 
                         c("e3_sex_Time1", "hs2_visit_age_years_Time1", "h_cohort_BIB", "h_cohort_EDEN", "h_cohort_KANC", "h_cohort_MOBA", "h_cohort_INMA")]
phenotype_n<-phenotype[rownames(phenotype) %in% rownames(proj_north$prot),
                       c("e3_sex_Time1", "hs2_visit_age_years_Time1", "h_cohort_BIB", "h_cohort_EDEN", "h_cohort_KANC")]
phenotype_s<-phenotype[rownames(phenotype) %in% rownames(proj_south$prot),
                       c("e3_sex_Time1", "hs2_visit_age_years_Time1", "h_cohort_INMA")]


#Codebook
codebook <- read.csv2("./db/exposome/CODEBOOK_ANALYSIS_AUGUSTO.csv")

#Load exposure data
exposome <- read.csv2("./results/ExWAS/exposome_filtered.csv", row.names = 1)

#Scale exposures
to_factor <- c("h_edumc_None", "h_fish_preg_Ter",  "h_fruit_preg_Ter", "h_legume_preg_Ter",
               "h_veg_preg_Ter", "h_dairy_preg_Ter", "h_meat_preg_Ter","e3_asmokyn_p_None",
               "e3_alcpreg_yn_None")

exposome[, !colnames(exposome) %in% to_factor] <- scale(exposome[, !colnames(exposome) %in% to_factor], center = T)
exposome[to_factor] <- lapply(exposome[to_factor], factor)

exposome_n <- exposome[rownames(exposome) %in% rownames(X_north$prot),]
exposome_s <- exposome[rownames(exposome) %in% rownames(X_south$prot),]
exposome_all <- exposome[rownames(exposome) %in% rownames(X_combined$prot),]

#Prepare df
lat_vars_n <- prepare_df(projections = proj_north, X_data = X_north, exposure_data = exposome_n, response = 4)
lat_vars_s <- prepare_df(projections =  proj_south, X_data = X_south, exposure_data = exposome_s, response = 4)
lat_vars_all <- prepare_df(projections =  proj_all, X_data = X_combined, exposure_data = exposome_all, response = 4)

all(rownames(lat_vars_n)==rownames(phenotype_n))
all(rownames(lat_vars_s)==rownames(phenotype_s))
all(rownames(lat_vars_all)==rownames(phenotype_all))

#ExWAS--------------------------------------------------------------------------
rgcca_validation <- readRDS("./results/RGCCA/model/rgcca_final.rds")
response <- rgcca_validation$call$response
ncomp <- rgcca_validation$call$ncomp

## Association between components and the exposures------------------------------
name_components <- paste0(rep(names(X_combined)[-response],times=ncomp[-response]))

##north
res_exposure_N <- list()
for (comp in name_components){
  form=paste(comp," ~1")
  res_exposure_N[[comp]] <-ExWAS_mixed(data=cbind(lat_vars_n, phenotype_n),
                                     expos_name = colnames(exposome_n),
                                     form=form) %>%  mutate(comp=comp)}

res_exposure_N <- suppressMessages(purrr::reduce(res_exposure_N, full_join)) %>% arrange(p)

##south
exposome_s <- exposome_s[, colSums(exposome_s != 0) > 0]
res_exposure_S <- list()
for (comp in name_components){
  form=paste(comp," ~1")
  res_exposure_S[[comp]] <-ExWAS_mixed(data= cbind(lat_vars_s, phenotype_s),
                                       expos_name = colnames(exposome_s),
                                       form=form) %>%  mutate(comp=comp)}

res_exposure_S <- suppressMessages(purrr::reduce(res_exposure_S, full_join)) %>% arrange(p)

##all
res_exposure_all <- list()
for (comp in name_components){
  form=paste(comp," ~1")
  res_exposure_all[[comp]] <-ExWAS_mixed(data=cbind(lat_vars_all, phenotype_all),
                                          expos_name = colnames(exposome_all),
                                          form=form) %>%  mutate(comp=comp)}

res_exposure_all <- suppressMessages(purrr::reduce(res_exposure_all, full_join)) %>% arrange(p)

#Save results
list_res_exposure <- list("res_exposure_N"=res_exposure_N,
                         "res_exposure_S"=res_exposure_S,
                         "res_exposure_all"=res_exposure_all)

write_xlsx(list_res_exposure, paste0(output_path,"/ExWAS_comp_exposome.xlsx"))
sheets <- excel_sheets(paste0(output_path,"/ExWAS_comp_exposome.xlsx"))
list_res_exposure <- lapply(sheets, read_excel, path = paste0(output_path,"/ExWAS_comp_exposome.xlsx"))
names(list_res_exposure)<-sheets

## Association between components and outcomes----------------------------------
##north
res_outcome_N <- list()
for (outcome in colnames(X_north$Y)){
  form = paste(outcome," ~ 1")
  res_outcome_N[[outcome]]  <-ExWAS_mixed(data=lat_vars_n,
                                          expos_name = name_components,
                                          form=form) %>% mutate(outcome=outcome)}

res_outcome_N <- suppressMessages(purrr::reduce(res_outcome_N, full_join)) %>% arrange(p)

##south
res_outcome_S <- list()
for (outcome in colnames(X_south$Y)){
  form = paste(outcome," ~ 1")
  res_outcome_S[[outcome]]  <-ExWAS_mixed(data=lat_vars_s,
                                          expos_name = name_components,
                                          form=form) %>% mutate(outcome=outcome)}

res_outcome_S <- suppressMessages(purrr::reduce(res_outcome_S, full_join)) %>% arrange(p)

##all
res_outcome_all  <- list()
for (outcome in colnames(X_combined$Y)){
  form = paste(outcome," ~ 1")
  res_outcome_all[[outcome]]  <-ExWAS_mixed(data=lat_vars_all,
                                          expos_name = name_components,
                                          form=form) %>% mutate(outcome=outcome)}

res_outcome_all <- suppressMessages(purrr::reduce(res_outcome_all, full_join)) %>% arrange(p)

#Save results
list_res_outcome <- list("res_outcome_N"=res_outcome_N,
                         "res_outcome_S"=res_outcome_S,
                         "res_outcome_all"=res_outcome_all)

write_xlsx(list_res_outcome, paste0(output_path,"/ExWAS_comp_outcome.xlsx"))
sheets <- excel_sheets(paste0(output_path,"/ExWAS_comp_outcome.xlsx"))
list_res_outcome <- lapply(sheets, read_excel, path = paste0(output_path,"/ExWAS_comp_outcome.xlsx"))
names(list_res_outcome)<-sheets

#Sankeyplot---------------------------------------------------------------------
#Preparar codebook
zdia <- codebook[codebook$Variable_name_TRANS == "hs_zdia_bp",]
zsys <- codebook[codebook$Variable_name_TRANS == "hs_zsys_bp",]

zdia_t1 <- transform(zdia, Variable_name_TRANS = "hs2_zdia_bp.v3_2017_Time1", 
                     description = paste0(zdia$description, " Time1"),
                     Label.short..e.g..for.figures.= "Diastolic BP SD score in childhood")
zdia_t2 <- transform(zdia, Variable_name_TRANS = "hs2_zdia_bp.v3_2017_Time2",
                     description = paste0(zdia$description, " Time2"),
                     Label.short..e.g..for.figures.="Diastolic BP SD score in adolescence")
zsys_t1 <- transform(zsys, Variable_name_TRANS = "hs2_zsys_bp.v3_2017_Time1",
                     description = paste0(zsys$description, " Time1"),
                     Label.short..e.g..for.figures.="Systolic BP SD score in childhood")
zsys_t2 <- transform(zsys, Variable_name_TRANS = "hs2_zsys_bp.v3_2017_Time2",
                     description = paste0(zsys$description, " Time2"),
                     Label.short..e.g..for.figures.="Systolic BP SD score in adolescence")

codebook <- rbind(codebook, zdia_t1, zdia_t2, zsys_t1, zsys_t2)


#Sankey plot North

sankey_res_N <- plotSankey(res_outcome = list_res_outcome$res_outcome_N, 
                           res_exposure = list_res_exposure$res_exposure_N,
                           name_outcomes = colnames(X_north$Y), 
                           variable_labels = codebook, 
                           path = "./results/Sankey_plot/North",
                           p_val = 0.05)

#Sankey plot south
sankey_res_S <-plotSankey(res_outcome = list_res_outcome$res_outcome_S,
                          res_exposure = list_res_exposure$res_exposure_S,
                          name_outcomes = colnames(X_south$Y), 
                          variable_labels = codebook, 
                          path = "./results/Sankey_plot/South",
                          p_val = 0.05)


#Sankey plot all
sankey_res_all <-plotSankey(res_outcome = list_res_outcome$res_outcome_all, 
                            res_exposure = list_res_exposure$res_exposure_all,
                            name_outcomes = colnames(X_combined$Y), 
                            variable_labels = codebook, 
                            path = "./results/Sankey_plot/all",
                            p_val = 0.05)



#LASSO (sharp)------------------------------------------------------------------
exposome <- read.csv2("./results/ExWAS/exposome_filtered.csv", row.names = 1)
to_dummy <- c("h_edumc_None", "h_fish_preg_Ter",  "h_fruit_preg_Ter", "h_legume_preg_Ter",
              "h_veg_preg_Ter", "h_dairy_preg_Ter", "h_meat_preg_Ter")
#Scale exposures
exposome[, !colnames(exposome) %in% to_factor] <- scale(exposome[, !colnames(exposome) %in% to_factor], center = T)
ids<-rownames(exposome)

exposome <- dummy_cols(exposome,
                       select_columns = to_dummy,
                       remove_selected_columns = T,
                       remove_first_dummy = T)
rownames(exposome)<-ids
exposome_n <- exposome[rownames(exposome) %in% rownames(X_north$prot),]
exposome_s <- exposome[rownames(exposome) %in% rownames(X_south$prot),]
exposome_all <- exposome[rownames(exposome) %in% rownames(X_combined$prot),]

#NORTH
res_exp_n<-list()
for (m in 1:length(proj_north)){
  penalty <- c(rep(1,ncol(exposome_n)),rep(0,ncol(phenotype_n)))
  
  res <- VariableSelection(xdata = cbind(exposome_n, phenotype_n),
                           ydata = proj_north[[m]],
                           family = "gaussian", 
                           Lambda = LambdaSequence(lmax =  1e-1, lmin =  1e-3, cardinal = 100),
                           seed = c(12345),
                           penalty.factor = penalty,
                           tau=0.8, 
                           K = 500, 
                           n_cat=3 , 
                           pi_list=seq(0.51, 0.99, by = 0.01), 
                           beep=NULL, 
                           verbose=TRUE)
  
  res_exp_n[[m]] <- res
  names(res_exp_n)[m]<-names(proj_north)[m]
}   

#SOUTH
res_exp_s<-list()
for (m in 1:length(proj_south)){
  penalty <- c(rep(1,ncol(exposome_s)),rep(0,ncol(phenotype_s)))
  
  res <- VariableSelection(xdata = cbind(exposome_s, phenotype_s),
                           ydata = proj_south[[m]],
                           family = "gaussian", 
                           Lambda = LambdaSequence(lmax =  1e-1, lmin =  1e-3, cardinal = 100),
                           seed = c(12345),
                           penalty.factor = penalty,
                           tau=0.8, 
                           K = 500, 
                           n_cat=3 , 
                           pi_list=seq(0.51, 0.99, by = 0.01), 
                           beep=NULL, 
                           verbose=TRUE)
  
  res_exp_s[[m]] <- res
  names(res_exp_s)[m]<-names(proj_south)[m]
}   


#ALL
res_exp_all<-list()
for (m in 1:length(proj_all)){
  penalty <- c(rep(1,ncol(exposome_all)),rep(0,ncol(phenotype_all)))
  
  res <- VariableSelection(xdata = cbind(exposome_all, phenotype_all),
                           ydata = proj_all[[m]],
                           family = "gaussian", 
                           Lambda = LambdaSequence(lmax =  1e-1, lmin = 1e-3, cardinal = 100),
                           seed = c(12345), 
                           penalty.factor = penalty,
                           tau=0.8, 
                           K = 500, 
                           n_cat=3 , 
                           pi_list=seq(0.51, 0.99, by = 0.01), 
                           beep=NULL, 
                           verbose=TRUE)
  
  res_exp_all[[m]] <- res
  names(res_exp_all)[m]<-names(proj_all)[m]
}   

#RESULTS 

#north
for (m in 1:length(res_exp_n)){
  name<- names(res_exp_n[m])
  stab_m <- res_exp_n[[m]]
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
  
  write.csv2(
    cbind("Names"=selected_names,
          "1_Selprop"=1-selected_selprop,
          "Beta"=selected_beta), row.names = FALSE,
    file = paste0("./results/ExWAS/LASSO/Results/n_exp_",name,".csv")
  )
  svg(paste0("./results/ExWAS/LASSO/Calibration_Plots/Cal_Plot_N_", name, ".svg"))
  CalibrationPlot(stab_m)
  dev.off()
  
}


#South
for (m in 1:length(res_exp_s)){
  name<- names(res_exp_s[m])
  stab_m <- res_exp_s[[m]]
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
  
  write.csv2(
    cbind("Names"=selected_names,
          "1_Selprop"=1-selected_selprop,
          "Beta"=selected_beta), row.names = FALSE,
    file = paste0("./results/ExWAS/LASSO/Results/s_exp_",name,".csv")
  )
  svg(paste0("./results/ExWAS/LASSO/Calibration_Plots/Cal_Plot_S_", name, ".svg"))
  CalibrationPlot(stab_m)
  dev.off()
  
}

#all
for (m in 1:length(res_exp_all)){
  name<- names(res_exp_all[m])
  stab_m <- res_exp_all[[m]]
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
  
  write.csv2(
    cbind("Names"=selected_names,
          "1_Selprop"=1-selected_selprop,
          "Beta"=selected_beta), row.names = FALSE,
    file = paste0("./results/ExWAS/LASSO/Results/all_exp_",name,".csv")
  )
  svg(paste0("./results/ExWAS/LASSO/Calibration_Plots/Cal_Plot_all_", name, ".svg"))
  CalibrationPlot(stab_m)
  dev.off()
  }


#PLOTS--------------------------------------------------------------------------
##SANKEY OUTCOME-COMP------------------------------------------------------------

df <- list_res_outcome$res_outcome_all

df <- df %>%
  mutate(outcome = case_when(
    outcome == "hs2_zdia_bp.v3_2017_Time1" ~ "Diastolic BP SD score in childhood",
    outcome == "hs2_zdia_bp.v3_2017_Time2" ~ "Diastolic BP SD score in adolescence",
    outcome == "hs2_zsys_bp.v3_2017_Time1" ~ "Systolic BP SD score in childhood",
    outcome == "hs2_zsys_bp.v3_2017_Time2" ~ "Systolic BP SD score in adolescence",
    TRUE ~ outcome))

df <- df[order(df$variable),]
df$variable <- c(rep(c("Prot-LC"),4), rep(c("Serum-LC"),4), rep(c("Urine-LC"),4))

df <- df %>% filter(p_corrected < 0.05)


nodes <- data.frame(name = unique(c(df$variable, df$outcome)), group = c("#5E0626", "#A7A7A7", "#3C6B66", "#F0B077", "#8CC5E3", "#EA801C", "#3594CC"))
nodes <- nodes[c(1,2,3,5,7,4,6),]

links <- df %>%
  mutate(
    source = match(variable, nodes$name) - 1,
    target = match(outcome, nodes$name) - 1,
    value = abs(beta),
    color = ifelse(beta < 0, "olivedrab", "coral")
  ) %>%
  select(source, target, value, color)

# Sankey plot
sankey <- sankeyNetwork(Links = links,
              Nodes = nodes, 
              Source = "source", 
              Target = "target", 
              Value = "value", 
              NodeID = "name", 
              units = "Beta", 
              fontSize = 20, 
              fontFamily = "serif", 
              nodeWidth = 20,
              sinksRight = F,
              iterations = 0, colourScale <- JS(
  "d3.scaleOrdinal()
  .domain(['#5E0626', '#A7A7A7', '#3C6B66', '#8CC5E3', '#3594CC', '#F0B077', '#EA801C', 'olivedrab', 'coral'])
  .range(['#5E0626', '#A7A7A7', '#3C6B66', '#8CC5E3', '#3594CC', '#F0B077', '#EA801C', 'olivedrab', 'coral'])"
), NodeGroup = "group", LinkGroup = "color")

saveWidget(sankey, file = "./results/ExWAS/Sankey_LC_outcome.html", selfcontained = TRUE)
webshot("./results/ExWAS/Sankey_LC_outcome.html", 
        file = "./results/ExWAS/Sankey_LC_outcome.tiff",
        selector = "body", delay = 1)

##FOREST ExWAS PLOT COMP-EXPOSURE-----------------------------------------------
df_N <- list_res_exposure$res_exposure_N %>% mutate(source = "North")
df_S <- list_res_exposure$res_exposure_S %>% mutate(source = "South")
df_all <- list_res_exposure$res_exposure_all %>% mutate(source = "Pooled")

df_combined <- rbind(df_N, df_S, df_all)
df_combined <- df_combined %>%  group_by(variable, comp, modality) %>%  filter(any(p < 0.05))
df_combined$source <- factor(df_combined$source, 
                             levels = c("North", "South", "Pooled"), 
                             labels = c("North", "South", "Pooled"))

codebook <- codebook %>% rename(Label = `Label.short..e.g..for.figures.`)
df_combined <- inner_join(df_combined, codebook[, 8:9], by = c("variable" = "Variable_name_TRANS"))

df_combined <- df_combined %>%
  mutate(Label = case_when(modality == 1 ~ paste(Label, "- Yes or No"),
                           modality == 2 ~ paste(Label, "intake - Medium vs low"),
                           modality == 3 ~ paste(Label, "intake - High vs low"),
                           TRUE ~ Label))
df_combined$Label[df_combined$Label=="Traffic__100m"]<-"Traffic (100m)"

write.csv2(df_combined, "./results/ExWAS/ExWAS_comp_exposome_fil.csv")

#PROT
colors <- c("North" = "#3A0417",
            "South" = "#A03358",  
            "Pooled" = "#5E0626") 

pdf("./results/ExWAS/Forest_plot/ExWAS_prot.pdf")
ggplot(df_combined[df_combined$comp=="prot",],
       aes(y = reorder(Label, beta), x = beta, xmin = `CI 2.5`, xmax = `CI 97.5`, shape = source, colour = source)) +
  geom_pointrange(size = 0.8, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = colors) + 
  labs(y = "Exposure", x = "Beta (CI 95%)", title = "Prot-LC") +
  guides(shape = guide_legend("Population"), color = guide_legend("Population")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.1) +  
  theme_minimal() +
  theme(text = element_text(size = 14))
dev.off()

#SERUM
colors <- c("North" = "#707070",
            "South" = "#D3D3D3",  
            "Pooled" = "#A7A7A7") 

pdf("./results/ExWAS/Forest_plot/ExWAS_serum.pdf")
ggplot(df_combined[df_combined$comp=="serum",],
       aes(y = reorder(Label, beta), x = beta, xmin = `CI 2.5`, xmax = `CI 97.5`, shape = source, , colour = source)) +
  scale_color_manual(values = colors) + 
  geom_pointrange(size = 0.8, position = position_dodge(width = 0.5)) +
  labs(y = "Exposure", x = "Beta (CI 95%)", title = "Serum-LC") +
  guides(shape = guide_legend("Population"), color = guide_legend("Population")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.1) +  
  theme_minimal() +
  theme(text = element_text(size = 14))
dev.off()

#URINE
colors <- c("North" = "#1F3F3C",
            "South" = "#72A09C",  
            "Pooled" = "#3C6B66") 

pdf("./results/ExWAS/Forest_plot/ExWAS_urine.pdf")
ggplot(df_combined[df_combined$comp=="urine",],
       aes(y = reorder(Label, beta), x = beta, xmin = `CI 2.5`, xmax = `CI 97.5`, shape = source, , colour = source)) +
  scale_color_manual(values = colors) + 
  geom_pointrange(size = 0.8, position = position_dodge(width = 0.5)) +
  labs(y = "Exposure", x = "Beta (CI 95%)", title = "Urine-LC") +
  guides(shape = guide_legend("Population"), color = guide_legend("Population")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.1) +  
  theme_minimal() +
  theme(text = element_text(size = 14))
dev.off()


##FOREST LASSO PLOT COMP-EXPOSURE-----------------------------------------------
all_exp_prot <- read.csv2("results/ExWAS/LASSO/Results/all_exp_prot.csv") %>% mutate(source = "all", comp = "prot")
all_exp_serum <- read.csv2("results/ExWAS/LASSO/Results/all_exp_serum.csv") %>% mutate(source = "all", comp = "serum")
all_exp_urine <- read.csv2("results/ExWAS/LASSO/Results/all_exp_urine.csv") %>% mutate(source = "all", comp = "urine")
n_exp_prot <- read.csv2("results/ExWAS/LASSO/Results/n_exp_prot.csv") %>% mutate(source = "N", comp = "prot")
n_exp_serum <- read.csv2("results/ExWAS/LASSO/Results/n_exp_serum.csv") %>% mutate(source = "N", comp = "serum")
n_exp_urine <- read.csv2("results/ExWAS/LASSO/Results/n_exp_urine.csv") %>% mutate(source = "N", comp = "urine")
s_exp_prot <- read.csv2("results/ExWAS/LASSO/Results/s_exp_prot.csv") %>% mutate(source = "S", comp = "prot")
s_exp_serum <- read.csv2("results/ExWAS/LASSO/Results/s_exp_serum.csv") %>% mutate(source = "S", comp = "serum")
s_exp_urine <- read.csv2("results/ExWAS/LASSO/Results/s_exp_urine.csv") %>% mutate(source = "S", comp = "urine")

df_lasso <- rbind(
  all_exp_prot, 
  all_exp_serum, 
  all_exp_urine, 
  n_exp_prot, 
  n_exp_serum, 
  n_exp_urine, 
  s_exp_prot, 
  s_exp_serum, 
  s_exp_urine)

df_lasso$Beta <- as.numeric(df_lasso$Beta)
df_lasso <- inner_join(df_lasso, codebook[, 8:9], by = c("Names" = "Variable_name_TRANS"))
write.csv2(df_lasso, "./results/ExWAS/LASSO_comp_exposome.csv")

#PLOT 95%????? COM FERHO??
#PROT
png("./results/ExWAS/Forest_plot/LASSO_prot.png", width = 1300, height = 900)
ggplot(df_lasso[df_lasso$comp=="prot",],
       aes(y = reorder(Label, Beta), x = Beta, xmin = Beta, xmax = Beta, shape = source, colour = source)) +
  scale_color_manual(values = colors) +
  geom_pointrange(size = 0.8, position = position_dodge(width = 0.5)) +
  labs(y = "Exposure", x = "Beta (CI 95%)", title = "LASSO Beta exposure-prot component") +
  guides(shape = guide_legend("Population"), color = guide_legend("Population")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.1) +
  theme_minimal() +
  theme(text = element_text(size = 14))
dev.off()

#SERUM
png("./results/ExWAS/Forest_plot/LASSO_serum.png", width = 1300, height = 900)
ggplot(df_lasso[df_lasso$comp=="serum",],
       aes(y = reorder(Label, Beta), x = Beta, xmin = Beta, xmax = Beta, shape = source, colour = source)) +
  scale_color_manual(values = colors) +
  geom_pointrange(size = 0.8, position = position_dodge(width = 0.5)) +
  labs(y = "Exposure", x = "Beta (CI 95%)", title = "LASSO Beta exposure-serum component") +
  guides(shape = guide_legend("Population"), color = guide_legend("Population")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.1) +
  theme_minimal() +
  theme(text = element_text(size = 14))
dev.off()

#URINE
png("./results/ExWAS/Forest_plot/LASSO_urine.png", width = 1300, height = 900)
ggplot(df_lasso[df_lasso$comp=="urine",],
       aes(y = reorder(Label, Beta), x = Beta, xmin = Beta, xmax = Beta, shape = source, colour = source)) +
  scale_color_manual(values = colors) +
  geom_pointrange(size = 0.8, position = position_dodge(width = 0.5)) +
  labs(y = "Exposure", x = "Beta (CI 95%)", title = "LASSO Beta exposure-urine component") +
  guides(shape = guide_legend("Population"), color = guide_legend("Population")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.1) +
  theme_minimal() +
  theme(text = element_text(size = 14))
dev.off()





#ExWAS childhood diet-LC________________________________________________________
exposome1 <- read.csv2("./db/exposome/data_HELIXsubcohort_18_09_23_FILTERED_F.csv")
rownames(exposome1) <- exposome1$HelixID
exposome1 <- exposome1[c("h_bf_None",
                         "h_bfdur_Ter",  
                         "hs_bakery_prod_Ter",  
                         "hs_beverages_Ter",  
                         "hs_break_cer_Ter",  
                         "hs_caff_drink_Ter",  
                         "hs_dairy_Ter",  
                         "hs_fastfood_Ter",  
                         "hs_org_food_Ter",  
                         "hs_proc_meat_Ter",  
                         "hs_readymade_Ter",  
                         "hs_total_bread_Ter",  
                         "hs_total_cereal_Ter",  
                         "hs_total_fish_Ter",  
                         "hs_total_fruits_Ter",  
                         "hs_total_lipids_Ter",  
                         "hs_total_meat_Ter",  
                         "hs_total_potatoes_Ter",  
                         "hs_total_sweets_Ter",  
                         "hs_total_veg_Ter",  
                         "hs_total_yog_Ter",  
                         "h_cereal_post_Log",  
                         "h_dairy_post_Log",  
                         "h_fastfood_post_Log",  
                         "h_fish_post_Log",  
                         "h_fruit_post_Log",  
                         "h_legume_post_Log",  
                         "h_meat_post_Log",  
                         "h_nonalc_post_Log",  
                         "h_sugar_post_Log",  
                         "h_veg_post_Log")]

diet_factors<-c("h_bf_None","h_bfdur_Ter", "hs_bakery_prod_Ter", "hs_beverages_Ter", "hs_break_cer_Ter",
              "hs_caff_drink_Ter", "hs_dairy_Ter", "hs_fastfood_Ter", "hs_org_food_Ter",
              "hs_proc_meat_Ter", "hs_readymade_Ter", "hs_total_bread_Ter", "hs_total_cereal_Ter", "hs_total_fish_Ter",
              "hs_total_fruits_Ter", "hs_total_lipids_Ter", "hs_total_meat_Ter", "hs_total_potatoes_Ter", 
              "hs_total_sweets_Ter", "hs_total_veg_Ter", "hs_total_yog_Ter")

exposome1[, !colnames(exposome1) %in% diet_factors] <- scale(exposome1[,!colnames(exposome1) %in% diet_factors])
exposome1[, colnames(exposome1) %in% diet_factors] <- lapply(exposome1[,diet_factors], as.factor)


exposome1_n <- exposome1[rownames(exposome1) %in% rownames(X_north$prot),]
exposome1_s <- exposome1[rownames(exposome1) %in% rownames(X_south$prot),]
exposome1_all <- exposome1[rownames(exposome1) %in% rownames(X_combined$prot),]

lat_vars_n1 <- prepare_df(projections = proj_north["urine"], X_data = X_north[3:4], exposure_data = exposome1_n, response = 2)
lat_vars_s1 <- prepare_df(projections = proj_south["urine"], X_data = X_south[3:4], exposure_data = exposome1_s, response = 2)
lat_vars_all1 <- prepare_df(projections = proj_all["urine"], X_data = X_combined[3:4], exposure_data = exposome1_all, response = 2)


res_diet_urine_n <- ExWAS_mixed(data=cbind(lat_vars_n1, phenotype_n),
                              expos_name = colnames(exposome1),
                              form=paste("urine"," ~1")) %>% mutate(Group = "North")

res_diet_urine_s <- ExWAS_mixed(data=cbind(lat_vars_s1, phenotype_s),
                                expos_name = colnames(exposome1),
                                form=paste("urine"," ~1")) %>% mutate(Group = "South")

res_diet_urine_all <- ExWAS_mixed(data=cbind(lat_vars_all1, phenotype_all),
                                expos_name = colnames(exposome1),
                                form=paste("urine"," ~1")) %>% mutate(Group = "Pooled")
df_res_diet_urine <- rbind(res_diet_urine_n, res_diet_urine_s, res_diet_urine_all) %>% arrange(p)
write_xlsx(df_res_diet_urine, "./results/ExWAS/ExWAS_UrineLC_dietvars_childhood.xlsx")

#LM pregnancy hg ~ fish_ter_____________________________________________________
#Load data
exposome <- read.csv2("./results/ExWAS/exposome_filtered.csv", row.names = 1)
phenotype <- readRDS("./db/pheno/final/bp_wide_validN5332023-10-16.rds")

rownames(phenotype) <- phenotype$HelixID
exposome_all <- exposome[rownames(exposome) %in% rownames(X_combined$prot),]
phenotype_all <- phenotype[rownames(phenotype) %in% rownames(X_combined$prot),]

all(rownames(exposome_all)==rownames(phenotype_all))
df <- cbind(exposome_all, phenotype_all)
df <- df[, !duplicated(names(df))]

#Linear regression
df$h_fish_preg_Ter <- as.factor(df$h_fish_preg_Ter)
lm_model <- lm(hs_hg_m_Log2 ~ h_fish_preg_Ter + h_cohort + e3_sex.x + hs2_visit_age_years_Time1, data = df)
summary(lm_model)
#summary(lm(hs_hg_m_Log2 ~ h_fish_preg_Ter, data = df))

#Boxplot

ggplot(df, aes(x = factor(h_fish_preg_Ter), y = hs_hg_m_Log2, fill = factor(h_fish_preg_Ter))) +
  geom_boxplot() +
  labs(title = "Mercury Levels by Fish Intake",
       x = "Fish Intake During Pregnancy",
       y = "Log2 Mercury Levels") +
  theme_bw() +
  theme(legend.position = "none") 


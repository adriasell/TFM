require(foreign)
require(nnet)
require(ggplot2)
require(reshape2)
library(dplyr)
library(tidyverse)
library(fastDummies)
library(ggeffects)

setwd("./TFM")

#Load data
pheno <- readRDS("./db/pheno/final/bp_wide_validN5332023-10-16.rds")
rownames(pheno)<-pheno$HelixID
pheno<-pheno[order(rownames(pheno)),]
pheno <- pheno[c("hs2_BPcat_v3_2017_bin_Time2", "HelixID")]

pheno$hs2_BPcat_v3_2017_bin_Time2 <- as.factor(pheno$hs2_BPcat_v3_2017_bin_Time2)

#Load projections 
projections <- readRDS("./results/RGCCA/projections.rds")
proj_north <- projections$`Projections N` %>% purrr::reduce(cbind) %>% data.frame()
colnames(proj_north) <- names(projections$`Projections N`)
proj_south <- projections$`Projections S` %>% purrr::reduce(cbind) %>% data.frame()
colnames(proj_south) <- names(projections$`Projections S`)
proj_all <-  projections$`Projections all` %>% purrr::reduce(cbind) %>% data.frame()
colnames(proj_all) <- names(projections$`Projections all`)

#dataset containing proj and BP trajectories 
pheno_all <- pheno[rownames(pheno) %in% rownames(proj_all),]
pheno_n <- pheno[rownames(pheno) %in% rownames(proj_north),]
pheno_s <- pheno[rownames(pheno) %in% rownames(proj_south),]

all(rownames(proj_all)==rownames(pheno_all))
all(rownames(proj_north)==rownames(pheno_n))
all(rownames(proj_south)==rownames(pheno_s))

proj_all$HelixID <- rownames(proj_all)
proj_north$HelixID <- rownames(proj_north)
proj_south$HelixID <- rownames(proj_south)

df_all <- merge(proj_all, pheno_all, by = "HelixID") %>% mutate(Population = "All")
df_n <- merge(proj_north, pheno_n, by = "HelixID") %>% mutate(Population = "North")
df_s <- merge(proj_south, pheno_s, by = "HelixID") %>% mutate(Population = "South")
 

#Log regression_________________________________________________________
log_all <- glm(hs2_BPcat_v3_2017_bin_Time2 ~ prot + serum + urine, data = df_all, family = binomial())
predictions.log_all_prot <- ggemmeans(log_all, terms = "prot [all]") %>%  rename("Altered" = predicted) %>%  data.frame()
predictions.log_all_serum <- ggemmeans(log_all, terms = "serum [all]") %>% rename("Altered" = predicted) %>% data.frame()
predictions.log_all_urine <- ggemmeans(log_all, terms = "urine [all]") %>% rename("Altered" = predicted) %>% data.frame()

summary(log_all)

#Urine odds ratio prob
1/exp(log_all$coefficients[4])

#
process_predictions <- function(predictions) {
  predictions %>%
    mutate(Normal = 1 - Altered) %>%
    pivot_longer(cols = c(Normal, Altered), names_to = "BP_trajectory", values_to = "probability") %>%
    mutate(conf.low = ifelse(BP_trajectory == "Normal", 1 - conf.low, conf.low),
           conf.high = ifelse(BP_trajectory == "Normal", 1 - conf.high, conf.high),
           BP_trajectory = as.factor(BP_trajectory))
}

predictions.log_all_prot <- process_predictions(predictions.log_all_prot)
predictions.log_all_serum <- process_predictions(predictions.log_all_serum)
predictions.log_all_urine <- process_predictions(predictions.log_all_urine)

predictions.log_all_prot$BP_trajectory<- as.factor(predictions.log_all_prot$BP_trajectory)
predictions.log_all_serum$BP_trajectory<- as.factor(predictions.log_all_serum$BP_trajectory)
predictions.log_all_urine$BP_trajectory<- as.factor(predictions.log_all_urine$BP_trajectory)

#Plots 
colors <- c("Normal"="#089912",
            "Altered"="#950C00")

ggplot(predictions.log_all_prot, aes(x = x, y = probability, fill = BP_trajectory)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = BP_trajectory), alpha = 0.2) + 
  labs(x = "Prot LC", y = "Probabilities",
       title = "Predicted Probabilities of having BP altered in adolescence") +
  scale_fill_manual(values = colors) +
  theme_minimal()


ggplot(predictions.log_all_serum, aes(x = x, y = probability, fill = BP_trajectory)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = BP_trajectory), alpha = 0.2) + 
  labs(x = "Serum LC", y = "Probabilities",
       title = "Predicted Probabilities of having BP altered in adolescence") +
  scale_fill_manual(values = colors) +
  theme_minimal()

pdf("./results/Log_regression/logreg_urine.pdf")
ggplot(predictions.log_all_urine, aes(x = x, y = probability, fill = BP_trajectory)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = BP_trajectory), alpha = 0.2) + 
  labs(x = "Urine LC", y = "Probabilities",
       title = "Predicted Probabilities of having BP altered in adolescence",
       fill = "Altered adolescence BP (>90 perc)") +
  scale_fill_manual(values = colors) +
  theme_minimal()
dev.off()

#North
log_n <- glm(hs2_BPcat_v3_2017_bin_Time2 ~ prot + serum + urine, data = df_n, family = binomial())
predictions.log_n_prot <- ggemmeans(log_n, terms = "prot [all]") %>% rename("Altered" = predicted) %>% data.frame()
predictions.log_n_serum <- ggemmeans(log_n, terms = "serum [all]") %>% rename("Altered" = predicted) %>% data.frame()
predictions.log_n_urine <- ggemmeans(log_n, terms = "urine [all]") %>% rename("Altered" = predicted) %>% data.frame()
summary(log_n)

predictions.log_n_prot <- process_predictions(predictions.log_n_prot)
predictions.log_n_serum <- process_predictions(predictions.log_n_serum)
predictions.log_n_urine <- process_predictions(predictions.log_n_urine)


ggplot(predictions.log_n_prot, aes(x = x, y = probability, fill = BP_trajectory)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = BP_trajectory), alpha = 0.2) + 
  labs(x = "Prot LC", y = "Probabilities", title = "Predicted Probabilities of having BP altered in adolescence") +
  scale_fill_manual(values = colors) +
  theme_minimal()

ggplot(predictions.log_n_serum, aes(x = x, y = probability, fill = BP_trajectory)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = BP_trajectory), alpha = 0.2) + 
  labs(x = "Serum LC", y = "Probabilities",
       title = "Predicted Probabilities of having BP altered in adolescence") +
  scale_fill_manual(values = colors) +
  theme_minimal()


ggplot(predictions.log_n_urine, aes(x = x, y = probability, fill = BP_trajectory)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = BP_trajectory), alpha = 0.2) + 
  labs(x = "Urine LC", y = "Probabilities",
       title = "Predicted Probabilities of having BP altered in adolescence") +
  scale_fill_manual(values = colors) +
  theme_minimal()


#South
log_s <- glm(hs2_BPcat_v3_2017_bin_Time2 ~ prot + serum + urine, data = df_s, family = binomial())
predictions.log_s_prot <- ggemmeans(log_s, terms = "prot [all]") %>% rename("Altered" = predicted) %>% data.frame()
predictions.log_s_serum <- ggemmeans(log_s, terms = "serum [all]") %>% rename("Altered" = predicted) %>% data.frame()
predictions.log_s_urine <- ggemmeans(log_s, terms = "urine [all]") %>% rename("Altered" = predicted) %>% data.frame()
summary(log_s)

predictions.log_s_prot <- process_predictions(predictions.log_s_prot)
predictions.log_s_serum <- process_predictions(predictions.log_s_serum)
predictions.log_s_urine <- process_predictions(predictions.log_s_urine)

ggplot(predictions.log_s_prot, aes(x = x, y = probability, fill = BP_trajectory)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = BP_trajectory), alpha = 0.2) + 
  labs(x = "Prot LC", y = "Probabilities", title = "Predicted Probabilities of having BP altered in adolescence") +
  scale_fill_manual(values = colors) +
  theme_minimal()

ggplot(predictions.log_s_serum, aes(x = x, y = probability, fill = BP_trajectory)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = BP_trajectory), alpha = 0.2) + 
  labs(x = "Serum LC", y = "Probabilities",
       title = "Predicted Probabilities of having BP altered in adolescence") +
  scale_fill_manual(values = colors) +
  theme_minimal()


ggplot(predictions.log_s_urine, aes(x = x, y = probability, fill = BP_trajectory)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = BP_trajectory), alpha = 0.2) + 
  labs(x = "Urine LC", y = "Probabilities",
       title = "Predicted Probabilities of having BP altered in adolescence") +
  scale_fill_manual(values = colors) +
  theme_minimal()

#Cases for each population
table(df_all$hs2_BPcat_v3_2017_bin_Time2)
 table(df_n$hs2_BPcat_v3_2017_bin_Time2)
table(df_s$hs2_BPcat_v3_2017_bin_Time2)


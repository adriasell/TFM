
#Table 2 Metrics________________________________________________________________
rgcca_val <- readRDS("./results/RGCCA/model_results/results_rgcca_validation.rds")
rgcca_N_all <- readRDS("./results/RGCCA/model_results/results_rgcca_trainN_testall.rds")

metric_N <- rgcca_val$metrics$train %>% mutate(Population = "N")
metric_S <- rgcca_val$metrics$test %>% mutate(Population = "S")
metric_all <- rgcca_N_all$`RGCCA pred (X_pred)`$metric$test %>% mutate(Population = "Pooled")

metrics_combined <- rbind(metric_N, metric_S, metric_all)
metrics_combined$Population <- factor(metrics_combined$Population,
                                      levels = c("N", "S", "Pooled"), 
                                      labels = c("North", "South", "Pooled"))

metrics_combined <- metrics_combined[order(rownames(metrics_combined)), ] %>% t() %>% data.frame()
metrics_combined1 <- rbind(metrics_combined["Population",], as.data.frame(lapply(metrics_combined[-5, ], function(x) round(as.numeric(x), 3)))) %>%
  as.data.frame()

rownames(metrics_combined1)[2:5] <- c("Diastolic BP SD score in childhood", "Systolic BP SD score in childhood",
                                     "Diastolic BP SD score in adolescence", "Systolic BP SD score in adolescence")


write.csv2(metrics_combined1,"./results/Paper_figures/metrics.csv")

#Supplementary figure____________________________________________________________
#Hierarchical clustering and correlation LC & biomarkers selected

library(corrplot)
library(pheatmap)

set.seed(1999)

#Load data----------------------------------------------------------------------
#Latent vars
latent_vars <- read.csv2("./results/RGCCA/model_results/lat_vars.csv", row.names = 1)
colnames(latent_vars) <- c("Prot-LC", "Serum-LC", "Urine-LC")

#Biomarkers
list_X <- readRDS("./results/RGCCA/list_X.rds")
X_combined <- list_X$X_combined
df_biomarkers <- cbind(X_combined$prot[c("IL6", "IL1beta", "HGF", "BAFF", "TNFalfa", "IL8")],
                       X_combined$serum["log.PC.aa.C38.3"],
                       X_combined$urine[c("p.cresol.sulfate", "X3.Indoxylsulfate")]) %>% scale()

#Projections
proj <- readRDS("./results/RGCCA/projections.rds")
projections <-proj$`Projections all` %>% purrr::reduce(cbind)
colnames(projections) <- c("Prot-LC", "Serum-LC", "Urine-LC")

#Plots--------------------------------------------------------------------------
#Corr plot LC
pdf("./results/Paper_figures/Extra/correlation_LC.pdf")
corrplot.mixed(cor(latent_vars), upper = "ellipse", lower = "number",
               tl.pos = "lt", tl.col = "black", tl.offset=1, tl.srt = 40, 
               tl.cex = 0.5, number.cex = 0.5)
dev.off()

#Corr plot selected biomarkers
pdf("./results/Paper_figures/Extra/correlation_biomarkers.pdf")
corrplot.mixed(cor(scale(df_biomarkers)), upper = "ellipse", lower = "number",
               tl.pos = "lt", tl.col = "black", tl.offset=1, tl.srt = 40, 
               tl.cex = 0.5, number.cex = 0.5)
dev.off()

#Hierarchical clustering LC
pdf("./results/Paper_figures/Extra/hclust_LC.pdf")
pheatmap(t(projections), show_colnames = F)
dev.off()

#Hierarchical clustering biomarkers+
pdf("./results/Paper_figures/Extra/hclust_biomarkers.pdf")
pheatmap(t(df_biomarkers),show_colnames = F)
dev.off()



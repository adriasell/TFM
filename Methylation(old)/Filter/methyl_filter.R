library(tidyverse)

setwd("./TFM")
helixid_n <- c(read.csv("./helixid_n.csv", row.names = 1)$x)
helixid_s <- c(read.csv("./helixid_s.csv", row.names = 1)$x)

#Load methylation data
omics<-readRDS("./results/denoising/denoising_methyl/methyl_denoised.RDS")

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

df<-data.frame(t(getBeta(omics)))
rownames(df) <- pData(omics)$HelixID.x
fil_df <- df[,colnames(df) %in% bibliography_cpgs, drop = F]
fil_df$HelixID <- rownames(fil_df)
fil_df <- fil_df %>% arrange(HelixID)
write.csv2(fil_df, "./results/denoising/denoising_methyl/fil_methyl_all.csv", row.names = T)

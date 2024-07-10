####  required packages
library(dplyr)
library(dbplyr)
library(survival)
library(ggplot2)
library(ggfortify)
library(jsonlite)
library(ggrepel)
library(ggforce)
library(plotly)
library(ggplot2)
library(xlsx)
library(survMisc)
library("stringr")
library("survminer")


####  set paths
melanodb_path <- '../database/MelanoDB.db'

#### load database
melanodb <- DBI::dbConnect(RSQLite::SQLite(), melanodb_path) 
src_dbi(melanodb) 

#### create dataframes

df_patients <- tbl(melanodb, 'clinical')

df_vars = data.frame(df_patients %>% select(sex, age, BOR, source))
df_vars$source <- sapply(df_vars$source, function(x) strsplit(strsplit(x, ',')[[1]][2], ':')[[1]][2])
df_vars$source <- sapply(df_vars$source, function(x) str_replace(x, '"', ''))

df_pfs = data.frame(df_patients %>% select(PFS_month, PFS_statut, source))
df_pfs$source <- sapply(df_pfs$source, function(x) strsplit(strsplit(x, ',')[[1]][2], ':')[[1]][2])
df_pfs$source <- sapply(df_pfs$source, function(x) str_replace(x, '"', ''))
df_pfs$PFS_statut <- sapply(df_pfs$PFS_statut, as.numeric)
df_test <- na.omit(df_pfs)

df_os = data.frame(df_patients %>% select(OS_month, OS_statut, source))
df_os$source <- sapply(df_os$source, function(x) strsplit(strsplit(x, ',')[[1]][2], ':')[[1]][2])
df_os$source <- sapply(df_os$source, function(x) str_replace(x, '"', ''))df_pfs$PFS_statut <- sapply(df_pfs$PFS_statut, as.numeric)
df_os$OS_statut[df_os$OS_statut == 'dead'] <- 1
df_os$OS_statut[df_os$OS_statut == 'alive'] <- 0
df_os$OS_statut <- sapply(df_os$OS_statut, as.numeric)


#### Statistics on survival analysis curves ####
df_pfs$source <- factor(df_pfs$source)

#logrank test
surv_diff <- survdiff(Surv(PFS_month, PFS_statut) ~ source, data = df_pfs) 
surv_diff$pvalue #pvalue = 0.00024

# pairwise test - observation of differential curves
res <- try(pairwise_survdiff(Surv(PFS_month, PFS_statut) ~ source, data = na.omit(df_pfs)))
symnum(res$p.value, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
       symbols = c("****", "***", "**", "*", "+", " "),
       abbr.colnames = FALSE, na = "")
df_pfs.fit <- ten(Surv(PFS_month, PFS_statut) ~ source, data = na.omit(df_pfs))
comp(df_pfs.fit, p=0, q=0, scores=1:6)

#logrank test removing Louveau & Van allen patients
df_pfs_reduced <- df_pfs[!(df_pfs$source %in% c('Baptiste Louveau', 'Eliezer M. Van Allen')),]
surv_diff_red <- survdiff(Surv(PFS_month, PFS_statut) ~ source, data = df_pfs_reduced, rho=0)
surv_diff_red$pvalue 


## OS
df_os$source <- factor(df_os$source)
#logrank test
surv_diff_os <- survdiff(Surv(OS_month, OS_statut) ~ source, data = df_os) 
surv_diff_os$pvalue #pvalue = 0.00027

# pairwise - observation of differential curves
res <- try(pairwise_survdiff(Surv(OS_month, OS_statut) ~ source, data = na.omit(df_os)))
symnum(res$p.value, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
       symbols = c("****", "***", "**", "*", "+", " "),
       abbr.colnames = FALSE, na = "")
df_os.fit <- ten(Surv(OS_month, OS_statut) ~ source, data = na.omit(df_os))
comp(df_os.fit, p=0, q=0, scores=1:6)
#logrank test - removing Louveau & Van allen patients
df_os_reduced <- df_os[!(df_os$source %in% c('"Baptiste Louveau', '"Yibing Yan')),]
surv_diff_os_red <- survdiff(Surv(OS_month, OS_statut) ~ source, data = df_os_reduced, rho=0)
surv_diff_os_red$pvalue 


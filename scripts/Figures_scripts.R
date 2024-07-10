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
library(xlsx)
library(survMisc)
library(BatchQC)
library(sva)
library(stats)
library(ComplexHeatmap)
library(tibble)
library(circlize)
library(tidyverse)
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

#### Pie chart & Boxplot plots ####
##### SOURCE Distribution
df_patients_effectives <- df_patients %>% select(patientID, drug, BOR, BRAF_mut, source) %>% group_by(source) %>% summarize(NumberOfPatients=n_distinct(patientID)) %>% as.data.frame()
get_author <- function(x){
  return(fromJSON(x)$author)
}
df_patients_effectives$source <- sapply(df_patients_effectives$source, get_author)
myPalette = c("#6d597a", "#d5bdaf", "#E76F51", "#C0DC7E", "#F4A261", "#e63946", "#2A9D8F", "#264653", "#E9C46A")
pieLabels = paste(df_patients_effectives$source, "\n", df_patients_effectives$NumberOfPatients)
pie(c(df_patients_effectives$NumberOfPatients), labels = pieLabels, border="white", col=myPalette, radius=1, cex=1.5)

###### DRUG Distribution
df_drugs <- df_patients %>% select(patientID, drug, BOR, BRAF_mut, source) %>% group_by(drug) %>% summarize(NumberOfPatients=n_distinct(patientID)) %>% as.data.frame()
myPalette <- c("#19D3F3", "#AB63FA", "#FF6692", "#FF6692", "#00CC96", "#FEA15A")
pieLabels = paste(df_drugs$drug, "\n", df_drugs$NumberOfPatients)
pie(c(df_drugs$NumberOfPatients), labels = pieLabels, border="white", col=myPalette, radius=1.05,cex=2)

###### BOR Distribution
df_dcr <- df_patients %>% select(patientID, drug, BOR, BRAF_mut, source) %>% group_by(BOR) %>% summarize(NumberOfPatients=n_distinct(patientID)) %>% as.data.frame()
myPalette <- c("#81B29A", "#F4F1DE", "#7469FC", "#E07A5F", "#3D405B", "#00CC96", "#F2CC8F")
pieLabels = paste(df_dcr$BOR, "\n", df_dcr$NumberOfPatients)
pie(c(df_dcr$NumberOfPatients), labels = pieLabels, border="white", col=myPalette, cex=2)

###### AGE Distribution
ggboxplot(df_vars, x = "source", y = "age", fill = "source",)+
  theme(legend.text=element_text(size=18), text=element_text(size=13.2))+
  scale_fill_manual(values=c("#264653", "#2A9D8F", "#E9C46A", "#F4A261", "#E76F51", "#d5bdaf", "#e63946", "#C0DC7E", "#6d597a"))

#### Survival plots ####

###### PFS 
km_source_fit <- survfit(Surv(PFS_month, PFS_statut) ~ source, data=df_pfs)
p <- ggsurvplot(km_source_fit, data = df_pfs, 
                surv.median.line = "hv",# Add medians survival
                size = 0.9,
                legend.title = "Source",
                ylab="Progression Free-Disease probability",
                xlab="Time in months",
                font.tickslab = c(16, face = "bold"),
                font.x = c(20),
                font.y = c(20),
                font.legend = c(20),
                tables.height = 0.2,
                palette = c("#E76F51", "#E9C46A", "#2A9D8F","#8000FF","#C0DC7E","#264653","#F4A261"),
                ggtheme=theme_bw(),
                pval = TRUE,
                conf.int = FALSE)
#risk.table = TRUE,)
p$plot <- p$plot + 
  scale_x_continuous(expand = c(0, 0, .05, 0)) +
  scale_y_continuous(expand = c(0, 0, .05, 0))
p

###### OS 
kmos_source_fit <- survfit(Surv(OS_month, OS_statut) ~ source, data=df_os)
q <- ggsurvplot(kmos_source_fit, data = df_os, 
                surv.median.line = "hv",# Add medians survival
                pval=TRUE,
                conf.int = FALSE,
                size = 0.8,
                legend.title = "Source",
                xlab="Time in months",
                font.tickslab = c(16, face = "bold"),
                font.x = c(20),
                font.y = c(20),
                font.legend = c(20),
                tables.height = 0.2,
                palette = c("#E76F51", "#2A9D8F", "#E63946","#D5BDAF","#C0DC7E","#264653","#F4A261"),
                ggtheme=theme_bw())
q$plot <- q$plot + 
  scale_x_continuous(expand = c(0, 0, .05, 0)) +
  scale_y_continuous(expand = c(0, 0, .05, 0))
q


#### Mutation landscape ####
df_snvs = read.csv('./inputs/melanodb_snv.csv')
df_snvs_infos = read.csv('./inputs/melanodb_snv_infos.csv')
df_snvs <- df_snvs %>% column_to_rownames('X')
mat <- as.matrix(df_snvs)
mat[,1] <- as.numeric(mat[,1])

df_snvs_infos <- filter(df_snvs_infos, patientID != "VS_Pat_29")
df_snvs_infos <- df_snvs_infos[, c("source")]
set.seed(1)
annot = HeatmapAnnotation(df = df_snvs_infos, col = list(type = c("Pauline Blateau, Jerome Solassol" =  rgb(red=0.14, green=0.27, blue=0.32), 
                                                                  "Federica Catalanotti, David B. Solit" = rgb(red=0.16, green=0.61, blue=0.56),
                                                                  "Eliezer M. Van Allen, Dirk Schadendorf"= rgb(red=0.96, green=0.81, blue=0.23),
                                                                  "Baptiste Louveau, Samia Mourah"= rgb(red=0.95, green=0.63, blue=0.38))))
color1 <- rgb(red=0.14, green=0.27, blue=0.32)
color2 <- rgb(red=0.16, green=0.61, blue=0.56)
color3 <- rgb(red=0.96, green=0.81, blue=0.23)
color4 <- rgb(red=0.95, green=0.63, blue=0.38)
#annot <- colorRamp2(c("Pauline Blateau, Jerome Solassol", "Federica Catalanotti, David B. Solit", "Eliezer M. Van Allen, Dirk Schadendorf", "Baptiste Louveau, Samia Mourah"), 
#                    c(color1, color2, color3, color4), space = "RGB")

patient_annotation = HeatmapAnnotation(df = df_snvs_infos, col = list(cell.source=brain.cell.cols.assigned))

Heatmap(mat, 
        show_row_names = TRUE,
        col = colorRamp2(c(0, 1), c("white", "black")),
        rect_gp = gpar(col="grey"),
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        show_heatmap_legend = FALSE,
        clustering_distance_columns = "euclidean",
        #clustering_method_columns ="ward.D",
        column_dend_gp = gpar(),
        row_dend_gp = gpar(),
        top_annotation = annot,
        column_names_gp = gpar(fontsize = 8))

#### Clinical features and treatment outcomes summary table ####

# Features to retrieve: total number of patients, sex, age, Lactate dehydrogenase, type of drug, best overall response, progression-free disease, overall survival, vital status 

df_patients_features <- df_patients %>% select(patientID, sex, age, LDH, drug, BOR, PFS_month, OS_month, OS_statut, source) %>% as.data.frame()
get_author <- function(x){
  return(fromJSON(x)$author)
}
df_patients_features$source <- sapply(df_patients_features$source, get_author)

to_excel <- function(dataframe){
  write.csv(dataframe,file=paste0('./Effectives_counts/', substitute(dataframe), '.csv'))
}

##### Global effectives

summary_tot_sex <- df_patients_features %>% group_by(sex) %>% summarize(NumberOfPatients=n_distinct(patientID)) %>% mutate(perc = 100*NumberOfPatients/sum(NumberOfPatients)) %>% as.data.frame()
to_excel(summary_tot_sex)

summary_tot_LDH <- df_patients_features %>% group_by(LDH) %>% summarize(NumberOfPatients=n_distinct(patientID)) %>% mutate(perc = 100*NumberOfPatients/sum(NumberOfPatients)) %>% as.data.frame()
to_excel(summary_tot_LDH)

summary_tot_drug <- df_patients_features %>% group_by(drug) %>% summarize(NumberOfPatients=n_distinct(patientID)) %>% mutate(perc = 100*NumberOfPatients/sum(NumberOfPatients)) %>% as.data.frame()
to_excel(summary_tot_drug)

summary_tot_BORR <- df_patients_features %>% group_by(BOR) %>% summarize(NumberOfPatients=n_distinct(patientID)) %>% mutate(perc = 100*NumberOfPatients/sum(NumberOfPatients)) %>% as.data.frame()
to_excel(summary_tot_BORR)

summary_tot_age <- df_patients_features %>% mutate(age = as.numeric(age)) %>% summarise(across(age, list(min=min, max=max, median=median)))
to_excel(summary_tot_age)

summary_tot_pfs <- df_patients_features %>% mutate(pfs = as.numeric(PFS_month)) %>% summarise(across(pfs, list(min=min, max=max, median=median)))
to_excel(summary_tot_pfs)

summary_tot_os <- df_patients_features %>% mutate(os = as.numeric(OS_month)) %>% summarise(across(os, list(min=min, max=max, median=median)))
to_excel(summary_tot_os)

summary_tot_os_stat <- df_patients_features %>% group_by(OS_statut) %>% summarize(NumberOfPatients=n_distinct(patientID)) %>% mutate(perc = 100*NumberOfPatients/sum(NumberOfPatients)) %>% as.data.frame()
to_excel(summary_tot_os_stat)


##### Effectives by studies
summary_sex <- df_patients_features %>% group_by(source, sex) %>% summarise(NumberOfPatients=n_distinct(patientID)) %>% mutate(perc = 100*NumberOfPatients/sum(NumberOfPatients)) %>% as.data.frame()
to_excel(summary_sex)

summary_LDH <- df_patients_features %>% group_by(source, LDH) %>% summarize(NumberOfPatients=n_distinct(patientID)) %>% mutate(perc = 100*NumberOfPatients/sum(NumberOfPatients)) %>% as.data.frame()
to_excel(summary_LDH)

summary_drug <- df_patients_features %>% group_by(source, drug) %>% summarize(NumberOfPatients=n_distinct(patientID)) %>% mutate(perc = 100*NumberOfPatients/sum(NumberOfPatients)) %>% as.data.frame()
to_excel(summary_drug)

summary_BORR <- df_patients_features %>% group_by(source, BOR) %>% summarize(NumberOfPatients=n_distinct(patientID)) %>% mutate(perc = 100*NumberOfPatients/sum(NumberOfPatients)) %>% as.data.frame()
to_excel(summary_BORR)

summary_age <- df_patients_features %>% group_by(source) %>% mutate(age = as.numeric(age)) %>% summarise(across(age, list(min=min, max=max, median=median)))
to_excel(summary_age)

summary_pfs <- df_patients_features %>% group_by(source) %>% mutate(pfs = as.numeric(PFS_month)) %>%  summarise(across(pfs, list(min=min, max=max, median=median)))
to_excel(summary_pfs)

summary_os <- df_patients_features %>% group_by(source) %>% mutate(os = as.numeric(OS_month)) %>% summarise(across(os, list(min=min, max=max, median=median)))
to_excel(summary_os)

summary_os_stat <- df_patients_features %>% group_by(source, OS_statut) %>% summarize(NumberOfPatients=n_distinct(patientID)) %>% mutate(perc = 100*NumberOfPatients/sum(NumberOfPatients)) %>% as.data.frame()
to_excel(summary_os_stat)

tables_effectives <- c(summary_tot_sex, 
                       summary_tot_LDH, 
                       summary_tot_drug, 
                       summary_tot_BORR, 
                       summary_tot_age, 
                       summary_tot_pfs, 
                       summary_tot_os, 
                       summary_tot_os_stat,
                       summary_sex, 
                       summary_LDH, 
                       summary_drug,
                       summary_BORR,
                       summary_age,
                       summary_pfs, 
                       summary_os,
                       summary_os_stat)

for (tab in tables_effectives)
{
  print(substitute(tab))
}


#### BatchQC plots #####
melanodb_data <- read.csv('./inputs/melanodb_ge.csv')
row.names(melanodb_data) <- melanodb_data$X
melanodb_data <- melanodb_data[-1]
melanodb_data <- na.omit(melanodb_data)

melanodb_data <- as.data.frame(melanodb_data)
colnames(melanodb_data) <- gsub("^X", "", colnames(melanodb_data))
colnames(melanodb_data) <- gsub("\\.", "-", colnames(melanodb_data))

#melanodb_data <- melanodb_data %>% select(-samples_to_remove)
melanodb_info <- read.csv('./inputs/melanodb_ge_infos.csv')

t <- data.matrix(melanodb_data)

batchQC(t, batch=melanodb_info$Batch, condition=melanodb_info$category,
        report_file="./outputs/batchqc_report.html", report_dir=".",
        report_option_binary="111111111",
        view_report=FALSE, interactive=TRUE, batchqc_output=TRUE)

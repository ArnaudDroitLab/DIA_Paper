
## LIBRARIES --------------------------------------------
library(dplyr)
library(plyr)
library(tidyverse)
library(data.table)

# Store your exports in a directory called "Input_1_Skyline/Results"
setwd("Z:/resultats/RD/RD139_LogicielsDIA/RD139_Papier_R_210108/")


### Skyline ----
## Load files from Skyline 
all_files <- list.files(path = "./Input_1_Skyline/Results/", pattern = "Library_Export.csv", full.names = T)
l <- lapply(all_files, fread, sep = ",") 
all_files <- gsub("./Input_1_Skyline/Results/|_Export.csv", "", all_files)

for(i in 1:length(all_files)) {
  colnames(l[[i]]) <- gsub(".*UPS1_", "UPS1_", colnames(l[[i]])) %>%
    gsub("inj1r", "inj1", .) %>%
    gsub("0_1fmol", "0.1fmol", .) %>%
    gsub("0_25fmol", "0.25fmol", .) %>%
    gsub("2_5fmol", "2.5fmol", .) %>%
    gsub(" ", ".", .)
  
  l[[i]] <- l[[i]] %>%
    mutate(Sequence = Peptide.Modified.Sequence.Three.Letter.Codes) %>%
    dplyr::select(-c(contains("Library.Dot"), contains("Average.Mass"))) %>%
    dplyr::select(Accession = Protein.Name, ModifiedSequence = Peptide.Modified.Sequence.Three.Letter.Codes,
                  Sequence, Precursor = Precursor.Charge,
                  contains("_0.1fmol_"), contains("_0.25fmol_"), contains("_1fmol_"), contains("_2.5fmol_"), 
                  contains("_5fmol_"), contains("_10fmol_"), contains("_25fmol_"), contains("_50fmol_")) 
}



dt <- rbindlist(l, use.names = TRUE, idcol = "Experiment")  %>%
  as.data.frame %>%
  mutate(Experiment = ifelse(Experiment == 1, all_files[1], 
                             ifelse(Experiment == 2, all_files[2], 
                             ifelse(Experiment == 3, all_files[3], all_files[4]))) )


## Filtering
df <- dt %>%
  filter(Accession != "Decoys",
         Accession != "Biognosys standards", 
         Accession != "P07339ups|CATD_HUMAN_UPS") 

Area <- seq(6, ncol(df), by = 2) 
Qval <- seq(7, ncol(df), by = 2) 

for (i in 1:length(Area)) { 
  df[,Area[i]] <- ifelse(df[,Qval[i]] <= 0.01, df[,Area[i]], NA)
}

df <- df %>% 
  dplyr::select(-contains("Detection.Q.Value")) 
colnames(df) <- gsub(".Total.Area", "", colnames(df))


# Reformate table
df <- df %>%
  mutate(specie = str_extract(Accession, "HUMAN|ECOLI"), 
         id = if_else(specie == "HUMAN", 
                      str_extract(Accession, "[A-Z][0-9]*"), 
                      str_extract(Accession, "\\|[^\\|]*\\|") %>% str_replace_all("\\|", "")), 
         Accession = ifelse(specie == "HUMAN", paste0(id, "_HUMAN"), paste0(id, "_ECOLI")))  %>%     
  dplyr::select(-id, -specie) %>%
  mutate(Sequence = gsub("\\[CAM\\]", "", Sequence) %>%
           gsub("\\[Oxi\\]", "", .) %>%
           gsub("\\[Dea\\]", "", .))


# Change modified sequence
transition <- read_csv("Input_1_Other/transition.csv") %>% unique
for (y in transition$old) {
  new <- filter(transition, old == y) %>% pull(new)
  stopifnot(length(new)==1)
  df$ModifiedSequence <- str_replace_all(df$ModifiedSequence, y, new)         
}


# Separate tables for each window size 
pep <- list()
for (i in 1:length(all_files)) {
  name <- unique(df$Experiment)[i]
  pep[[i]] <- df %>%
    filter(Experiment == name)
}


# Normalization 
dat_Norm = dat_median_all = factor_norm = dat_Norm = list()
for (i in 1:length(all_files)) {
  dat_Norm[[i]] <- pep[[i]] %>%
    dplyr::select(contains("UPS1_"))
  dat_median_all[[i]] <- apply(dat_Norm[[i]], 2, function(x) { median(x, na.rm = T) })
  factor_norm [[i]] <- dat_median_all[[i]] / median(dat_median_all[[i]])
  dat_Norm[[i]] <- t(apply(dat_Norm[[i]], 1, "/", factor_norm [[i]])) %>% as.data.frame
}


# Export for BoxPlot
for (i in 1:length(all_files)) {
  write.table(pep[[i]][-c(1:5)], paste0("./Output_2_DataTreatment/Normalisation/", all_files[i], "_BeforeNorm.txt"))
  write.table(dat_Norm[[i]], paste0("./Output_2_DataTreatment/Normalisation/", all_files[i], "_AfterNorm.txt"))
}


# Missing values imputation
dat_quant = dat_Norm_Imp = list()
for (i in 1:length(all_files)) {
  
  dat_quant[[i]] <- apply(dat_Norm[[i]], 2, function(x) {quantile(x, probs = 0.01, na.rm = T)})  
  dat_Norm_Imp[[i]] <- dat_Norm[[i]]
  
  for(j in 1:length(dat_quant[[i]])) {
    dat_Norm_Imp[[i]][is.na(dat_Norm_Imp[[i]][,j]),j] <- dat_quant[[i]][j]
  }
  
}



# Add row_to_filter columns for each triplicate 
x <- c("_0.1fmol_", "_0.25fmol_", "_1fmol_", "_2.5fmol_","_5fmol_", "_10fmol_", "_25fmol_", "_50fmol_")

for (i in 1:length(all_files)) {
  
  
  dat_Norm_Imp[[i]] <- cbind(dat_Norm_Imp[[i]], matrix(nrow = nrow(dat_Norm_Imp[[i]]), ncol=length(x)))
  colnames(dat_Norm_Imp[[i]])[25:ncol(dat_Norm_Imp[[i]])] <- c("RtoF_0.1fmol", "RtoF_0.25fmol", "RtoF_1fmol", "RtoF_2.5fmol", "RtoF_5fmol", "RtoF_10fmol", "RtoF_25fmol", "RtoF_50fmol")
  
  for (j in 1:length(x)) {     
    dat_Norm_Imp[[i]][24+j] <- apply( dat_Norm[[i]], 1, function(row) 
    { length(which(is.na( row[grepl(x[j], colnames(dat_Norm[[i]]))] ) == F)) >= 3 } ) 
    
  }
  
  # Row to filter - IDs
  dat_Norm_Imp[[i]]$RtoF_ID <- apply(dat_Norm[[i]], 1, function(x) ifelse( sum(is.na(x)) >= 24, FALSE, TRUE) )

} 


# Missing values count
MV = forMV = list()
for ( i in 1:length(all_files)) {
  forMV[[i]] <- dat_Norm[[i]] %>%
    mutate(RtoF_ID = dat_Norm_Imp[[i]]$RtoF_ID) %>%
    filter(RtoF_ID == TRUE)
  MV[[i]] <- sum(is.na(forMV[[i]])) / ( nrow(forMV[[i]]) * ncol(forMV[[i]])) *100
}  

do.call(rbind.data.frame, MV) %>%
  rename_all(function(.) {"MissVal"}) %>%
  mutate(Exp = all_files) %>%
  write.table(., paste0("./Output_2_DataTreatment/MissingValues/", "Skyline_Library_MissingValues.txt"))

# Precursor table to export
dat = list()
for (i in 1:length(all_files)) {
  dat[[i]] <- cbind(pep[[i]][1:5], dat_Norm_Imp[[i]]) %>%
    filter(RtoF_ID == TRUE)
  write.table(dat[[i]], paste0("./Input_2_DataTreatment/", all_files[i], ".txt"), sep = "\t")
}




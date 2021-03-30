
## LIBRARIES --------------------------------------------
library(dplyr)
library(plyr)
library(tidyverse)
library(data.table)

# Store your exports in a directory called "Input_1_Spectronaut/Results"
setwd("Z:/resultats/RD/RD139_LogicielsDIA/RD139_Papier_R_210108/")

### Spectronaut ----

### 1.Fasta ----
## Load files from Spectronaut 
all_files <- list.files(path = "./Input_1_Spectronaut/Results", pattern = "Fasta_Export.csv", full.names = T)
l <- lapply(all_files, fread, sep = ",", dec = ".") 

all_files <- gsub("./Input_1_Spectronaut/Results/|_Export.csv", "", all_files)
dt <- rbindlist(l, use.names = TRUE, idcol = "Experiment")  %>%
  mutate(Experiment = ifelse(Experiment == 1, all_files[1], 
                             ifelse(Experiment == 2, all_files[2], 
                             ifelse(Experiment == 3, all_files[3], all_files[4]))) )


# To Reformate Accession 
Name_Acc <- read.csv("Input_1_Other/UPS1_Accession_Name.txt", sep="\t", header = T, stringsAsFactors = F)
# To Reformate PTM Sequence 
transition <- read.csv("Input_1_Other/transition.csv") %>% unique

# Filtering
df <- dt %>%
  as.data.frame %>%
  dplyr::select(Experiment, R.FileName, PG.ProteinAccessions, PG.Qvalue, EG.ModifiedSequence, PEP.StrippedSequence, EG.PrecursorId, EG.Qvalue, FG.Quantity) %>%
  filter(FG.Quantity > 10,
         PG.ProteinAccessions != "CATD_HUMAN_UPS") %>%
  mutate(R.FileName = gsub(".*UPS1_", "UPS1_", R.FileName) %>%
           gsub("_inj1r", "_inj1", .) %>%
           gsub("0_1fmol", "0.1fmol", .) %>%
           gsub("0_25fmol", "0.25fmol", .) %>%
           gsub("2_5fmol", "2.5fmol", .),
         PG.ProteinAccessions = ifelse(!str_detect(PG.ProteinAccessions, "HUMAN"), paste0(PG.ProteinAccessions, "_ECOLI"), gsub("_HUMAN_UPS", "", PG.ProteinAccessions)) ) %>%
  mutate(PG.ProteinAccessions = ifelse(!str_detect(PG.ProteinAccessions, "_ECOLI"), Name_Acc$Accession[match(PG.ProteinAccessions, Name_Acc$Name)], PG.ProteinAccessions),
         EG.ModifiedSequence = gsub("_", "", EG.ModifiedSequence), 
         EG.PrecursorId = str_sub(EG.PrecursorId, -1, -1)) %>%  
  dplyr::filter(EG.PrecursorId >= 2, 
                EG.PrecursorId <= 4) %>%
  dplyr::select(-c(EG.Qvalue, PG.Qvalue))          



# Change modifications 
for (y in transition$old) {
  new <- filter(transition, old == y) %>% pull(new)
  stopifnot(length(new)==1)
  df$EG.ModifiedSequence <- str_replace_all(df$EG.ModifiedSequence, y, new)         
}


# Reorder columns
df <- df %>%  
  spread(R.FileName, FG.Quantity) %>%     
  dplyr::select(Experiment, Accession = PG.ProteinAccessions, ModifiedSequence = EG.ModifiedSequence, Sequence = PEP.StrippedSequence, Precursor = EG.PrecursorId,  contains("_0.1fmol_"), contains("_0.25fmol_"), contains("_1fmol_"), contains("_2.5fmol_"),  contains("_5fmol_"), contains("_10fmol_"), contains("_25fmol_"), contains("_50fmol_")) 

# Separate tables for each window size
pep <- list()
for (i in 1:length(all_files)) {
  name <- unique(df$Experiment)[i]
  pep[[i]] <- df %>%
    filter(Experiment == name)
}


# No normalization performed for Spectronaut
dat_Norm = list()
for (i in 1:length(all_files)) {
  dat_Norm[[i]] <- pep[[i]] %>%
    dplyr::select(contains("UPS1_"))
}


# Export for BoxPlot
for (i in 1:length(all_files)) {
  write.table(dat_Norm[[i]], paste0("./Output_2_DataTreatment/Normalisation/", all_files[i], "_NormSoft.txt"))
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
  write.table(., paste0("./Output_2_DataTreatment/MissingValues/", "Spectronaut_Fasta_MissingValues.txt"))

# Precursor table to export
dat = list()
for (i in 1:length(all_files)) {
  dat[[i]] <- cbind(pep[[i]][1:5], dat_Norm_Imp[[i]]) %>%
    filter(RtoF_ID == TRUE)
  write.table(dat[[i]], paste0("./Input_2_DataTreatment/", all_files[i], ".txt"), sep = "\t")
}

rm(list = ls())






### 2.Library ----
## Load files from Spectronaut - Library
all_files <- list.files(path = "./Input_1_Spectronaut/Results", pattern = "Library_Export.csv", full.names = T)
l <- lapply(all_files, fread, sep = ",", dec = ".") 

all_files <- gsub("./Input_1_Spectronaut/Results/|_Export.csv", "", all_files)
dt <- rbindlist(l, use.names = TRUE, idcol = "Experiment")  %>%
  mutate(Experiment = ifelse(Experiment == 1, all_files[1], 
                             ifelse(Experiment == 2, all_files[2], 
                             ifelse(Experiment == 3, all_files[3], all_files[4]))) )


# To Reformate Accession 
Name_Acc <- read.csv("Input_1_Other/UPS1_Accession_Name.txt", sep="\t", header = T, stringsAsFactors = F)
# To Reformate PTM Sequence 
transition <- read_csv("Input_1_Other/transition.csv") %>% unique


# Filtering
df <- dt %>%
  as.data.frame %>%
  dplyr::select(Experiment, R.FileName, PG.ProteinAccessions, PG.Qvalue, EG.ModifiedSequence, PEP.StrippedSequence, EG.PrecursorId, EG.Qvalue, FG.Quantity) %>%
  dplyr::filter(FG.Quantity > 10) %>%
  mutate(R.FileName = gsub(".*UPS1_", "UPS1_", R.FileName) %>%
           gsub("_inj1r", "_inj1", .) %>%
           gsub("0_1fmol", "0.1fmol", .) %>%
           gsub("0_25fmol", "0.25fmol", .) %>%
           gsub("2_5fmol", "2.5fmol", .),
         PG.ProteinAccessions = ifelse(!str_detect(PG.ProteinAccessions, "HUMAN"), paste0(PG.ProteinAccessions, "_ECOLI"), gsub("_HUMAN_UPS", "", PG.ProteinAccessions)) ) %>%
  mutate(PG.ProteinAccessions = ifelse(!str_detect(PG.ProteinAccessions, "_ECOLI"), Name_Acc$Accession[match(PG.ProteinAccessions, Name_Acc$Name)], PG.ProteinAccessions),
         EG.ModifiedSequence = gsub("_", "", EG.ModifiedSequence), 
         EG.PrecursorId = str_sub(EG.PrecursorId, -1, -1)) %>%  
  dplyr::filter(EG.PrecursorId >= 2, 
                EG.PrecursorId <= 4, 
                nchar(PEP.StrippedSequence) <= 30) %>%
  dplyr::select(-c(EG.Qvalue, PG.Qvalue))



# Change modifications 
for (y in transition$old) {
  new <- filter(transition, old == y) %>% pull(new)
  stopifnot(length(new)==1)
  df$EG.ModifiedSequence <- str_replace_all(df$EG.ModifiedSequence, y, new)         
}


# Reorder columns
df <- df %>%  
  spread(R.FileName, FG.Quantity) %>%     
  dplyr::select(Experiment, Accession = PG.ProteinAccessions, ModifiedSequence = EG.ModifiedSequence, Sequence = PEP.StrippedSequence, Precursor = EG.PrecursorId,  contains("_0.1fmol_"), contains("_0.25fmol_"), contains("_1fmol_"), contains("_2.5fmol_"),  contains("_5fmol_"), contains("_10fmol_"), contains("_25fmol_"), contains("_50fmol_")) 



# Separate tables for each window size
pep <- list()
for (i in 1:length(all_files)) {
  name <- unique(df$Experiment)[i]
  pep[[i]] <- df %>%
    filter(Experiment == name)
}


# No normalization performed for Spectronaut
dat_Norm = list()
for (i in 1:length(all_files)) {
  dat_Norm[[i]] <- pep[[i]] %>%
    dplyr::select(contains("UPS1_"))
}


# Export for BoxPlot
for (i in 1:length(all_files)) {
  write.table(dat_Norm[[i]], paste0("./Output_2_DataTreatment/Normalisation/", all_files[i], "_NormSoft.txt"))
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
  write.table(., paste0("./Output_2_DataTreatment/MissingValues/", "Spectronaut_Library_MissingValues.txt"))

# Precursor table to export
dat = list()
for (i in 1:length(all_files)) {
  dat[[i]] <- cbind(pep[[i]][1:5], dat_Norm_Imp[[i]]) %>%
    filter(RtoF_ID == TRUE)
  write.table(dat[[i]], paste0("./Input_2_DataTreatment/", all_files[i], ".txt"), sep = "\t")
}

rm(list = ls())



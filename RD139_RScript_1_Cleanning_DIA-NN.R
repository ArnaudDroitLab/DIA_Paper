
## LIBRARIES --------------------------------------------
library(dplyr)
library(plyr)
library(tidyverse)
library(data.table)


# Store your exports in a directory called "Input_1_DIANN/Results"
setwd("Z:/resultats/RD/RD139_LogicielsDIA/RD139_Papier_R_210108/")


## DIA-NN ----
### 1.Fasta ----

# Load files for DIA-NN - Fasta
all_files <- list.files(path = "./Input_1_DIANN/Results/", pattern = "_Fasta_Export.tsv", full.names = T)
l <- lapply(all_files, fread, sep = "\t") 

all_files <- gsub("./Input_1_DIANN/Results/|_Export.tsv", "", all_files)
dt <- rbindlist(l, use.names = TRUE, idcol = "Experiment")  %>%
  mutate(Experiment = ifelse(Experiment == 1, all_files[1], 
                          ifelse(Experiment == 2, all_files[2], 
                          ifelse(Experiment == 3, all_files[3], all_files[4]))) )


# Filtering
df <- dt %>%  
  as.data.frame %>%
  mutate(Protein.Group = gsub("\\;.*", "", Protein.Group)) %>%
  filter(Protein.Group != "P08311ups|CATG_HUMAN_UPS") %>%
  dplyr::select(Experiment, File.Name, Protein.Group, Protein.Names, Modified.Sequence, Stripped.Sequence, Precursor.Charge, Precursor.Normalised, Lib.Q.Value, Q.Value, PG.Q.Value) %>%
  mutate(File.Name = gsub("\\\\", "", File.Name) %>%
           gsub("_inj1r", "_inj1", .) %>%
           gsub("C:UsersproteoDesktopRD139_Feb2021mzMLRD139_Mixed_|C:UsersproteoDesktopRD139_Feb2021mzMLRD139_Narrow_|C:UsersproteoDesktopRD139_Feb2021mzMLRD139_Overlap_|C:UsersproteoDesktopRD139_Feb2021mzMLRD139_Wide_|.mzML", "", .) %>%
           gsub("0_1fmol", "0.1fmol", .) %>%
           gsub("0_25fmol", "0.25fmol", .) %>%
           gsub("2_5fmol", "2.5fmol", .))  %>%
  filter(Precursor.Charge != 1,
        # Lib.Q.Value <= 0.01,    already filtered 
        # Q.Value <= 0.01,        already filtered 
         PG.Q.Value <= 0.01) %>%
  dplyr::rename(ModifiedSequence = Modified.Sequence, Accession = Protein.Group, Sequence = Stripped.Sequence, Precursor = Precursor.Charge) %>%
  mutate(specie = ifelse(str_detect(Protein.Names, "ECOLI"), "ECOLI", "HUMAN"),
         id = if_else(specie == "HUMAN", 
                      str_extract(Accession, "[A-Z][0-9]*"), 
                      Accession),
         Accession = paste(id, specie, sep = "_")) %>% 
  dplyr::select(-c(Protein.Names, Lib.Q.Value, Q.Value, PG.Q.Value, specie, id))


# NA replacement instead of 0
df$Precursor.Normalised[df$Precursor.Normalised == 0 ] <- NA

# Reorder columns
df <- df %>%
  spread(File.Name, Precursor.Normalised) %>%  
  dplyr::select(Experiment : Precursor,
         contains("_0.1fmol_"), contains("_0.25fmol_"), contains("_1fmol_"), contains("_2.5fmol_"),
         contains("_5fmol_"), contains("_10fmol_"), contains("_25fmol_"), contains("_50fmol_"))

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


# No normalization performed for DIA-NN
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
  dat_Norm_Imp[[i]]$RtoF_ID <- apply(dat_Norm[[i]], 1, function(x) ifelse( sum(is.na(x)) >= 24, FALSE,TRUE))

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
write.table(., paste0("./Output_2_DataTreatment/MissingValues/", "DIA-NN_Fasta_MissingValues.txt"))

# Precursor table to export
dat = list()
for (i in 1:length(all_files)) {
  dat[[i]] <- cbind(pep[[i]][1:5], dat_Norm_Imp[[i]]) %>%
    filter(RtoF_ID == TRUE)
  write.table(dat[[i]], paste0("./Input_2_DataTreatment/", all_files[i], ".txt"), sep = "\t")
}

rm(list = ls())


### 2.Library from Skyline ----

# Load files from DIA-NN - Library
all_files <- list.files(path = "./Input_1_DIANN/Results/", pattern = "_Library_Export.tsv", full.names = T)
l <- lapply(all_files, fread, sep = "\t") 

all_files <- gsub("./Input_1_DIANN/Results/|_Export.tsv", "", all_files)
dt <- rbindlist(l, use.names = TRUE, idcol = "Experiment")  %>%
  as.data.frame %>%
  mutate(Experiment = ifelse(Experiment == 1, all_files[1], 
                          ifelse(Experiment == 2, all_files[2], 
                          ifelse(Experiment == 3, all_files[3], all_files[4]))) )


# Filtering
df <- dt %>%  
  as.data.frame %>%
  dplyr::select(Experiment, File.Name, Protein.Group, Modified.Sequence, Stripped.Sequence, 
                Precursor.Charge, Precursor.Normalised, Q.Value, PG.Q.Value) %>%
  mutate(File.Name = gsub("\\\\", "", File.Name) %>%
           gsub("_inj1r", "_inj1", .) %>%
           gsub("C:UsersproteoDesktopRD139_Feb2021mzMLRD139_Mixed_|C:UsersproteoDesktopRD139_Feb2021mzMLRD139_Narrow_|C:UsersproteoDesktopRD139_Feb2021mzMLRD139_Overlap_|C:UsersproteoDesktopRD139_Feb2021mzMLRD139_Wide_|.mzML", "", .) %>%
           gsub("0_1fmol", "0.1fmol", .) %>%
           gsub("0_25fmol", "0.25fmol", .) %>%
           gsub("2_5fmol", "2.5fmol", .))  %>%
  filter(Protein.Group != "Biognosys standards",
       # Q.Value <= 0.01,        already filtered 
         PG.Q.Value <= 0.01) %>%
  dplyr::rename(ModifiedSequence = Modified.Sequence, Accession = Protein.Group, Sequence = Stripped.Sequence, Precursor = Precursor.Charge) %>%
  mutate(specie = ifelse(str_detect(Accession, "ECOLI"), "ECOLI", "HUMAN"),
         id = if_else(specie == "HUMAN", 
                      str_extract(Accession, "[A-Z][0-9]*"), 
                      str_extract(Accession, "\\|[^\\|]*\\|") %>% str_replace_all("\\|", "")),
         Accession = paste(id, specie, sep = "_")) %>% 
  dplyr::select(-c(Q.Value, PG.Q.Value, specie, id))


# NA replacement instead of 0
df$Precursor.Normalised[df$Precursor.Normalised == 0 ] <- NA

# Reorder columns
df <- df %>% 
  spread(File.Name, Precursor.Normalised) %>%  
  select(Experiment : Precursor,
         contains("_0.1fmol_"), contains("_0.25fmol_"), contains("_1fmol_"), contains("_2.5fmol_"),
         contains("_5fmol_"), contains("_10fmol_"), contains("_25fmol_"), contains("_50fmol_")) 

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

# No nomalization for DIA-NN
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
  write.table(., paste0("./Output_2_DataTreatment/MissingValues/", "DIA-NN_Library_MissingValues.txt"))


# Precursor table to export
dat = list()
for (i in 1:length(all_files)) {
  dat[[i]] <- cbind(pep[[i]][1:5], dat_Norm_Imp[[i]])  %>%
    filter(RtoF_ID == TRUE)
  write.table(dat[[i]], paste0("./Input_2_DataTreatment/", all_files[i], ".txt"), sep = "\t")
}

rm(list = ls())

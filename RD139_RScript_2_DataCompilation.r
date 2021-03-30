## LIBRARIES --------------------------------------------
library(dplyr)
library(plyr)
library(tidyverse)
library(data.table)
library(pROC)
library(MLmetrics) 


## WORKING DIRECTORY ----
setwd("Z:/resultats/RD/RD139_LogicielsDIA/RD139_Papier_R_210108/")


## LOAD FILES ----
all_files <- list.files(path = "./Input_2_DataTreatment/", pattern = ".txt", full.names = T)
l <- lapply(all_files, fread, sep = "\t") 

dt <- rbindlist(l, use.names = T) %>%
  as.data.frame %>%
  dplyr::select(-V1) %>%
  dplyr::rename(Experience = Experiment)

Exp <- unique(dt$Experience)
dat = list()

for(i in 1:length(Exp)) {
  dat[[i]] <- dt %>%
    filter(Experience == Exp[i])
}


## NUMBER OF IDs / QUANT ----

# Empty tables
names = c("Experience", "0.1fmol", "0.25fmol", "1fmol", "2.5fmol", "5fmol", "10fmol", "25fmol", "50fmol")

ID_prot_Eco <- ID_prot_UPS <- ID_pep_Eco <- ID_pep_UPS <- matrix(ncol = 9, nrow = length(Exp)) %>%
  as.data.frame %>%
  rename_all(function(.) {names}) %>%
  mutate(Experience = Exp)

# Fill df
x = x_UPS = x_Eco = list()

for (i in 1:length(Exp)) {
  
  x[[i]] <- dat[[i]] %>%
    mutate(Specie = ifelse(str_detect(Accession, "HUMAN"), "HUMAN", "Ecoli")) %>%
    dplyr::select(Accession, Sequence, Specie, starts_with("RtoF_"))   
  x_UPS[[i]] <- x[[i]] %>%
    filter(Specie == "HUMAN")
  x_Eco[[i]] <- x[[i]] %>%
    filter(Specie == "Ecoli")
  
  for(j in 1:8) {
    
    # Peptides
    ID_pep_UPS[i, j+1] <- length(unique(x_UPS[[i]]$Sequence[which(x_UPS[[i]][,j+3] == TRUE )]))   
    ID_pep_Eco[i, j+1] <- length(unique(x_Eco[[i]]$Sequence[which(x_Eco[[i]][,j+3] == TRUE )]))  
    
    # Proteins
    ID_prot_UPS[i, j+1] <- length(unique(x_UPS[[i]]$Accession[which(x_UPS[[i]][,j+3] == TRUE )]))
    ID_prot_Eco[i, j+1] <- length(unique(x_Eco[[i]]$Accession[which(x_Eco[[i]][,j+3] == TRUE )]))  
    
  }
  
  # Add number of ID values (seen once on the 24 samples)
  # Eco 
  ID_pep_Eco$ID[i] <- length(unique(x_Eco[[i]]$Sequence))
  ID_prot_Eco$ID[i] <- length(unique(x_Eco[[i]]$Accession))

  # UPS1
  ID_pep_UPS$ID[i] <- length(unique(x_UPS[[i]]$Sequence))
  ID_prot_UPS$ID[i] <- length(unique(x_UPS[[i]]$Accession))
  
}


# Add mean of quantification for Eco 
ID_pep_Eco$Quant <- apply(ID_pep_Eco[,2:9], 1, mean)
ID_pep_Eco$sd <- apply(ID_pep_Eco[,2:9], 1, sd)
ID_prot_Eco$Quant <- apply(ID_prot_Eco[,2:9], 1, mean)
ID_prot_Eco$sd <- apply(ID_prot_Eco[,2:9], 1, sd)


# Export
# Identifications Eco
write.table(ID_pep_Eco, "Output_2_DataTreatment/SuppTable_IDs_Quant_pep_Eco.txt")
write.table(ID_prot_Eco, "Output_2_DataTreatment/SuppTable_IDs_Quant_prot_Eco.txt")

# Identifications UPS1
write.table(ID_pep_UPS, "Output_2_DataTreatment/SuppTable_IDs_Quant_pep_UPS.txt")
write.table(ID_prot_UPS, "Output_2_DataTreatment/SuppTable_IDs_Quant_prot_UPS.txt")



## REPRODUCIBILITY - CV ----

tmp = list()
tmp_pep = list()
tmp_prot = list()
Seq = Acc = list()

Conc = c("UPS1_0.1fmol_", "UPS1_0.25fmol_", "UPS1_1fmol_", "UPS1_2.5fmol_", 
         "UPS1_5fmol_", "UPS1_10fmol_", "UPS1_25fmol_", "UPS1_50fmol_" )
RtF = c("RtoF_0.1fmol", "RtoF_0.25fmol", "RtoF_1fmol", "RtoF_2.5fmol",
        "RtoF_5fmol", "RtoF_10fmol", "RtoF_25fmol", "RtoF_50fmol")

CV_pep <- as.data.frame(matrix(ncol = 9, nrow = 0)) %>%
  dplyr::rename_all(function(.) {c("Accession", "Sequence", "Rep1", "Rep2", "Rep3", "Specie", "CV", "Concentration", "Experience")} )

CV_prot <- as.data.frame(matrix(ncol = 8, nrow = 0)) %>%
  dplyr::rename_all(function(.) {c("Accession", "Rep1", "Rep2", "Rep3", "Specie", "Experience", "CV", "Concentration")} )



for ( i in 1:length(Exp)) {

  # For the 8 concentrations
  tmp[[1]] = tmp[[2]] = tmp[[3]] = tmp[[4]] = tmp[[5]] = tmp[[6]] = tmp[[7]] = tmp[[8]] = dat[[i]]
  
  
  for (x in 1:8) {
    
    # Keep only one concentration column (intensities) and precursors detected 3 times over the 3 replicates in this group
    tmp[[x]] <- tmp[[x]] %>%
      dplyr::select(Accession, Sequence, contains(Conc[x]), contains(RtF[x]))  %>%
      dplyr::filter((.[,6]) == TRUE)
     
    Acc[[x]] <- tmp[[x]]$Accession
    Seq[[x]] <- tmp[[x]]$Sequence
    
 
    # Aggregate by peptide stripped sequence
    tmp_pep[[x]] <- tmp[[x]] %>%
      dplyr::select(starts_with("UPS1_")) %>%
      aggregate(x = .,  by = list(Acc[[x]], Seq[[x]]),  FUN = "sum", na.rm = T) %>%
      rename_all(function(.) {c("Accession", "Sequence", "Rep1", "Rep2", "Rep3")}) %>%
      mutate(Specie = ifelse(str_detect(Accession, "ECOLI"), "Ecoli", "HUMAN"),
             CV = apply(.[,3:5], 1, function(z) { sd(z) / mean(z)*100 }) ,
             Concentration = gsub("RtoF_", "", RtF[x]), 
             Experience = Exp[i])
   
    # Aggreagte by accession number
    tmp_prot[[x]] <- tmp[[x]] %>%
      dplyr::select(-c(Accession, Sequence, starts_with("RtoF_"))) %>%
      aggregate(x = .,  by = list(Acc[[x]]),  FUN = "sum", na.rm = T) %>%
      mutate(Specie = ifelse(str_detect(Group.1, "ECOLI"), "Ecoli", "HUMAN"),
             Experience = Exp[i]) %>%
      dplyr::rename_all(function(.) {c( "Accession", "Rep1", "Rep2", "Rep3", "Specie", "Experience")}) %>%
      na.omit %>%
      mutate(CV = apply(.[,2:4], 1, function(z) { sd(z) / mean(z) * 100 }) , 
             Concentration = gsub("RtoF_", "", RtF[x])) 
    
    
  } # End of loop - concentrations
  
  # Peptide
  CV_x <- do.call(rbind.data.frame, tmp_pep) 
  CV_pep <- rbind(CV_pep, CV_x)
  # Protein
  CV_y <- do.call(rbind.data.frame, tmp_prot) 
  CV_prot <- rbind(CV_prot, CV_y)

} # End of loop - Experience


# Export 
CV_pep %>% 
  dplyr::select(Experience, Specie, CV, Concentration) %>%
  write.table(., "Output_2_DataTreatment/SuppTable_CV_pep.txt")
CV_prot %>% 
  dplyr::select(Experience, Specie, CV, Concentration) %>%
  write.table(., "Output_2_DataTreatment/SuppTable_CV_prot.txt")






## LINEARITY - R2 ----

tmp = tmp2 = list()
df_lin = reg = list()

# Regression function
regression <- function(df_lin) {
  reg_fun <- lm(formula = df_lin$Int ~ df_lin$Concentration) 
  slope <- round(coef(reg_fun)[2], 3)  
  intercept <- round(coef(reg_fun)[1], 3) 
  R2 <- round(as.numeric(summary(reg_fun)[8]), 3)
  R2.Adj <- round(as.numeric(summary(reg_fun)[9]), 3)
  c(slope, intercept, R2, R2.Adj)
}


# For linearities - from 1fmol to 50fmol
for ( i in 1:length(Exp)) {
  
  # Keep only the 6 highest UPS1 concentrations (from 1fmol to 50fmol)
  tmp[[i]] <- dat[[i]] %>%
    dplyr::select( contains( "UPS1_") ) %>%
    dplyr::select( -contains("_0.1fmol_"), -contains("_0.25fmol_"))
  
  # Sum of precursor intensities by accession, if 3 quant values are detectable in the 3 replicates in the 50fmol condition
  tmp[[i]] <- aggregate(x = tmp[[i]][which(dat[[i]]$RtoF_50fmol == T),], 
                        by = list(as.character(dat[[i]]$Accession[which(dat[[i]]$RtoF_50fmol == T)])), 
                        FUN = "sum", na.rm = T) %>%
    gather(Sample, Int, -Group.1) %>%
    mutate(Specie = ifelse(str_detect(Group.1, "HUMAN"), "UPS1", "Ecoli"),
           Concentration = gsub("UPS1_|_inj.*", "", Sample) %>%
             factor(., levels = c("1fmol", "2.5fmol", "5fmol", "10fmol", "25fmol", "50fmol"))) 
  
  # Mean of the 3 replicates and log2 transform
  tmp2[[i]] <- aggregate(data = tmp[[i]], x = tmp[[i]]$Int, 
                          by = list(tmp[[i]]$Concentration , tmp[[i]]$Group.1), 
                          FUN = mean, na.rm = TRUE) %>%
    dplyr::select(Concentration = Group.1, Accession = Group.2, Int = x) %>%
    mutate(Experience = Exp[i],
           Specie = ifelse(str_detect(Accession, "HUMAN"), "UPS1", "Ecoli"),
           Int = log(Int, 2)) %>%
    spread(Concentration, Int)
  
  # Keep uPS1 and log10(concentration)
  df_lin[[i]] <-  tmp2[[i]] %>%
    dplyr::filter(Specie == "UPS1") %>%
    gather(Sample, Int, -Accession, -Experience, -Specie) %>%
    mutate(Concentration = gsub("fmol", "", Sample) %>% as.numeric %>% log10)
  
  # Regression data
  reg[[i]] <- ddply(df_lin[[i]], "Accession", regression) 
  colnames(reg[[i]]) <- c("Accession", "Slope", "Intercept", "R2", "R2.Adj")
  reg[[i]]$Experience <- Exp[i]
  
} 

# Export Linearity data
do.call(rbind.data.frame, reg)  %>%
  write.table(.,  "Output_2_DataTreatment/SuppTable_Linearities.txt")



## ACCURACY - MAPE, SENSITIVITY - AUC ----
# Ratios to do 
tab <- read.table("Input_1_Other/Tab_ratio.csv", 
                  sep = ";", header = T, stringsAsFactors = F)

# Create empty tables to compile data
# AUC
AUC_all <- matrix(nrow = nrow(tab), ncol = length(Exp)+1) %>%
  as.data.frame %>%
  dplyr::rename_all(function(.) {c("Ratio", Exp)} ) %>%
  mutate(Ratio = tab$Ratio)

# Sensitivity - FDP
names = c(tab$Ratio, "Experience", "Type")
X <- data.frame( matrix(ncol = length(names), nrow = length(Exp)*4) ) %>% 
  dplyr::rename_all(function(.) {names} ) %>%
  mutate(Experience = rep(Exp, each = 4),
         Type = rep(c("True Positive",  "True Negative", "False Positive", "False Negative"), 
                    time = length(Exp))) 

Seq_TP <- seq(1, nrow(X), by=4)
Seq_TN <- seq(2, nrow(X), by=4)
Seq_FP <- seq(3, nrow(X), by=4)
Seq_FN <- seq(4, nrow(X), by=4)

# MAPE
MAPE_UPS_all <- matrix(nrow = 0, ncol = 3) %>% 
  as.data.frame %>%
  dplyr::rename_all(function(.) {c("Experience", "log2Ratio", "MAPE")} )

# Intermediate tables
tmp  =  list()
Acc = list()
ROC = list()
AUC = list()
forMAPE = list()
MAPE_UPS = MAPE_UPS2 = list()
export = list()


for (i in 1:length(Exp)) {
  
  # For the 28 pairwise comparisons
  tmp[[1]] = tmp[[2]] = tmp[[3]] = tmp[[4]] = tmp[[5]] = tmp[[6]] = tmp[[7]] = tmp[[8]] = tmp[[9]] = tmp[[10]] = tmp[[11]] = tmp[[12]] = tmp[[13]] = tmp[[14]] = tmp[[15]] = tmp[[16]] = tmp[[17]] = tmp[[18]] = tmp[[19]] = tmp[[20]] = tmp[[21]] = tmp[[22]] = tmp[[23]] = tmp[[24]] = tmp[[25]] = tmp[[26]] = tmp[[27]] = tmp[[28]] = dat[[i]] 
  
  for(x in 1:nrow(tab)) {  
    
    # Keep only pairwise columns intensities and precursors detected 3 times over the 3 replicates in one of the group of UPS1 concentrations.
    tmp[[x]] <- tmp[[x]] %>%
      dplyr::select(Accession, contains(as.character(tab$g1[x])), contains(as.character(tab$g2[x])), 
                    as.character(tab$RtF_g1[x]), as.character(tab$RtF_g2[x]) ) %>%
      dplyr::filter((.[,8]) == TRUE | (.[,9]) == TRUE ) %>%
      dplyr::select(-starts_with("RtoF_"))
    
    Acc[[x]] <- tmp[[x]]$Accession
    
    # Sum of precursor intensities by accession number to obtain protein level quantification
    # Ratio, and z-score calculation
    tmp[[x]] <- tmp[[x]] %>%
      dplyr::select(-Accession) %>%
      aggregate(x = .,  by = list(Acc[[x]]),  FUN = "sum", na.rm = T) %>%
      mutate(Exp = Exp[[i]],
             Specie = ifelse(str_detect(Group.1, "HUMAN"), "HUMAN", "ECOLI"),
             mean_g1 = apply(.[2:4], 1 , mean), 
             mean_g2 = apply(.[5:7], 1 , mean),
             Ratio = mean_g1 / mean_g2,
             log2Ratio = log(Ratio, 2), 
             ExpRatio = ifelse(Specie == "HUMAN", tab$ExpRatio[x], 1) %>% as.numeric,
             log2ExpRatio = log(ExpRatio, 2),
             zscore = as.numeric(scale(log2Ratio, center = TRUE, scale = TRUE)))
    
    # Welch t-test
    Pval <- numeric(nrow(tmp[[x]]))
    for(a in 1:nrow(tmp[[x]])) {
      Pval[a] <- t.test(tmp[[x]][a,2:4], tmp[[x]][a,5:7])$p.value 
    }
    
    # Calcul of adjusted p-value (Benjamini Hochberg)
    # Determination of True Positive, True Negative, False Positive and False Negative proteins
    tmp[[x]] <- tmp[[x]] %>%
      mutate(pval = Pval,
             qval = p.adjust(pval , method = "BH"),
             log10pval = -log(pval, 10), 
             Reg = ifelse( ((abs(zscore) > 1.96) & (qval < 0.05)), "Reg", "NoReg"),
             TP_FN = ifelse( (Reg == "Reg" & Specie == "HUMAN"), "True Positive",
                     ifelse( (Reg == "Reg" & Specie == "ECOLI"), "False Positive", 
                     ifelse( (Reg == "NoReg" & Specie == "HUMAN"), "False Negative", "True Negative" ))) )
    
    # ROC 
    ifelse(length(unique(tmp[[x]]$Specie)) == 2, 
           ROC[[x]] <- tmp[[x]] %>% 
             roc(Specie, qval, levels = c("ECOLI", "HUMAN"), direction = ">"),
           ROC[[x]] <-  NA)
    
    # AUC  
    ifelse(class(ROC[[x]]) %in% c("logical"),
           AUC[[x]] <- NA,
           AUC[[x]] <- auc(ROC[[x]])  %>% as.numeric)
    
    # TP - TN - FP - FN
    X[Seq_TP[i], x] <- length( which(tmp[[x]]$TP_FN == "True Positive") )
    X[Seq_TN[i], x] <- length( which(tmp[[x]]$TP_FN == "True Negative") )
    X[Seq_FP[i], x] <- length( which(tmp[[x]]$TP_FN == "False Positive") )
    X[Seq_FN[i], x] <- length( which(tmp[[x]]$TP_FN == "False Negative") )
    
    
    # MAPE UPS
    ifelse(length(unique(tmp[[x]]$Specie)) == 2, 
           forMAPE[[x]] <- tmp[[x]] %>% 
             mutate(RATIO = tab$Ratio[x]) %>% 
             filter(Specie == "HUMAN") %>%
             dplyr::select(Group.1, RATIO, log2Ratio, log2ExpRatio), 
           forMAPE[[x]] <- matrix(ncol = 4, nrow = 1) %>%
             as.data.frame %>%
             dplyr::rename_all(function(.) {c("Group.1", "RATIO", "log2Ratio", "log2ExpRatio")}) %>%
             mutate(RATIO = tab$Ratio[x]))
    
    forMAPE1 <- do.call(rbind.data.frame, forMAPE) %>%
      mutate(RATIO = factor(RATIO, levels = tab$Ratio)) %>%
      spread(Group.1, log2Ratio) 
    
    MAPE_UPS[[x]] <- forMAPE1[x,] %>%
      gather(Protein, log2RExp, -RATIO, -log2ExpRatio) %>%
      dplyr::select(RATIO = RATIO, log2Rth = log2ExpRatio, log2RExp = log2RExp) %>% 
      na.omit
    
    ifelse(nrow(MAPE_UPS[[x]]) == 0, 
           MAPE_UPS2[[x]] <- NA, 
           MAPE_UPS2[[x]] <- MAPE(MAPE_UPS[[x]]$log2RExp, MAPE_UPS[[x]]$log2Rth))
    
  } # End of loop - Ratios
  
  
  # Compile all AUC from all experience in a single table
  z <- AUC %>% unlist %>% as.vector
  AUC_all[, i+1] <- z
  
  # Compile all MAPE from all experience in a single table
  MAPE_UPS3 <- matrix(ncol = 3, nrow = length(MAPE_UPS2)) %>%
    as.data.frame %>%
    dplyr::rename_all( function(.) {c("Experience", "Ratio", "MAPE")} )  %>%
    mutate(Experience = Exp[[i]],
           Ratio =  factor(tab$Ratio, levels = as.character(tab$Ratio) ),
           MAPE = MAPE_UPS2 %>% unlist) %>%
    mutate(MAPE = MAPE * 100)
  MAPE_UPS_all <- rbind(MAPE_UPS_all, MAPE_UPS3)
  
} # End of loop - Experience


# Sensitivity - FDP table
X <- X %>%
  gather(Ratio, Count, -c(Experience, Type)) %>%
  mutate(Type = gsub(" ", "", Type),
         Ratio = factor(Ratio, levels = as.character(tab$Ratio))) %>%
  spread(Type, Count) %>%
  mutate(Sensitivity = TruePositive / (TruePositive + FalseNegative) * 100,
         Specificity = (1 - (TrueNegative / (FalsePositive + TrueNegative))) * 100, 
         FDP = FalsePositive / (TruePositive + FalsePositive) * 100)

# Export
# AUC
write.table(AUC_all, "Output_2_DataTreatment/SuppTable_AUC.txt")

# MAPE UPS 
MAPE_UPS_all %>%
  unique %>%
  spread(Ratio, MAPE) %>%
  write.table(., "Output_2_DataTreatment/SuppTable_MAPE_UPS_all.txt")

# Specificity - FDP
write.table(X, "Output_2_DataTreatment/SuppTable_TP-TN-FP-FN_qval005_zscore.txt")


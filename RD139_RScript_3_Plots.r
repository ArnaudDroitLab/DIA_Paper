## LIBRARIES --------------------------------------------
library(ggplot2)
library(gridExtra)
library(ggpubr)



## WORKING DIRECTORY --
setwd("Z:/resultats/RD/RD139_LogicielsDIA/RD139_Papier_R_210108/")

## COLORS -----------------------------------------------
Col <- data.frame(Software = c("DIA-NN","DIA-Umpire", "OpenSWATH", "Skyline", "ScaffoldDIA", "Spectronaut" ),
                  Color = c("#1b6bb5", "#9e107b", "#5e5b56",  "#17962a", "#ca431d", "#ff9d00" ))

set_color <- function(df) {
  tibble(Software = df %>% pull(Software) %>% unique) %>%
    left_join(Col, by = "Software") %>%
    mutate(Software = factor(Software, levels = levels(df$Software))) %>%
    arrange(Software)
}

## ORDER ------------------------------------------------
Lev_Acq = c("Narrow", "Overlap", "Mixed", "Wide")
Lev_Conc = c("0.1fmol", "0.25fmol", "1fmol", "2.5fmol", "5fmol", "10fmol", "25fmol", "50fmol")


## 1 - QUANTIFICATIONS - PEPTIDES ----------------------------
x_Eco <- read.csv("./Output_2_DataTreatment/SuppTable_IDs_Quant_pep_Eco.txt", header = T, sep = " ", stringsAsFactors = F) %>%
  dplyr::select(Experience, Quant, sd) %>%
  separate(Experience, c("Software","Acquisition", "Extraction"), sep = "_") %>%
  mutate(Type = "Peptides ", 
         Software = factor(Software), 
         Acquisition = factor(Acquisition, levels = Lev_Acq))

x_UPS <-  read.csv("./Output_2_DataTreatment/SuppTable_IDs_Quant_pep_UPS.txt", header = T, sep =" ") %>%
  dplyr::select(-ID) %>%
  gather(Concentration, Count, -Experience) %>%
  separate(Experience, c("Software","Acquisition", "Extraction"), sep = "_") %>%
  mutate(Concentration = gsub("X", "", Concentration) %>% 
           gsub("fmol", " fmol", .) %>% 
           factor(., levels = unique(.)), 
         Type = "Peptides", 
         Software = factor(Software), 
         Acquisition = factor(Acquisition, levels = Lev_Acq)) 

# Eco fasta
df <- x_Eco %>% 
  filter(Extraction == "Fasta")
colors <- set_color(df) %>% mutate(Color = as.character(Color))
p1 <- ggplot(df, aes(Type, Quant, color = Software, group = Software)) +
  geom_point(size = 0.8) +
  geom_errorbar(aes(ymin = Quant - sd, ymax = Quant + sd), width = 0.2, alpha = 0.8, size = 0.4) +
  facet_grid(rows = vars(Type), cols = vars(Acquisition)) +
  theme_minimal() +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Fasta\nE.coli", limits = c(3500, 19000)) +
  scale_color_manual(values = colors$Color) +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1, vjust = 0.25, color = "white"),
        axis.text.y = element_text(size = 6.5), 
        legend.position = "none", 
        legend.title = element_blank(),
        strip.text.y = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold")) 

# UPS fasta
df <- x_UPS %>%
  filter(Extraction == "Fasta")
colors <- set_color(df) %>% mutate(Color = as.character(Color))
p2 <- ggplot(df, aes(Concentration, Count, color = Software, group = Software))  +
  geom_point(size = 0.8) +
  geom_line(size = 0.5) +
  facet_grid(rows = vars(Type), cols = vars(Acquisition)) +
  theme_minimal() +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Fasta\nUPS1", limits = c(0, 710)) +
  scale_color_manual(values = colors$Color) +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1, vjust = 0.25),
        axis.text.y = element_text(size = 6.5),
        legend.position = "none", 
        legend.title = element_blank(),
        strip.text.y = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold")) 

# Eco library
df <- x_Eco %>% 
  filter(Extraction == "Library") 
colors <- set_color(df) %>% mutate(Color = as.character(Color))
p3 <-ggplot(df, aes(Type, Quant, color = Software, group = Software)) +
  geom_point(size = 0.8) +
  geom_errorbar(aes(ymin = Quant - sd, ymax = Quant + sd), width = 0.2, alpha = 0.8, size = 0.4) +
  facet_grid(rows = vars(Type), cols = vars(Acquisition)) +
  theme_minimal() +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Library\nE.coli", limits = c(3500, 19000)) +
  scale_color_manual(values = colors$Color) +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1, vjust = 0.25, color = "white"),
        axis.text.y = element_text(size = 6.5), 
        legend.position = "none", 
        legend.title = element_blank(),
        strip.text.y = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold")) 

# UPS library
df <- x_UPS %>%
  filter(Extraction == "Library") 
colors <- set_color(df) %>% mutate(Color = as.character(Color))
p4 <- ggplot(df, aes(Concentration, Count, color = Software, group = Software))  +
  geom_point(size = 0.8) +
  geom_line(size = 0.5) +
  facet_grid(rows = vars(Type), cols = vars(Acquisition)) +
  theme_minimal() +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Library\nUPS1", limits = c(0, 710)) +
  scale_color_manual(values = colors$Color) +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1, vjust = 0.25),
        axis.text.y = element_text(size = 6.5),
        legend.position = "none", 
        legend.title = element_blank(),
        strip.text.y = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold"))

# Legend 
df <- x_Eco
colors <- set_color(df) %>% mutate(Color = as.character(Color))
leg <- ggplot(df, aes(Type, Quant, color = Software, group = Software)) +
  geom_point() +
  geom_line(size = 0.5) +
  scale_color_manual(values = colors$Color) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "bottom", 
        legend.text = element_text(size = 8))
leg <- get_legend(leg)  
leg <- as_ggplot(leg)



pdf("Output_2_DataTreatment/Figures/1_Quant_Peptides.pdf", height = 5, width = 8)
grid.arrange(grobs = list(p1, p2, p3, p4, leg),
             ncol = 2, nrow = 3, 
             widths = c(2, 5),  
             heights = c(3, 3, 1),
             layout_matrix = rbind(c(1 ,2), 
                                   c(3, 4),
                                   c(5, 5)), 
             top = "PEPTIDES\n")
dev.off()


## 1 - QUANTIFICATIONS - PROTEINS ----------------------------
x_Eco <- read.csv("./Output_2_DataTreatment/SuppTable_IDs_Quant_prot_Eco.txt", header = T, sep = " ", stringsAsFactors = F) %>%
  dplyr::select(Experience, Quant, sd) %>%
  separate(Experience, c("Software","Acquisition", "Extraction"), sep = "_") %>%
  mutate(Type = "Proteins ", 
         Software = factor(Software), 
         Acquisition = factor(Acquisition, levels = Lev_Acq))

x_UPS <-  read.csv("./Output_2_DataTreatment/SuppTable_IDs_Quant_prot_UPS.txt", header = T, sep =" ") %>%
  dplyr::select(-ID) %>%
  gather(Concentration, Count, -Experience) %>%
  separate(Experience, c("Software","Acquisition", "Extraction"), sep = "_") %>%
  mutate(Concentration = gsub("X", "", Concentration) %>% 
           gsub("fmol", " fmol", .) %>% 
           factor(., levels = unique(.)), 
         Type = "Peptides", 
         Software = factor(Software), 
         Acquisition = factor(Acquisition, levels = Lev_Acq)) 

# Eco fasta
df <- x_Eco %>% 
  filter(Extraction == "Fasta")
colors <- set_color(df) %>% mutate(Color = as.character(Color))
p1 <- ggplot(df, aes(Type, Quant, color = Software, group = Software)) +
  geom_point(size = 0.8) +
  geom_errorbar(aes(ymin = Quant - sd, ymax = Quant + sd), width = 0.2, alpha = 0.8, size = 0.4) +
  facet_grid(rows = vars(Type), cols = vars(Acquisition)) +
  theme_minimal() +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Fasta\nE.coli", limits = c(700, 2200)) +
  scale_color_manual(values = colors$Color) +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1, vjust = 0.25, color = "white"),
        axis.text.y = element_text(size = 6.5), 
        legend.position = "none", 
        legend.title = element_blank(),
        strip.text.y = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold")) 

# UPS fasta
df <- x_UPS %>%
  filter(Extraction == "Fasta")
colors <- set_color(df) %>% mutate(Color = as.character(Color))
p2 <- ggplot(df, aes(Concentration, Count, color = Software, group = Software))  +
  geom_point(size = 0.8) +
  geom_line(size = 0.5) +
  facet_grid(rows = vars(Type), cols = vars(Acquisition)) +
  theme_minimal() +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Fasta\nUPS1", limits = c(0, 50)) +
  scale_color_manual(values = colors$Color) +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1, vjust = 0.25),
        axis.text.y = element_text(size = 6.5),
        legend.position = "none", 
        legend.title = element_blank(),
        strip.text.y = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold")) 

# Eco library
df <- x_Eco %>% 
  filter(Extraction == "Library") 
colors <- set_color(df) %>% mutate(Color = as.character(Color))
p3 <- ggplot(df, aes(Type, Quant, color = Software, group = Software)) +
  geom_point(size = 0.8) +
  geom_errorbar(aes(ymin = Quant - sd, ymax = Quant + sd), width = 0.2, alpha = 0.8, size = 0.4) +
  facet_grid(rows = vars(Type), cols = vars(Acquisition)) +
  theme_minimal() +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Library\nE.coli", limits = c(700, 2200)) +
  scale_color_manual(values = colors$Color) +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1, vjust = 0.25, color = "white"),
        axis.text.y = element_text(size = 6.5),
        legend.position = "none", 
        legend.title = element_blank(),
        strip.text.y = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold")) 

# UPS library
df <- x_UPS %>%
  filter(Extraction == "Library") 
colors <- set_color(df) %>% mutate(Color = as.character(Color))
p4 <- ggplot(df, aes(Concentration, Count, color = Software, group = Software))  +
  geom_point(size = 0.8) +
  geom_line(size = 0.5) +
  facet_grid(rows = vars(Type), cols = vars(Acquisition)) +
  theme_minimal() +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Library\nUPS1", limits = c(0, 50)) +
  scale_color_manual(values = colors$Color) +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1, vjust = 0.25),
        axis.text.y = element_text(size = 6.5),
        legend.position = "none", 
        legend.title = element_blank(),
        strip.text.y = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold"))

# Legend 
df <- x_Eco
colors <- set_color(df) %>% mutate(Color = as.character(Color))
leg <- ggplot(df, aes(Type, Quant, color = Software, group = Software)) +
  geom_point() +
  geom_line(size = 0.5) +
  scale_color_manual(values = colors$Color) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "bottom", 
        legend.text = element_text(size = 8))
leg <- get_legend(leg)  
leg <- as_ggplot(leg)


pdf("Output_2_DataTreatment/Figures/1_Quant_Proteins.pdf", height = 5, width = 8)
grid.arrange(grobs = list(p1, p2, p3, p4, leg),
             ncol = 2, nrow = 3, 
             widths = c(2, 5),  
             heights = c(3, 3, 1),
             layout_matrix = rbind(c(1 ,2), 
                                   c(3, 4),
                                   c(5, 5)), 
             top = "PROTEINS\n")
dev.off()


## 2 - REPRODUCIBILITY - CV - PEPTIDES -----------------------
## Peptides Eco ----
CV_Eco_pep <- read.csv("./Output_2_DataTreatment/SuppTable_CV_pep.txt",
                       header = T, sep = " ", stringsAsFactors = T) %>%
  filter(Specie == "Ecoli") %>%
  dplyr::select(-Specie)
Exp <- CV_Eco_pep$Experience
CV_Eco_pep <- CV_Eco_pep %>%  
  dplyr::select(CV) %>%
  aggregate(., by = list(Exp), FUN = mean, na.rm = T) %>%
  mutate(Experience = Group.1) %>%
  separate(Group.1, c("Software", "Acquisition", "Extraction"), sep = "_")  %>%
  mutate(Software = factor(Software), 
         Acquisition = factor(Acquisition, levels = Lev_Acq) )
write.table(CV_Eco_pep, "./Output_2_DataTreatment/SuppTable_CV_pep_Eco.txt")


## Peptides UPS1 ----
CV_UPS_pep <- read.csv("./Output_2_DataTreatment/SuppTable_CV_pep.txt",
                       header = T, sep = " ", stringsAsFactors = T) %>%
  filter(Specie == "HUMAN") %>%
  dplyr::select(-Specie) %>%
  mutate(Combined = paste(Experience, Concentration, sep = "_")) 
Comb <- CV_UPS_pep$Combined
CV_UPS_pep <- CV_UPS_pep %>%  
  dplyr::select(CV) %>%
  aggregate(., by = list(Comb), FUN = mean, na.rm = T) %>%
  separate(Group.1, c("Software", "Acquisition", "Extraction", "Concentration"), sep = "_")  %>%
  mutate(Software = factor(Software), 
         Acquisition = factor(Acquisition, levels = Lev_Acq), 
         Concentration = factor(Concentration, levels = Lev_Conc))

x <- CV_UPS_pep %>%
  mutate(Exp = paste(Software, Extraction, Acquisition, sep = "_")) %>%
  dplyr::select(-c(Software, Acquisition, Extraction)) %>%
  spread(Concentration, CV)
write.table(x, "./Output_2_DataTreatment/SuppTable_CV_pep_UPS.txt")

# Eco fasta
df <- CV_Eco_pep  %>% 
  filter(Extraction == "Fasta") %>%
  mutate(Type = "Peptide")
colors <- set_color(df) %>% mutate(Color = as.character(Color))
p1 <- ggplot(df, aes(Type, CV, color = Software, group = Software)) +
  geom_point(size = 0.8) +
  facet_grid(rows = vars(Type), cols = vars(Acquisition)) +
  theme_minimal() +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Fasta\nE.coli", limits = c(5, 30)) +
  scale_color_manual(values = colors$Color) +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1, vjust = 0.25, color = "white"),
        axis.text.y = element_text(size = 6.5), 
        legend.position = "none", 
        legend.title = element_blank(),
        strip.text.y = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold")) 

# UPS fasta
df <- CV_UPS_pep  %>% 
  filter(Extraction == "Fasta") 
colors <- set_color(df) %>% mutate(Color = as.character(Color))
p2 <- ggplot(df, aes(Concentration, CV, color = Software, group = Software))  +
  geom_point(size = 0.8) +
  geom_line(size = 0.5) +
  facet_grid( cols = vars(Acquisition)) +
  theme_minimal() +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Fasta\nUPS1", limits = c(0, 70)) +
  scale_color_manual(values = colors$Color) +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1, vjust = 0.25),
        axis.text.y = element_text(size = 6.5),
        legend.position = "none", 
        legend.title = element_blank(),
        strip.text.y = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold")) 

# Eco library
df <- CV_Eco_pep  %>% 
  filter(Extraction == "Library") %>%
  mutate(Type = "Peptide")
colors <- set_color(df) %>% mutate(Color = as.character(Color))
p3 <- ggplot(df, aes(Type, CV, color = Software, group = Software)) +
  geom_point(size = 0.8) +
  facet_grid(rows = vars(Type), cols = vars(Acquisition)) +
  theme_minimal() +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Library\nE.coli", limits = c(5, 30)) +
  scale_color_manual(values = colors$Color) +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1, vjust = 0.25, color = "white"),
        axis.text.y = element_text(size = 6.5), 
        legend.position = "none", 
        legend.title = element_blank(),
        strip.text.y = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold")) 

# UPS library
df <- CV_UPS_pep  %>% 
  filter(Extraction == "Library") 
colors <- set_color(df) %>% mutate(Color = as.character(Color))
p4 <- ggplot(df, aes(Concentration, CV, color = Software, group = Software))  +
  geom_point(size = 0.8) +
  geom_line(size = 0.5) +
  facet_grid( cols = vars(Acquisition)) +
  theme_minimal() +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Library\nUPS1", limits = c(0, 70)) +
  scale_color_manual(values = colors$Color) +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1, vjust = 0.25),
        axis.text.y = element_text(size = 6.5),
        legend.position = "none", 
        legend.title = element_blank(),
        strip.text.y = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold")) 

# Legend 
df <- CV_Eco_pep %>% mutate(Type = "Peptides")
colors <- set_color(df) %>% mutate(Color = as.character(Color))
leg <- ggplot(df, aes(Type, CV, color = Software, group = Software)) +
  geom_point() +
  geom_line(size = 0.5) +
  scale_color_manual(values = colors$Color) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "bottom", 
        legend.text = element_text(size = 8))
leg <- get_legend(leg)  
leg <- as_ggplot(leg)


pdf("Output_2_DataTreatment/Figures/2_Reproducibility_CV_Peptides.pdf", height = 5, width = 8)
grid.arrange(grobs = list(p1, p2, p3, p4, leg),
             ncol = 2, nrow = 3, 
             widths = c(2, 5),  
             heights = c(3, 3, 1),
             layout_matrix = rbind(c(1 ,2), 
                                   c(3, 4),
                                   c(5, 5)), 
             top = "PEPTIDES\n")
dev.off()


## 2 - REPRODUCIBILITY - CV - PROTEINS -----------------------
## Proteins Eco ----
CV_Eco_prot <- read.csv("./Output_2_DataTreatment/SuppTable_CV_prot.txt",
                        header = T, sep = " ", stringsAsFactors = T) %>%
  filter(Specie == "Ecoli") %>%
  dplyr::select(-Specie)
Exp <- CV_Eco_prot$Experience
CV_Eco_prot <- CV_Eco_prot %>%  
  dplyr::select(CV) %>%
  aggregate(., by = list(Exp), FUN = mean, na.rm = T) %>%
  mutate(Experience = Group.1) %>%
  separate(Group.1, c("Software", "Acquisition", "Extraction"), sep = "_")  %>%
  mutate(Software = factor(Software), 
         Acquisition = factor(Acquisition, levels = Lev_Acq) )
write.table(CV_Eco_prot, "./Output_2_DataTreatment/SuppTable_CV_prot_Eco.txt")

## Proteins UPS1 ----
CV_UPS_prot <- read.csv("./Output_2_DataTreatment/SuppTable_CV_prot.txt",
                        header = T, sep = " ", stringsAsFactors = T) %>%
  filter(Specie == "HUMAN") %>%
  dplyr::select(-Specie) %>%
  mutate(Combined = paste(Experience, Concentration, sep = "_")) 
Comb <- CV_UPS_prot$Combined
CV_UPS_prot <- CV_UPS_prot %>%  
  dplyr::select(CV) %>%
  aggregate(., by = list(Comb), FUN = mean, na.rm = T) %>%
  separate(Group.1, c("Software", "Acquisition", "Extraction", "Concentration"), sep = "_")  %>%
  mutate(Software = factor(Software), 
         Acquisition = factor(Acquisition, levels = Lev_Acq), 
         Concentration = factor(Concentration, levels = Lev_Conc))

x <- CV_UPS_prot %>%
  mutate(Exp = paste(Software, Extraction, Acquisition, sep = "_")) %>%
  dplyr::select(-c(Software, Acquisition, Extraction)) %>%
  spread(Concentration, CV)
write.table(x, "./Output_2_DataTreatment/SuppTable_CV_prot_UPS.txt")

# Eco fasta
df <- CV_Eco_prot  %>% 
  filter(Extraction == "Fasta") %>%
  mutate(Type = "Protein")
colors <- set_color(df) %>% mutate(Color = as.character(Color))
p1 <- ggplot(df, aes(Type, CV, color = Software, group = Software)) +
  geom_point(size = 0.8) +
  facet_grid(rows = vars(Type), cols = vars(Acquisition)) +
  theme_minimal() +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Fasta\nE.coli", limits = c(0, 25)) +
  scale_color_manual(values = colors$Color) +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1, vjust = 0.25, color = "white"),
        axis.text.y = element_text(size = 6.5), 
        legend.position = "none", 
        legend.title = element_blank(),
        strip.text.y = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold")) 

# UPS fasta
df <- CV_UPS_prot  %>% 
  filter(Extraction == "Fasta") 
colors <- set_color(df) %>% mutate(Color = as.character(Color))
p2 <- ggplot(df, aes(Concentration, CV, color = Software, group = Software))  +
  geom_point(size = 0.8) +
  geom_line(size = 0.5) +
  facet_grid( cols = vars(Acquisition)) +
  theme_minimal() +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Fasta\nUPS1", limits = c(0, 65)) +
  scale_color_manual(values = colors$Color) +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1, vjust = 0.25),
        axis.text.y = element_text(size = 6.5),
        legend.position = "none", 
        legend.title = element_blank(),
        strip.text.y = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold")) 

# Eco library
df <- CV_Eco_prot  %>% 
  filter(Extraction == "Library") %>%
  mutate(Type = "Peptide")
colors <- set_color(df) %>% mutate(Color = as.character(Color))
p3 <- ggplot(df, aes(Type, CV, color = Software, group = Software)) +
  geom_point(size = 0.8) +
  facet_grid(rows = vars(Type), cols = vars(Acquisition)) +
  theme_minimal() +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Library\nE.coli", limits = c(0, 25)) +
  scale_color_manual(values = colors$Color) +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1, vjust = 0.25, color = "white"),
        axis.text.y = element_text(size = 6.5), 
        legend.position = "none", 
        legend.title = element_blank(),
        strip.text.y = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold")) 

# UPS library
df <- CV_UPS_prot  %>% 
  filter(Extraction == "Library") 
colors <- set_color(df) %>% mutate(Color = as.character(Color))
p4 <- ggplot(df, aes(Concentration, CV, color = Software, group = Software))  +
  geom_point(size = 0.8) +
  geom_line(size = 0.5) +
  facet_grid( cols = vars(Acquisition)) +
  theme_minimal() +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Library\nUPS1", limits = c(0, 65)) +
  scale_color_manual(values = colors$Color) +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1, vjust = 0.25),
        axis.text.y = element_text(size = 6.5),
        legend.position = "none", 
        legend.title = element_blank(),
        strip.text.y = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold")) 

# Legend 
df <- CV_Eco_prot %>% mutate(Type = "Proteins")
colors <- set_color(df) %>% mutate(Color = as.character(Color))
leg <- ggplot(df, aes(Type, CV, color = Software, group = Software)) +
  geom_point() +
  geom_line(size = 0.5) +
  scale_color_manual(values = colors$Color) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "bottom", 
        legend.text = element_text(size = 8))
leg <- get_legend(leg)  
leg <- as_ggplot(leg)

pdf("Output_2_DataTreatment/Figures/2_Reproducibility_CV_Proteins.pdf", height = 5, width = 8)
grid.arrange(grobs = list(p1, p2, p3, p4, leg),
             ncol = 2, nrow = 3, 
             widths = c(2, 5),  
             heights = c(3, 3, 1),
             layout_matrix = rbind(c(1 ,2), 
                                   c(3, 4),
                                   c(5, 5)), 
             top = "PROTEINSS\n")
dev.off()


## 3 - LINEARITY - R2 ----------------------------------------
Lin <- read.csv("./Output_2_DataTreatment/SuppTable_Linearities.txt", header = T, sep =" ", stringsAsFactors = T)

# Fasta
df <- Lin %>%
  separate(Experience, c("Software", "Acquisition", "Extraction"), sep = "_") %>%
  mutate(Software = factor(Software), 
         Acquisition = factor(Acquisition, levels = Lev_Acq)) %>%
  filter(Extraction == "Fasta") 

colors <- set_color(df) %>% mutate(Color = as.character(Color))

mu <- Lin %>%
  ddply(., "Experience", summarise, grp.mean = round(mean(R2), 3)) %>%
  separate(Experience, c("Software", "Acquisition", "Extraction"), sep = "_") %>%
  mutate(Software = factor(Software),
         Acquisition = factor(Acquisition)) %>%
  filter(Extraction == "Fasta") 

p1 <- ggplot(df, aes(R2,  fill = Software)) +
  geom_density() +
  geom_vline(data = mu, aes(xintercept = grp.mean),  linetype = "dashed", size = 0.7)  +
  geom_text(data = mu,  size = 3, inherit.aes = F, aes(x = 0.7, y = 15, label = grp.mean)) +
  facet_wrap( Acquisition  ~ Software, scales = "free", ncol = 4) +
  scale_x_continuous(limits = c(0.6, 1.01)) +
  labs(x = parse(text = "r^{2}"), y = "Density" ) +
  scale_fill_manual(values = colors$Color, label = colors$Software) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.5),
        panel.border = element_rect(color = 'black', fill = NA),
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Fasta") 


# Library
df <- Lin %>%
  separate(Experience, c("Software", "Acquisition", "Extraction"), sep = "_") %>%
  mutate(Software = factor(Software), 
         Acquisition = factor(Acquisition, levels = Lev_Acq)) %>%
  filter(Extraction == "Library") 

colors <- set_color(df) %>% mutate(Color = as.character(Color))

mu <- Lin %>%
  ddply(., "Experience", summarise, grp.mean = round(mean(R2), 3)) %>%
  separate(Experience, c("Software", "Acquisition", "Extraction"), sep = "_") %>%
  mutate(Software = factor(Software),
         Acquisition = factor(Acquisition, levels = Lev_Acq)) %>%
  filter(Extraction == "Library") 

p2 <- ggplot(df, aes(R2,  fill = Software)) +
  geom_density() +
  geom_vline(data = mu, aes(xintercept = grp.mean),  linetype = "dashed", size = 0.7)  +
  geom_text(data = mu,  size = 3, inherit.aes = F, aes(x = 0.7, y = 15, label = grp.mean)) +
  facet_wrap( Acquisition  ~ Software, scales = "free", ncol = 5) +
  scale_x_continuous(limits = c(0.6, 1.01)) +
  labs(x = parse(text = "r^{2}"), y = "Density" ) +
  scale_fill_manual(values = colors$Color, label = colors$Software) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.5),
        panel.border = element_rect(color = 'black', fill = NA),
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Library") 


pdf("Output_2_DataTreatment/Figures/3_Linearity-r2.pdf", height = 8.5, width = 12.5)
annotate_figure(
  ggarrange(p1, p2, 
            ncol = 2, nrow = 1, 
            widths = c(4, 5))
)
dev.off()

## 4 - ACCURACY - MAPE ----------------------------------------
tab <- read.table("Input_1_Other/Tab_ratio.csv", sep = ";", header = T, stringsAsFactors = F)

df <- read.csv("Output_2_DataTreatment/SuppTable_MAPE_UPS_all.txt", sep = " ") %>%
  dplyr::rename_all(function(.) {c("Experience", tab$Ratio)} ) %>%
  gather(Ratio, MAPE, -Experience) %>%
  separate(col = Experience, into = c("Software", "Acquisition", "Extraction"), sep = "_") %>%
  mutate(Ratio = factor(Ratio, levels = as.character(tab$Ratio)),
         Software = factor(Software),
         Acquisition = factor(Acquisition, levels = Lev_Acq),
         R1 = Ratio %>% as.character %>% gsub("fmol.*", "", .) %>% as.numeric, 
         R2 = Ratio %>% as.character %>% gsub("fmol", "", .) %>% gsub(".*\\/", "", .) %>% as.numeric, 
         ExpR = paste0(R1, "  /  ", R2,  "      (", R1 / R2, ")") %>% factor(., levels = unique(.)))

colors <- set_color(df) %>% mutate(Color = as.character(Color))

pdf("Output_2_DataTreatment/Figures/4_Accuracy_MAPE.pdf", height = 7, width = 9)
ggplot(df, aes(ExpR, MAPE, group = Software, color = Software, fill = Software)) +
  geom_point(size = 1) +
  geom_line() + 
  scale_color_manual(values = colors$Color) +
  scale_fill_manual(values = colors$Color) +  
  scale_y_continuous(name = "", limits = c(0, 126)) +
  scale_x_discrete(name = "") +
  theme_minimal() +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25, size = 8), 
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9),        
        legend.position = "bottom", 
        legend.title = element_blank()) +
  facet_grid(Acquisition ~ Extraction)
dev.off()

## 5 - SENSITIVITY - AUC ----------------------------------------
tab <- read.table("Input_1_Other/Tab_ratio.csv", sep = ";", header = T, stringsAsFactors = F)

df <- read.csv("Output_2_DataTreatment/SuppTable_AUC.txt", sep = " ") %>%
  gather(Experience, AUC, -Ratio) %>%
  separate(col = Experience, into = c("Software", "Acquisition", "Extraction"), sep = "_") %>%
  mutate(Software = gsub("\\.", "\\-", Software) %>% as.factor,
         Acquisition = factor(Acquisition, levels = Lev_Acq), 
         Ratio = factor(Ratio, levels = unique(Ratio)),
         R1 = Ratio %>% as.character %>% gsub("fmol.*", "", .) %>% as.numeric, 
         R2 = Ratio %>% as.character %>% gsub("fmol", "", .) %>% gsub(".*\\/", "", .) %>% as.numeric, 
         ExpR = paste0(R1, "  /  ", R2,  "      (", R1 / R2, ")") %>% factor(., levels = unique(.)))

colors <- set_color(df) %>% mutate(Color = as.character(Color))

pdf("Output_2_DataTreatment/Figures/5_AUC.pdf", height = 7, width = 9)
ggplot(df, aes(ExpR, AUC, group = Software, color = Software)) +
  geom_point(size = 1) +
  geom_line() + 
  scale_color_manual(values = colors$Color) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25, size = 8), 
        axis.text.y = element_text(size = 8),
        legend.position = "bottom", 
        legend.title = element_blank()) +
  facet_grid(Acquisition ~ Extraction)
dev.off()


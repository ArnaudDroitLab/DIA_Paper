# DIA_Paper

## 3 R scripts provided:

- 1_Cleanning allows filtering, perform (or not) intensity normalization, impute
  missing values and perform tables standardization across all software tools.
- 2_DataCompilation allows to compute the number of identifications and
  quantifications, coefficients of variation, calculate r2 for linearity
  curves, calculate MAPE for ratio accuracy and sensitivity and FDP for
  differential expression analysis.
- 3_Plots allows to draw all plots 

## General comments

- Please name all your software tool exports as followed:
  "SoftwareTool_AcquisitionMode_ExtractionMode" (example :
"Skyline_Narrow_Library.csv") and store them in a folder called
“Intput_1_Software/Results” (exemple : Input_1_Skyline/Results”)
- At the end of 1_Cleanning scripts, all the tables are in a standardized
  format and are containing the following columns: 
    - Experience : type "Software_ExtractionMode_AcquisitionMode"
    - Accession : "Accession_HUMAN" or "Accession_ECOLI"
    - ModifiedSequence :
	- xx[Carbamidomethyl]Cxx     
	- xx[Oxidation]Mxx      
	- xx[Deamidation]Nxx   
	- xx[Deamidation]Qxx     
    - Sequence : Stripped sequence
    - Precursor : Precursor charge
    - 24 columns called "UPS1_xfmol_injx", containing quantification values for
      each precursor
    - 8 columns called "RtoF_xfmol". If "TRUE", 3 quant values (over the 3
      replicates) were extracted by the tool, before missing value imputation.
    - 1 column called "RtoF_ID". If "TRUE", 1 quant value was extracted in at
      least on of the 24 samples, before missing value imputation. 
- Create a folder called “Input_1_Other” containing the files “transition.csv”,
  “UPS1_Accession_Name.txt” and “Tab_ratio.csv”.
- Create a folder called « Input_2_DataTreatment ». All the pre-processed
  tables from “1_Cleanning” steps will be stored in this folder. 
- Create a folder called “Output_2_DataTreatment”, containing 3 sub-folders
  called “Figures”, MissingValues” and “Normalisation”. All the tables
generated with the “2_DataCompilation” script will be generated in these
folders, and will be used to draw plots with the “3_Plots” script. 
>>>>>>> 4e80b3626f509433e9671d940118f417f86f2e03

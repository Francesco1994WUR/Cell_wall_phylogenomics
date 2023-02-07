#Detect copy number of cell wall genes across species

library(data.table)
library(readxl)
library(agricolae)
library(dplyr)
library(tidyr)


setwd("")

CW_Genes_all <- fread("", na.strings = ";")


#Load genome taxonomy and keep only angiosperm genes
Genome_taxonomy <- read_excel("")


#Annotate CW gene function
Function_files <- list.files(path = "", pattern = ".txt")

CW_Genes <- CW_Genes_all[0,]

for(f in 1:length(Function_files)){
  
  Current_file <- Function_files[f]
  
  #Load current function file
  setwd("")
  Curr_file <- fread(Current_file)
  
  #Subset CW genes belonging to current function file
  Genes_curr_file <- CW_Genes_all[CW_Genes_all$Gene %in% Curr_file$Gene,]
  
  #Annotate function name
  Function <- gsub(".txt", "", Current_file)
  Genes_curr_file$Broad_function <- Function
  
  #rbind to main dataframe
  CW_Genes <- rbind(CW_Genes, Genes_curr_file)
  
  #Remove unnecessary objects
  rm(Curr_file, Current_file, Genes_curr_file, PFAM_ann)
  
}

rm(CW_Genes_all, f, Function, Function_files)


#Set unique groups for analyses
Plant_groups <- unique(Genome_taxonomy$Group)


#Keep only CW genes whose species has 75% BUSCO
CW_Genes <- CW_Genes[CW_Genes$Species %in% Genome_taxonomy$Species_ID,]


#Count copy number of each CW gene function, per species
CN_functions <- as.data.frame(table(CW_Genes[,c("Broad_function", "Species")]))
CN_functions$Broad_function <- as.character(CN_functions$Broad_function)
CN_functions$Species <- as.character(CN_functions$Species)


#Save file
setwd("")
fwrite(CN_functions, "", sep = "\t")




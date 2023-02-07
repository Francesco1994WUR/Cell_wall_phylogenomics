#Identify clusters of CW genes along plant chromosomes

#Load packages
library(readxl)
library(data.table)
require("ff")
library(dplyr)


#Load CW genes

setwd("")

Gene_func_files <- list.files(path = ".", pattern = ".txt")


#Load BED files with progressive gene number on chromosome
BED_coord <- fread("")
colnames(BED_coord) <- c("Gene", "Chrom", "Start", "End", "Progressive_numb_on_chrom")


#Now we iterate on each CW gene file, we annotate coordinates, and we detect clustered genes

for(f in 1:length(Gene_func_files)){
  
  #Load CW gene file f
  setwd("")
  CW_file <- fread(Gene_func_files[f], na.strings = "NaNaN")
  
  #Annotate gene coordinates and progressive number on chrom
  CW_file <- left_join(CW_file, BED_coord, by = "Gene")
  
  #Order CW file by species, chrom, and gene number on chrom
  CW_file <- CW_file %>% arrange(Species, Chrom, Progressive_numb_on_chrom)
  
  #Calculate intergene dist
  For_diff_lat <- CW_file$Progressive_numb_on_chrom
  For_diff_lat <- For_diff_lat[2:length(For_diff_lat)]
  For_diff_lat <- c(For_diff_lat,NA)
  
  Dist_lat <- abs(CW_file$Progressive_numb_on_chrom-For_diff_lat)
  
  #Annotate intergene dist on CW file
  CW_file$Dist_lat <- Dist_lat
  
  #See at which row chromosomes change and put those cells as NAs
  CW_file[which(CW_file$Chrom != dplyr::lag(CW_file$Chrom))-1,13] <- NA
  
  #Set at NAs all the rows where diff between consecutive genes is > 3 (3 genes in between consecutive genes is threshold for clustering)
  CW_file[which(CW_file$Dist_lat > 3),13] <- NA
  
  #Fill distance of genes next to genes whose distance is still a number after filtering, with the previous number in distance column
  Row_numb_NAs <- which(!is.na(CW_file$Dist_lat))
  Diff <- c(diff(Row_numb_NAs, lag = 1),NA)
  
  Row_numb_NAs <- as.data.frame(cbind(Row_numb_NAs, Diff))
  
  Row_numb_NAs <- Row_numb_NAs[!Row_numb_NAs$Diff == 1,]
  
  Row_numb_NAs <- Row_numb_NAs$Row_numb_NAs
  
  Row_numb_NAs <- Row_numb_NAs+1
  
  CW_file[Row_numb_NAs,13] <- 1
  
  
  #Replace all numbers defining clusters with 1
  CW_file[which(!is.na(CW_file$Dist_lat)),13] <- 1
  CW_file[which(is.na(CW_file$Dist_lat)),13] <- 2
  
  
  CW_file$Clusters <- ifelse(x <- CW_file$Dist_lat == 1, cumsum(c(head(x, 1), tail(x, -1) - head(x, -1) == 1)), NA)
  
  
  Curr_function <- gsub(".txt", "", Gene_func_files[f])
  CW_file$Clusters <- paste(Curr_function, "Cluster", CW_file$Clusters, sep = "_")
  
  
  CW_file[which(grepl("_Cluster_NA", CW_file$Clusters)),14] <- NA
  
  
  #Keep relevant columns
  CW_file <- CW_file[,c(1,2,3,5,6,8,9,10,11,12,14)]
  
  
  #Save file
  setwd("")
  Curr_file <- Gene_func_files[f]
  fwrite(CW_file, file = Curr_file, sep = "\t")
  
}



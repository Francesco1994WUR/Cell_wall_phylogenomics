#Principal Component Analysis of gene copy number across angiosperms

library(factoextra)
library(readxl)
library(ggplot2)
library(gridExtra)
library(ggbiplot)
library(ggpubr)
library(dplyr)
library(openxlsx)
library(multicon)
library(tibble)



#Import gene copy number data per gene family
setwd("")

Data_all <- as.data.frame(read_excel(""))


#Import info on species taxonomy
Genome_taxonomy <- read_excel("")


#Compute PCAs using prcomp functions
PCA <- prcomp(Data_all[,1:151], center = TRUE, scale. = TRUE)


#Perform Horn method to know significant eigenvalues (=PCs to retain; sims is number of bootstraps used)
Eig_check <- get_eigenvalue(PCA)
fviz_eig(PCA, addlabels = TRUE, ylim = c(0, 50))


#PCA plot
PCA_plot <- ggbiplot(PCA,
         choices = c(1,2),
         ellipse=FALSE,
         groups=Data_all$Broad_clade,
         var.axes = FALSE) +
  theme(aspect.ratio = 1,
        legend.background = element_blank(),
        legend.direction = 'vertical',
        legend.position = 'right',
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "white"),
        panel.background = element_rect(fill = "grey 93"),
        panel.border = element_rect(colour = "black", fill = NA)) + 
  scale_color_manual(values=c("black", "black", "black", "black")) +
  scale_fill_manual(values = c("red", "dark grey", "yellow", "orange")) + # just offset by one to show
  geom_point(size = 3.5, shape = 21, aes(fill = groups, color = groups)) +
  scale_x_continuous(name="PC1 (39%)", limits=c(-2.2,3.5), breaks = c(-2,-1,0,1,2,3)) +
  scale_y_continuous(name="PC2 (7%)", limits=c(-3.5,1.3), breaks = c(-3,-2,-1,0,1))

PCA_plot


setwd("C:/Users/panca001/OneDrive - Wageningen University & Research/WUR/PhD/03_Data/12_Phylogenomics_cell_wall/07_PCA_Gene_PFAM_CN/")
ggsave(PCA_plot, filename="PCA_CN_plot_new.tiff", height=7, width=7, units="in", dpi=300, limitsize = FALSE)  


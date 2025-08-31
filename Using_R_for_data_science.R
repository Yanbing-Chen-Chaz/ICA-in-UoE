#put all files in a directory
#library packages I need, if the user didn't install them, install them
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
library(tidyverse)
library(pheatmap)

#set it as working directory and import data
setwd("C:/Users/75695/Desktop/RDS ICA") #change it to your directory
data_all <- read.csv("data_all.csv")
gene_annotation <- read.csv("gene_annotation.csv")
sample_annotation <- read.csv("sample_annotation.csv")
genelist <- read.table("genelist_3.txt")

#Use dplyr to convert the data to log scaled data
#Choose all columns having values, and check is there any NA values
cols_to_log <- c("A", "B","C", "D", "E", "F", "G", "H", "I", "J", "K", "L")
values <- data_all[, cols_to_log] 
if (any(is.na(values))) {
  stop("Error: NA values detected in `data_all`.")
}

# Check if any value <= 0, in that case should use log2(x+1) to avoid error
if (any(values <= 0)) {
  print("Some columns contain zero or negative values. Use log2(. + 1).")
  data_all <- data_all %>%
  mutate(across(all_of(cols_to_log), ~ log2(. + 1)))
 } else {
  print("No zero or negative values. Use log2(.).")
  data_all <- data_all %>%
  mutate(across(all_of(cols_to_log), ~ log2(.)))
 }

# Choose genes in genelist as my_data
my_data <- data_all %>%
  filter(X %in% genelist$V1)

#Choose annotations in genelist as listed_gene_annotation, use LongName to rename the rownames
listed_gene_annotation <- gene_annotation %>%
  filter(X %in% genelist$V1) 
rownames(my_data) <- listed_gene_annotation$LongName

#To plot, make a data frame without the column "X"
data_to_plot <- my_data[, !colnames(data_all) %in% "X"]

# Create annotation for rows (genes)
row_annotation <- data.frame(Type = listed_gene_annotation$Type)
rownames(row_annotation) <- listed_gene_annotation$LongName

# Create annotation for columns (samples)
col_annotation <- data.frame(Group = sample_annotation$TreatmentGroup)
rownames(col_annotation) <- sample_annotation$SampleName

#Make the first heatmap
pheatmap(
  data_to_plot,
  annotation_row = row_annotation,
  annotation_col = col_annotation,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Heatmap with Gene and Sample Clustering",
  cellwidth = 15,
  cellheight = 10,
  fontsize_row = 10,
  angle_col = 0
)

#Make the heatmap with only genes clustered 
pheatmap(
  data_to_plot,
  annotation_row = row_annotation,
  annotation_col = col_annotation,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "Heatmap with Gene Clustering",
  cellwidth = 15,
  cellheight = 10,
  fontsize_row = 10,
  angle_col = 0
)

# Create synthetic data
names <- c(rep("group1", 120) , rep("group2", 120) , rep("group3", 120), rep("group4", 120))
value <- unlist(data_to_plot, use.names = FALSE)
data <- data.frame(names,value)

#Draw a boxplot
boxplot(data$value ~ data$names, data = data,
        col = c("lightpink", "lightgreen", "lightblue", "lightcyan"),
        xlab ="Samples",ylab ="Log scaled gene expression level",
        main = " Gene expression differences between treatment groups")

# Represent points on the plot
stripchart(data$value ~ data$names,
           data = data,
           method = "jitter",
           pch = 19,
           col = 2:5,
           vertical = TRUE,
           add = TRUE)

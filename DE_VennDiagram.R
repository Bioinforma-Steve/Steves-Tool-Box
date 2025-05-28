#!/usr/bin/env Rscript

###### This script can be run from the command line
### Example execution: 
#### Rscript DE_VennDiagram.R -g genes_of_interest.txt -f group1.csv,group2.csv,group3.csv -e BD,BH,BL -o VennDiagram.pdf

library(optparse)

option_list = list(
  make_option(c("-g", "--genes"), type="character", default=NULL, 
              help="Text file containing genes of interest separated by \n", metavar="character"),
  make_option(c("-f", "--files"), type="character", default=NULL, 
              help="Comma-separated list filenames for Differential Expression result CSVs", metavar="character"),
  make_option(c("-e", "--experiments"), type="character", default=NULL, 
              help="Comma-separated list of names for comparison", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="VennDiagram.pdf", 
              help="output file name for pdf [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#read in genes of interest
goi <- read.csv(opt$genes,
                   header = F)

CSVs <- as.list(strsplit(opt$files, ",")[[1]])
ex_lst <- as.list(strsplit(opt$experiments, ",")[[1]])

#initialize a lst for dataframes
df_lst <- list()

#read in list of DGE results to dataframes and store as a list
for (csv in CSVs) {
  df <- read.csv(csv)
  #Subset data frame using the genes of interest
  df <- subset(df, Geneid %in% goi$V1)
  df <- rownames(df[df$pvalue < 0.05 & abs(df$log2FoldChange) > 1, ])
  df_lst[[csv]] <- df
}

#Venn Diagram
library(VennDiagram)


#build the venn diagram
vd <- venn.diagram(df_lst, 
                   filename = NULL,
                   height = 6, width = 6,
                   scaled = TRUE,
                   category.names = ex_lst)

#Write the plot to a pdf
library(grDevices)

pdf(file=opt$out)
    grid.draw(vd)
dev.off()
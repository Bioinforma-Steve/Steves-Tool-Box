#!/usr/bin/env Rscript

###### This script can be run from the command line
### Example execution:
#### Rscript KEGG_Enrichment_ORA.R -O fgr -g genes_of_interest.txt -f DE_results.csv -o Output_Directory

#### Example execution with filtering DE results enabled 
#### Rscript KEGG_Enrichment_ORA.R -O fgr -g genes_of_interest.txt -f DE_results.csv -F TRUE -o Output_Directory


rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Load relevant libraries
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(clusterProfiler) # for PEA analysis
library(grDevices) #Write the plot to a pdf
library(optparse)
library(grid)
library(ggplot2)
library(enrichplot)

print("packages loaded")

option_list = list(
  make_option(c("-O", "--organism"), type="character", default=NULL, 
              help="KEGG organism code (ex: fgr)", metavar="character"),
  make_option(c("-g", "--genes"), type="character", default=NULL, 
              help="Text file containing genes of interest separated by \n", metavar="character"),
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Filename of Differential Expression results (.csv)", metavar="character"),
  make_option(c("-F", "--Filter"), type="logical", default=FALSE, 
              help="Filter out non-significantly expressed genes (threshold p < 0.05)", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="KEGG_Enrichment_Out", 
              help="Name or Path to output directory [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#check for out option, if not provided, make the default
if (opt$out == "KEGG_Enrichment_Out") {
  dir.create("KEGG_Enrichment_Out", mode = "0777")
  print("No directory provided for output, created default")
  out <- "KEGG_Enrichment_Out"
} else { 
  if (!file.exists(opt$out)) { #If out is provided, check to see if it's created
    print("Provided directory does not exist, created it!")
    dir.create(opt$out) #if not, create it
    out <- opt$out
  }
  else {
    out <- opt$out
  }
}

#Load DE results
df <- read.csv(opt$file)

###### Over-representation Analysis########

#Subset DE results using genes of interest
#read in genes of interest
goi <- read.csv(opt$genes,
                header = F)

df.sub <- subset(df, Geneid %in% goi$V1)

# Annotate according to differential expression
df.sub <- df.sub %>% mutate(diffexpressed = case_when(
  log2FoldChange > 0 & padj < 0.05 ~ 'UP',
  log2FoldChange < 0 & padj < 0.05 ~ 'DOWN',
  padj > 0.05 ~ 'NO'
))

# Remove non-significant genes
if (opt$Filter) {
  print("Filtering enabled, non-significantly expressed genes filtered out")
  df.sub <- subset(df.sub, diffexpressed != 'NO')
}

enrich <- enrichKEGG(gene = df.sub$Geneid, organism = opt$organism)

#Make a dotplot
p1 <- dotplot(enrich, showCategory = 15)

pdf(file=paste(out, sep="/", "dotplot.pdf"))
grid.draw(p1)
dev.off()

# Barplot
bp <- barplot(enrich, showCategory = 20) 
bp2 <- mutate(enrich, qscore = -log(p.adjust, base = 10)) %>% barplot(x = "qscore")

pdf(file=paste(out, sep="/",  "Barplot.pdf"))
grid.draw(bp)
dev.off()

pdf(file=paste(out, sep="/",  "Barplot2.pdf"))
grid.draw(bp2)
dev.off()

#Cnetplot
d <- setNames(df.sub$log2FoldChange, df.sub$Geneid)
cnet <- cnetplot(enrich, categorySize = "pvalue", foldChange=d)

pdf(file=paste(out, sep="/", "CNetPlot.pdf"))
grid.draw(cnet)
dev.off()

# calculate pairwise similarities of the enriched terms using Jaccard’s similarity index
enrichres2 <- pairwise_termsim(enrich)

# Enrichment map 
em <- emapplot(enrichres2)

pdf(file=paste(out, sep="/", "Enrichment_Map.pdf"))
grid.draw(em)
dev.off()

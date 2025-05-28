#!/usr/bin/env Rscript

###### This script can be run from the command line
### Example execution:
#### Rscript GO_Enrichment_Analysis.R -g genes_of_interest.txt -f DE_results.csv -d BP -t2n TERM2NAME.tsv -o Output_Directory

#### Example execution with filtering DE results enabled 
#### Rscript GO_Enrichment_Analysis.R -g genes_of_interest.txt -f DE_results.csv -d BP -t2n TERM2NAME.tsv -F TRUE -o Output_Directory

#### This script requires three look up tables containing the GO IDs and their associated terms to be in the same directory where the
#### the script is executed. I built them using the GO.db library. I've separated them by database type.
#### They are labeled: godb_BP.txt, godb_MF.txt, and godb_CC.txt

rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Load relevant libraries
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(pheatmap)
library(clusterProfiler) # for PEA analysis
library(enrichplot) # for visualisations
library(ggupset) # for visualisations
library(grDevices) #Write the plot to a pdf
library(optparse)
library(grid)

print("packages loaded")

option_list = list(
  make_option(c("-g", "--genes"), type="character", default=NULL, 
              help="Text file containing genes of interest separated by \n", metavar="character"),
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Filename of Differential Expression results (.csv)", metavar="character"),
  make_option(c("-d", "--GO_Database"), type="character", default=NULL, 
              help="Select GO Databse: Biological Process = BP; Molecular Function = MF; Cellular Component = CF", metavar="character"),
  make_option(c("-t", "--TERM2GENE"), type="character", default=NULL, 
              help="TSV containing GO IDs in first column and corresponding Gene IDs in the second", metavar="character"),
  make_option(c("-F", "--Filter"), type="logical", default=FALSE, 
              help="Filter out non-significantly expressed genes (threshold p < 0.05)", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="GO_Enrichment_Out", 
              help="Name or Path to output directory [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#check for out option, if not provided, make the default
if (opt$out == "GO_Enrichment_Out") {
  dir.create("GO_Enrichment_Out", mode = "0777")
  print("No directory provided for output, created default")
  out <- "GO_Enrichment_Out"
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


#read in genes of interest
goi <- read.csv(opt$genes,
                header = F)

#Load DE results
df <- read.csv(opt$file)

#Subset DE results using genes of interest
df <- subset(df, Geneid %in% goi$V1)

# Annotate according to differential expression
df <- df %>% mutate(diffexpressed = case_when(
  log2FoldChange > 0 & padj < 0.05 ~ 'UP',
  log2FoldChange < 0 & padj < 0.05 ~ 'DOWN',
  padj > 0.05 ~ 'NO'
))

# Remove non-significant genes
if (opt$Filter) {
  print("Filtering enabled, non-significantly expressed genes filtered out")
  df <- subset(df, diffexpressed != 'NO')
}

#make GO ID-Gene ID mapping table
#read in the gene annotation file
goid2gene <- read.table(opt$TERM2GENE, header = T,
                        sep = "\t", quote = "", stringsAsFactors = F)

#Check to see if the first col of the TERM2GENE file is "go_id", if not change it
if (colnames(goid2gene)[1] != "go_id") {
  colnames(goid2gene)[1] <- "go_id"
}

#read in the GO annotation files by type
if (opt$GO_Database == "BP") {
  godb <- read.table("godb_BP.txt", header = T,
                     sep = "\t", quote = "", stringsAsFactors = F)
}
if (opt$GO_Database == "MF") {
  godb <- read.table("godb_MF.txt", header = T,
                     sep = "\t", quote = "", stringsAsFactors = F)
}
if (opt$GO_Database == "CC") {
  godb <- read.table("godb_CC.txt", header = T,
                     sep = "\t", quote = "", stringsAsFactors = F)
}

#separate genes based on their GO type
goid2gene <- goid2gene %>% filter(go_id %in% godb$go_id)

print("Sanity Checkpoint 1")

#make enrich objects for our genes of interest
enrich <- enricher(gene = df$Geneid,
                   TERM2GENE = goid2gene,
                   TERM2NAME = godb)
data.frame(enrich)
#Make a dotplot
p1 <- dotplot(enrich, showCategory = 10)

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
d <- setNames(df$log2FoldChange, df$Geneid)
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

# Treeplot
tp <- treeplot(enrichres2)

pdf(file=paste(out, sep="/", "Treeplot.pdf"))
grid.draw(tp)
dev.off()

print("Enrichment analysis completed!")

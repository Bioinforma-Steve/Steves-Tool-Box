import argparse
import os
import pickle as pkl
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

#Initialize the 
parser = argparse.ArgumentParser(description="A program to construct heatmaps using pyDeSeq2 derived differential expression results for genes of interest")
parser.add_argument("--DDS", help="A dds pickle object from a DeSeq differential expression analysis")
parser.add_argument("--goi", help="Name of the text file containing the genes of interest")
parser.add_argument("--name", help="Name of the heatmap file ex: heatmap.pdf")
parser.add_argument("--out", help="Path of directory where to send the output")

args = parser.parse_args()

#Set out path
OUTPUT_PATH = args.out

#load the pickle file containing the deseq results
file = open(args.DDS,'rb')
dds = pkl.load(file)
print("-------------DeSeqDataSet loaded-------------")

#Extract the normalized counts and add them to the dds object as a log transformed layer
dds.layers['log1p'] = np.log1p(dds.layers['normed_counts'])

#Transpose transformed counts and create a dataframe
grapher = pd.DataFrame(dds.layers['log1p'].T,
                       index=dds.var_names, columns=dds.obs_names)
print("-------------Graph dataframe created-------------")

#Load text file containing the genes of interest
txt = open(args.goi, "r")
goi = txt.read()
print("-------------Genes of Interest loaded-------------")

#Convert the genes of interest into a list
goi = goi.split('\n')
print("-------------Gene list created-------------")

#Close the genes of interest text file
txt.close()

#subset the differential expression results using the genes of interest
grapher_sub = grapher.loc[grapher.index.isin(goi)]
print("-------------Graph dataframe subset-------------")

#Build the plot
hm = sns.clustermap(grapher_sub, z_score=0, cmap = 'RdYlBu_r')
hm.ax_heatmap.set_xlabel('Treatment')
hm.ax_heatmap.set_ylabel('Gene')

#Write the plot to file
hm.fig.savefig(os.path.join(OUTPUT_PATH, args.name), format='pdf')
print("Heatmap created, program finished!")
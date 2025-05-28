"""Convert a GFF to a BED file.
Usage:
    gff_2_BED.py --GFF Input.gff --BED Out.bed --out /PATH/TO/OUT
    
GFF is the input file to convert

BED is the output

The script will output a BED file in the following format: Chromosome/Scaffold Start End Gene ID
 
"""

import argparse
import os
from BCBio import GFF
import gffpandas.gffpandas as gffpd
import pandas as pd

#Initialize the arguments
parser = argparse.ArgumentParser(description="A program to update a Funannotate produced GFF with a GFF containing BLAST2GO annotations from OmicsBox")
parser.add_argument("--GFF", help="A GFF file to be updated")
parser.add_argument("--BED", help="A GFF file used to update the other")
parser.add_argument("--out", help="Path for directory where output is deposited")

args = parser.parse_args()

out = args.out

#load in GFF to be converted
annotation = gffpd.read_gff3(args.GFF)
print("GFF loaded")

#Get only the genes
filtered_df = annotation.filter_feature_of_type(['gene'])

#convert the attributes to columns
attr_to_columns = filtered_df.attributes_to_columns()

#Get the needed columns
bed = attr_to_columns[['seq_id', 'start', 'end', 'ID']]

print("writing {} to {}".format(args.BED, out))

#write the BED to out
bed.to_csv(os.path.join(out, args.BED), sep='\t', header=False, index=False)
print("program completed")

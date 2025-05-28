"""Update a Funannotate produced GFF with a GFF containing BLAST2GO annotations from OmicsBox.
Usage:
    gff_2_gff_CDS_attributes_updater.py --GFF1 file.gff --GFF2 file.gff --name New_GFF.gff3 --out /PATH/TO/OUT
    
GFF1 that contains CDS features with incomplete annotations

GFF2 contains updated attributes for CDS features in GFF1

The script will output a new GFF file combining the structure of GFF1 with more comprehsive CDS attributes from GFF2
 
"""


import argparse
import pandas as pd
import os
from Bio import SeqIO
from BCBio import GFF
import gffpandas.gffpandas as gffpd

#Initialize the arguments
parser = argparse.ArgumentParser(description="A program to update a Funannotate produced GFF with a GFF containing BLAST2GO annotations from OmicsBox")
parser.add_argument("--GFF1", help="A GFF file to be updated")
parser.add_argument("--GFF2", help="A GFF file used to update the other")
parser.add_argument("--name", help="Name of the updated GFF file ex: File_updated.gff3")
parser.add_argument("--out", help="Path for directory where output is deposited")

args = parser.parse_args()

out = args.out

#load in GFF to be updated
annotation = gffpd.read_gff3(args.GFF1)
print("First GFF loaded")

#load in GFF with updated data
annotation2=gffpd.read_gff3(args.GFF2)
print("Second GFF loaded")

#Split the attributes column by the ";" up to the Description
annotation2_att_listed = annotation2.df['attributes'].str.split(";", n=2)

#Change "Description" to "product" for each attribute
for index, row in annotation2_att_listed.items():
    for item in range(len(row)):
        if 'Description=' in row[item]:
            row[item] = row[item].replace('Description=', 'product=')

#create a dictionary from gff annotations where keys are IDs and values are attributes starting with "Description"
CDS_dict = {i[1][0].replace('ID=', '') : ";".join(i[1][1:]) for i in annotation2_att_listed.items()}

#create a different dictionary to update the mRNA product attribute
mRNA_dict = {i[1][0].replace('ID=', '') : i[1][1].replace('product=', '') for i in annotation2_att_listed.items()}

# Iterate through the DataFrame rows and update CDS attributes by 
# appending the Blast2Go attributes to funannotate attributes
for index, row in annotation.df.iterrows():
    if row['type'] == "CDS":
        lst = row['attributes'].split(';')
        #print(lst)
        for item in lst:
            if 'Parent=' in item:
                id_value = item[7:]
                if lst[-1] == '':
                    lst[-1] = CDS_dict[id_value]
                else:
                    lst.append(CDS_dict[id_value])  
        att_str = ";".join(lst)
        annotation.df.loc[index, 'attributes'] = att_str + ";"
print("GFF attributes updated")

# Iterate through the DataFrame rows update the mRNA product description
for index, row in annotation.df.iterrows():
    if row['type'] == "mRNA":
        lst = row['attributes'].split(';')
        for item in range(len(lst)):
            if 'ID=' in lst[item]:
                id_value = lst[item][3:]
                #print(id_value)
            if 'product=' in lst[item] and id_value in mRNA_dict:
                lst[item] = lst[item].replace(lst[item][8:], mRNA_dict[id_value])
                #print(item)
        att_str = ";".join(lst)
        annotation.df.loc[index, 'attributes'] = att_str
print("GFF mRNA attributes updated")

#join the directory and output name
new_gff = os.path.join(out, args.name)
print("writing {} to {}".format(args.name, out))

#Write the new GFF to file
annotation.to_gff3(new_gff)
print("---Program Finished---")





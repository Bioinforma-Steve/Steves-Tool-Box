#!/usr/bin/env python
# coding: utf-8

"""
Construct a table from antiSMASH output that summarizes the detected clusters. 
For each cluster, the table contains genomic coordinates, BGC type, cluster type, 
top secondary metabolite (if applicable), the number of loci in cluster, 
and the gene models in the cluster.

Script successfully tested with versions 6 and 7. Not yet tested with antiSMASH 8.

Usage:

python antiSMASH_JSON-GBK_2_ClusterTable.py --input <antiSMASH Output Directory> --JSON <antiSMASH.json> --out <output.csv>
    
"""

import argparse
import re
import json
import pandas as pd
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

#Initialize the 
parser = argparse.ArgumentParser(description="A program to generate an antiSMASH cluster table in CSV format")
parser.add_argument("--input", type=str, help="path of the folder containing antiSMASH output, omit last '/'")
parser.add_argument("--JSON", type=str, help="JSON file output by antiSMASH")
parser.add_argument("--out", type=str, default="antiSMASH_Cluster_Table.csv", help="Name for the output cluster table")

args = parser.parse_args()

def getClusterAreas(JSON):
    lst = []
    for record in JSON['records']:
        rec_id = record['id']
        for area in record['areas']:
            lst.append([rec_id,area['start'],area['end'],area['products']])

    df = pd.DataFrame(lst, columns=['record','start','end','BGC_Type'])
    return(df)

def getTopSM(JSON):
    top_SM_lst = []

    for record in JSON['records']:
        rec_id = record['id']
        if 'antismash.modules.clusterblast' not in record['modules']:
            continue
        for result in record['modules']['antismash.modules.clusterblast']['knowncluster']['results']:
            if len(result['ranking']) < 1:    
                region = result['region_number']
                top_SM = ''
                top_SM_lst.append([rec_id,region,top_SM])
            else:
                region = result['region_number']
                top_SM = result['ranking'][0][0]['description']
                Cluster_Type = result['ranking'][0][0]['cluster_type']
                top_SM_lst.append([rec_id, region,Cluster_Type,top_SM])

    df2 = pd.DataFrame(top_SM_lst, columns=['record','region','Cluster_Type','Top_SM'])
    return(df2)

def mergeDFs(df1,df2):
    merged_df = pd.merge(df1, df2, how='inner', left_index=True, right_index=True, suffixes=('', '_drop'))
    merged_df.drop([col for col in merged_df.columns if 'drop' in col], axis=1, inplace=True)
    return(merged_df)

def getBGCLoci(gbk):
    CDS = []
    loc_list = []
    num_loc = 0

    for feature in gbk.features:
        if feature.type == "CDS":    
            if 'locus_tag' in feature.qualifiers:
                locus_tag = feature.qualifiers['locus_tag'][0]
                locus_tag = [locus_tag]
            #write to a new line using f-strings
            line = gbk.id,locus_tag
            CDS.append(line)

    CDS_df = pd.DataFrame(CDS, columns = ['Record','loci'])
    f = lambda x: list(sorted(set(z for y in x for z in y)))
    CDS_df = CDS_df.groupby(['Record'])['loci'].agg(f).reset_index()
    CDS_df['num_loc'] = [len(locus) for locus in CDS_df['loci']]
    #print(CDS_df)
    
    return(CDS_df)

def get_contig(STR):
    if re.search(r"contig|c0+", STR):
        STR = STR[0:-14]
        pattern = r"^(.*?)([_Ccontig]+)0*(\d+)$"
        match = re.match(pattern, STR)
        num = int(match.group(3))
        return num
    else:
        return STR[0:-14]

input_dir = args.input + '/'
regionGBK_files = [f for f in os.listdir(input_dir) if 'region' in f and f.endswith('.gbk')]
print("-------------GBKs Found-------------")

dat = []

#search for and loop through antiSMASH GBKs
for f in regionGBK_files:
    gbk = SeqIO.read(open(input_dir + f, 'r'),"genbank")
    gbk_loci = getBGCLoci(gbk)
    gbk_loci.insert(0,'Region', f'{f[-5:-4]}')
    gbk_loci.insert(1, 'contig', get_contig(f))
    dat.append(gbk_loci)
    df_GBKs = pd.concat(dat)

df_GBKs = df_GBKs.sort_values(by=['contig'], ascending = True)
df_GBKs = df_GBKs.reset_index(drop=True)

print(df_GBKs[:50])
print('-------------End of GBK Sanity Check-------------')

#bring in JSON
with open(os.path.join(input_dir, args.JSON)) as my_jf:
    json_text = json.load(my_jf)
    df1 = getClusterAreas(json_text)
    df2 = getTopSM(json_text)
    df_JSON = mergeDFs(df1,df2)
    df_JSON = df_JSON.reset_index(drop=True)
print('-------------JSON Loaded-------------')

print(df_JSON[:50])
print('-------------End of JSON Sanity Check-------------')

BGC_df = pd.merge(df_JSON,df_GBKs, how='inner', left_index=True, right_index=True)
BGC_df = BGC_df.drop('Record', axis=1)
BGC_df = BGC_df[['record','start','end','BGC_Type','Cluster_Type','region','Top_SM','num_loc','loci']]
print(BGC_df[:20])
print('-------------Cluster Table Constructed-------------')

BGC_df.to_csv(args.out, index = False)

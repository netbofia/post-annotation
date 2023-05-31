#! /usr/bin/env python3


######################################



import argparse
import glob
from pathlib import Path
import json
import pandas as pd
from tabulate import tabulate

#parser = argparse.ArgumentParser(prog="featureMapper", description="Generates a feature map from the miRCat hairpin files")
#parser.add_argument("-m", "--mircat", help="path to miRCat directory")
#args = parser.parse_args()

#mircat_folder = args.mircat
mircat_folder = "/home/brunocosta/Downloads/TEST_miRNAPlantPortal/mircat"
mircat_folder_path = Path(mircat_folder)


def tab(dataframe):
    # https://www.geeksforgeeks.org/display-the-pandas-dataframe-in-table-style/
    print(tabulate(dataframe, headers='keys', tablefmt='psql'))


hairpinFiles = None
if mircat_folder_path.is_dir():
    #Remove fixed in production
    mircat_files = glob.glob(mircat_folder+'/*_mirbase_noncons_output_filteredfixed.csv')

header="Chromosome,Start,End,Orientation,Abundance,Sequence,sRNA length,# Genomic Hits,Hairpin Length,Hairpin % G/C content,Minimum Free Energy,Adjusted MFE,miRNAstar,miRBaseID,pVal".split(",")
df_list = [pd.read_csv(file, skiprows=1, header=None, sep=",", on_bad_lines="warn") for file in mircat_files]
df = pd.concat(df_list)
df.columns = header
print(df.shape)
#Should be 75246-9 headers =75237 checks out fine
df = df.drop_duplicates()
#75237 --> 72506
print(df.shape)

#TODO unpack stars?
# recursively join miRNAs on same precursor

tab(df[df['miRNAstar'] != "NO"])
starDf = pd.DataFrame()
for starRow in df[df['miRNAstar'] != "NO"].itertuples():
    for miRNAstar in starRow.miRNAstar.strip().split(" "):
        starDf.loc[len(starDf), ["mature", "star","starAbundance" , "chromosome","start","end"]] = [starRow.Sequence, miRNAstar.split("(")[0], miRNAstar.split(")")[0].split("(")[1], starRow.Chromosome, starRow.Start, starRow.End]

tab(starDf)
# 3336 rows
# 1169 unique
print(starDf.duplicated().value_counts())
#Unique mature : star = 1988


## fasta of stars not in all_seq.tsv

df2 = pd.read_table("~/Downloads/TEST_miRNAPlantPortal/all_seq-sRNA-SuLop.tsv")
fw = open("absent_star_seq.tsv", "w")
for star in starDf.loc[~starDf.star.isin(df2.sequence), "star"].unique():
    mockExpression="\t".join(["0"] * 9)
    fw.write(f'{star}\tnovel-star\t{mockExpression}\n')


fw.flush()
fw.close()

for row in starDf.itertuples():
    print(star)
    print(len(df.loc[(df.Sequence == row.mature) & (df.Chromosome == row.chromosome ) & (df.Start == row.start) & (df.End == row.end) & (df.miRNAstar.str.contains(row.star)), :]))



#merge mature:star
#merge mature on same miRNA
#add star seq
#rename


features = None
scaffoldName = "NW_019827908.1"

file = open(file, "r")
header = file.readline().strip().split("\t")

if filePath.is_file() :
    features = [dict(zip(header,line.strip().split("\t"))) for line in file.readlines() if line.startswith(scaffoldName)]

i=0
export=[]
for i,feature in enumerate(features):
    export.append( {"name":f'miRNA{i}',"data":[{"chr":features[i]['chromossome'].split(" ",1)[0], "start": features[i]["Start"], "end": features[i]["End"],"type":"miRNA","sequence":features[i]["Sequence"],"MFE":features[i]["Minimum Free Energy"], "description": features[i]["chromossome"],"color":"red"  }]} )





fw=open("/home/brunocosta/git/sRNA-Portal-workflow/public/javascripts/features/qsuber-miRNA.json","w")
fw.write(json.dumps(export,indent=4))
fw.flush()
fw.close()



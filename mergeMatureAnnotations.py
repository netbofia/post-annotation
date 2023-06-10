#! /usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
pd.set_option('display.max_columns', 6)

mismatch = 2
df = pd.read_table(f'Result-m{mismatch}.tsv')

df.columns = ["target", "query", "start", "end", "orientation", "edits"]

df = df.astype({'edits': 'int'})
df = df.astype({'start': 'int'})
df = df.astype({'end': 'int'})

df['queryName'] = df['query'].str.split('_').str[1]
df['querySeq'] = df['query'].str.split('_').str[0]
df['targetName'] = df['target'].str.split('_').str[1]
df['targetSeq'] = df['target'].str.split('_').str[0]
df['family'] = ""


def add_name_to_query(dataf, query_sequence, counter):
    query_slice = dataf.loc[dataf["querySeq"] == query_sequence]
    for query_row in query_slice.itertuples():
        if query_row.family == "":
            dataf.loc[query_row.Index, "family"] = counter
            target_sequence = query_row.targetSeq
            dataf = add_name_to_query(dataf, target_sequence, counter)
    return dataf


novel = len(df[df["queryName"] == "novel"]["querySeq"].unique())

family = 0

# Get all sequences of novel miRNAs
for querySequence in df[df["queryName"] == "novel"]["querySeq"].unique():
    for mapping in df[df["querySeq"] == querySequence].itertuples():
        targetSequence = mapping.targetSeq
        if mapping.family == "":
            df.loc[df["querySeq"] == querySequence, "family"] = family
            df = add_name_to_query(df, targetSequence, family)
            print(family)
    ## Null + current count(0 indexed) + 1 -->family=0; unique={0,1}; len=2 means that a family was added in run
    if len(df["family"].unique()) > family+1:
        family += 1
print(df['accession'])



### Possible new stuff

piv = pd.pivot_table(df, "query", "accession", "queryName", aggfunc="count")['novel']
piv = pd.pivot_table(df, "queryName", "family", aggfunc="count")

## To make fig
dfNoDup = df[["query", "queryName", "querySeq", "family"]].drop_duplicates()
piv = pd.pivot_table(dfNoDup, "queryName", "family", aggfunc="count")

y = list(piv.values)
y.pop()
x = list(piv.index)
x.pop()
plt.title(f'Novel miRNAs that have been grouped\n into the same accession with mismatch=2')
plt.xlabel("accessions")
plt.ylabel("# mature miRNA")
plt.plot(x, y)
plt.savefig("Novel count same accession.png")

# Multiple conditionals
# df.loc[(df["accession"] != "") & (df["queryName"]=="novel"),:]


df2 = pd.read_table("~/Downloads/TEST_miRNAPlantPortal/all_seq-sRNA-SuLop.tsv")
df2.index = df2['sequence']

# 31372
dfNovelAll = df2[df2['name'] == "novel"]

# Import Excel
# pip install openpyxl
xlsx = pd.read_excel("~/Documents/Lab-People/Susana Lopes/Artigo sRNA/Supplementary Tables.xlsx",sheet_name=2,header=3)


# Exclusion attempt
expCols = [lib for lib in df2.columns if lib.startswith("Lib")]

labels = [["meta"] * 2 + ["abs"] * len(expCols) + ["bool"] * len(expCols), ["sequence", "name"] + expCols * 2]
tuples = list(zip(*labels))
cols = pd.MultiIndex.from_tuples(tuples, names=["Type", "Name"])

for lib in expCols:
    dfNovelAll.loc[lib+"_bool"] = 0
    dfNovelAll.loc[dfNovelAll[lib] > 0, lib+"_bool"] = 1

dfNovelAll.columns = cols


dfNovelAll.loc[:, 'mapped'] = False
dfNovelAll.loc[dfNovelAll['sequence'].isin(df.querySeq.unique()), 'mapped'] = True
dfNovelAll.loc[:, 'susana'] = False
susanaSeq = xlsx[xlsx.names == "novel"].Sequence.unique()
dfNovelAll.loc[dfNovelAll.sequence.isin(susanaSeq), "susana"] = True

# pip install matplotlib-venn
from matplotlib_venn import venn2

venn2([set(dfNovelAll[dfNovelAll.mapped].index), set(dfNovelAll[dfNovelAll.susana].index)], ("Mapped", "Susana"))
plt.show()

from matplotlib_venn import venn3

venn3([set(dfNovelAll.index), set(dfNovelAll[dfNovelAll.mapped].index), set(dfNovelAll[dfNovelAll.susana].index)], ("All", "Mapped", "Susana"))
plt.title("Novel miRNAs Qsuber dataset")
plt.savefig("Venn diagram all novel Susana.png")

dfNovelAll.loc[:, "family"] = None

mapped = df.loc[df.queryName == "novel", ['querySeq', "family"]].drop_duplicates()


for row in mapped.itertuples():
    dfNovelAll.loc[row.querySeq, 'family'] = row.family

novelCounter = df[df.family != ""].family.astype("int").max()
for row in dfNovelAll.itertuples():
    if row.family is None:
        novelCounter += 1
        dfNovelAll.loc[row.sequence, 'family'] = novelCounter

#From 31372 to 25077
dfNovelAll.to_excel("Novel_with_family.xlsx")



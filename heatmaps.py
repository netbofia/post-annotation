import pandas as pd
# pip install seaborn
import seaborn as sns
from matplotlib import pyplot as plt

sts7 = pd.read_excel("~/Documents/Lab-People/Susana Lopes/Artigo sRNA/Supplementary Tables.xlsx", sheet_name=9, header=3)
cons = sts7.loc[sts7['names'] != "novel"]
cons.index = cons.names
data = cons.loc[:, [col for col in cons.columns if col.startswith("Lib")]]
data = data[data.index != "tasi"]




# Identify columns
expCols = [lib for lib in data.columns if lib.startswith("Lib")]


labels = [["TDP"] * 3 + ["DX"] * 3, expCols]
tuples = list(zip(*labels))
cols = pd.MultiIndex.from_tuples(tuples, names=["Stage", "Library"])
data.columns = cols



# Start a new DataFrame with the row sums
families = data.index.str.split("_").str[0]
firstFamily = families.unique()[0]
hmData = pd.DataFrame([data[data.index.str.split("_").str[0] == firstFamily].sum()], columns=cols)
hmData.index = [firstFamily]

for family in data.index.str.split("_").str[0].unique()[1::]:
    hmData.loc[len(hmData)] = data[data.index.str.split("_").str[0] == family].sum()

hmData.index = data.index.str.split("_").str[0].unique()

hmData.loc[(hmData.index == "mir165") | (hmData.index == "mir166"), :].sum()
# Add row merged of mir165/166
hmData.loc[len(hmData)] = hmData.loc[(hmData.index == "mir165") | (hmData.index == "mir166"), :].sum()
newIndex = list(hmData.index[0:-1])
newIndex.append("mir165/166")
hmData.index = newIndex
hmData = hmData.drop(["mir165", "mir166"])

col_colors = ["r", "r", "r", "g", "g", "g"]

import math
hmData = hmData.apply(lambda x: x.apply(lambda y: math.log(y+1)))
sns.clustermap(hmData, metric="canberra", method="average", dendrogram_ratio=0.20, cmap="inferno")
plt.savefig("Fig 4(a).tif", dpi=300, pil_kwargs={"compression": "tiff_lzw"})

hmData = hmData.apply(lambda x: x.apply(lambda y: math.log(y+1)))
hmData = hmData.drop(["mir-92", "mir-1260"])
sns.clustermap(hmData, metric="canberra", method="average", dendrogram_ratio=0.20, cmap="inferno")
plt.savefig("Fig 4(b).tif", dpi=300, pil_kwargs={"compression": "tiff_lzw"})


data = data.apply(lambda x: x.apply(lambda y: math.log(y+1)))
sns.clustermap(data, metric="canberra", method="average", dendrogram_ratio=0.20, cmap="inferno")
plt.savefig("Fig S2(a).tif", dpi=300, pil_kwargs={"compression": "tiff_lzw"})

data = data.apply(lambda x: x.apply(lambda y: math.log(y+1)))
data = data.drop(["mir-92_1", "mir-1260_1"])
sns.clustermap(data, metric="canberra", method="average", dendrogram_ratio=0.20, cmap="inferno")
plt.savefig("Fig S2(b).tif", dpi=300, pil_kwargs={"compression": "tiff_lzw"})

allSeq = sts7.loc[:, [col for col in sts7.columns if col.startswith("Lib")]]
allSeq.index = sts7.Sequence
allSeq.columns = cols
allSeq = allSeq.apply(lambda x: x.apply(lambda y: math.log(y+1)))
sns.clustermap(allSeq, metric="canberra", method="average", dendrogram_ratio=(0.2,0.2), cmap="inferno", figsize=(9, 70), cbar_kws={"shrink": 0.05, "pad": 0.5, "orientation": "vertical"})
#plt.show()
plt.savefig("Fig S3.tif", dpi=300,  pil_kwargs={"compression": "tiff_lzw"})




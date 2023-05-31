#! /usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
from tabulate import tabulate
import xlsxwriter
#for i in vs/2/*;do cat $i/result >> Result-m2.tsv;done

mismatch = 2

df = pd.read_table(f'Result-m{mismatch}.tsv')
print(df.head())
df.columns = ["target", "query", "start", "end", "orientation", "edits"]

df = df.astype({'edits': 'int'})
df = df.astype({'start': 'int'})
df = df.astype({'start': 'int'})

df['queryName'] = df['query'].str.split('_').str[1]
df['targetName'] = df['target'].str.split('_').str[1]

#Filter Cons

pivot = pd.pivot_table(df, "query", "targetName", "queryName", aggfunc="count")
pivot.to_csv(f"ConsMapped-Edit{mismatch}.tsv", sep="\t")
#writer = pd.ExcelWriter('ConsMappingEdit0.xlsx', engine='xlsxwriter')
#pivot.to_excel(writer, sheet_name="Sheet1")
print(tabulate(pivot, headers='keys'))

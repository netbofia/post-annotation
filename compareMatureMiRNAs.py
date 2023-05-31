#! /usr/bin/env python3

# Paper: https://academic.oup.com/bioinformatics/article/24/13/1530/238234
# DOCS: https://bioinf.eva.mpg.de/patman/patman-1.2.html
import argparse
import os
import shutil
import re
import subprocess

parser = argparse.ArgumentParser(prog="compareMatureMiRNA", description="Compare all miRNAs against them selves")
parser.add_argument("-f", "--file", help="TSV with expression")
parser.add_argument("-r", "--report", help="Generate a report", action="store_true")
parser.add_argument("-p", "--plot", help="Generate plots", action="store_true")
args = parser.parse_args()


class SeqComparator:
    def __init__(self, args):
        self.patman = "/home/brunocosta/.Software/patman-1.2.2/patman"
        self.reportFile = "Report.tsv"
        self.file = args.file
        self.args = args

    def append_star_sequences(self):
        print("Not implemented")

    def merge(self):
        tsv = open(self.file, "r")
        sequences = [line.strip().split("\t")[0:2] for line in tsv.readlines()]
        for mismatch in range(0, 3):
            for idx, seq in enumerate(sequences[1::]):
                path = f'vs/{mismatch}/{seq[0]}'
                os.makedirs(path, exist_ok=True)
                fw = open(f'{path}/{seq[0]}.fasta', "w")
                fw.write(f'>{seq[0]}_{seq[1]}\n{seq[0]}')
                fw.flush()
                fw.close()
                fwAll = open(f'{path}/all.fasta', "w")
                fwAll.write("\n".join([f'>{s[0]}_{s[1]}\n{s[0]}' for s in sequences[1::] if s != seq]))
                fwAll.flush()
                fwAll.close()
                print(f'Mismatch:{mismatch} - {idx + 1}/{len(sequences) - 1}')
                patman=self.patman
                subprocess.run(f'{patman} -e {mismatch} -P {path}/{seq[0]}.fasta -D {path}/all.fasta -o {path}/result',
                               shell=True, check=True)
                result = open(f'{path}/result', "r")
                resultSize = os.stat(f'{path}/result')
                if resultSize.st_size == 0:
                    shutil.rmtree(path)

    def report(self):
        report = []
        print("Starting")
        for mismatch in range(0, 3):
            print(f'Mismatch:{mismatch}')
            sequences = os.listdir(f'vs/{mismatch}')
            for sequence in sequences:
                path = f'vs/{mismatch}/{sequence}'
                result = open(f'{path}/result', "r")
                tempScore = {"sequence": sequence, "query": "", "target": "", "len": len(sequence), "queryName": "",
                             "matches": "", "mismatch": mismatch}
                for idx, line in enumerate(result.readlines()):
                    parsedLine = dict([item for item in zip(["target", "query", "start", "end", "orientation", "edits"],
                                                            line.strip().split("\t"))])
                    if re.search("novel", parsedLine["query"]):
                        tempScore["query"] = "novel"
                    cons = re.search("mir[0-9]+_*[0-9]*", parsedLine["query"])
                    if cons:
                        tempScore["query"] = "cons"
                        tempScore["queryName"] = parsedLine["query"][cons.span(0)[0]:cons.span(0)[1]]
                    if re.search("tasi", parsedLine["query"]):
                        tempScore["query"] = "tasi"
                    # targetType
                    if re.search("novel", parsedLine["target"]):
                        if tempScore["target"] == "" or tempScore["target"] == "novel":
                            tempScore["target"] = "novel"
                        elif tempScore["target"] == "cons":
                            tempScore["target"] = "mix"
                    elif re.search("mir", parsedLine["target"]):
                        if tempScore["target"] == "" or tempScore["target"] == "cons":
                            tempScore["target"] = "cons"
                        elif tempScore["target"] == "novel":
                            tempScore["target"] = "mix"
                    elif re.search("tasi", parsedLine["target"]):
                        if tempScore["target"] == "" or tempScore["target"] == "tasi":
                            tempScore["target"] = "tasi"
                        elif tempScore["target"] == "novel" or tempScore["target"] == "cons":
                            tempScore["target"] = "mix"
                    else:
                        tempScore["target"] = parsedLine["target"]
                    tempScore["matches"] = idx + 1
                report.append(tempScore)

        fw = open(self.reportFile, "w")
        fw.write("\t".join(["sequence", "query", "target", "len", "queryName", "matches", "mismatch"]) + "\n")
        for line in report:
            fw.write("\t".join([str(i) for i in line.values()]) + "\n")
        fw.flush()
        fw.close()

    def plot(self):
        import matplotlib.pyplot as plt
        import pandas as pd
        df = pd.read_table(self.reportFile)
        print(df.head())
        df = df.astype({'len': 'int'})
        df = df.astype({'mismatch': 'int'})
        df = df.astype({'matches': 'int'})

        ###Overview
        title = f'Each series is a type of query with the specified number of possible mismatches'
        pt = pd.pivot_table(df, "matches", ["len"], ["mismatch", "query"], "sum")
        consCol = pt.columns.get_level_values(1) == "cons"
        print(consCol)
        pt.columns = ["-".join([str(k) for k in key]) for key in pt.keys()]
        print(pt)
        # pt.index = ["_".join([str(el) for el in i]) for i in pt.index]
        plt.title(title)
        plt.plot(pt)
        plt.legend(list(pt.keys()))
        plt.savefig(title + ".png")
        plt.clf()

        title += " Cons only"
        plt.title(title)
        print(pt.iloc[:, consCol])
        plt.plot(pt.iloc[:, consCol])
        plt.legend(list(pt.iloc[:, consCol].keys()))
        plt.xlabel("Length of Sequence")
        plt.ylabel("Total number of matches with other mature")
        plt.show()
        plt.savefig(f'{title}.png')

        for mismatch in range(0, 3):
            print(f'Mismatch:{mismatch}')
            subDf = df[df["mismatch"] == mismatch]
            # print(pd.crosstab(subDf['len'],subDf['mismatch']))
            plt.plot(pd.crosstab(subDf['len'], subDf['matches']))

            title = f'Mismatch{mismatch}-len vs mismatch.png'
            plt.title(title)

            plt.xlabel("Length of Sequence")
            plt.ylabel("Total number of matches with other mature")
            plt.savefig(title)
            plt.clf()

            ## Cons match count
            title = f'Mismatch{mismatch}-Conserved grouped by family, num of mappings'

            subDf = subDf[subDf["queryName"].notnull()]

            subDf['family'] = subDf["queryName"].str.split("_").str[0]
            print(pd.crosstab(subDf['family'], subDf['matches'], margins=True)["All"].head())
            crossTable = pd.crosstab(subDf['family'], subDf['matches'], margins=True)["All"]
            crossTable = crossTable.drop(["All"])
            # plt.plot(crossTable.values, crossTable.keys())
            # plt.plot(crossTable)
            fig, ax = plt.subplots()
            ax.set_xticklabels(crossTable.keys(), rotation=90)
            plt.plot(crossTable)
            plt.title(title)
            plt.ylabel("Total number of matches with other mature")
            plt.savefig(title + ".png")

            # for length in df["len"].unique():
            #    print(length)
            #    subDf=subDf[subDf["len"]==length]
            #    print(subDf.head())
            #    plt.plot(subDf['mismatch'])
            #    plt.show()

    def run(self):
        if self.args.report:
            self.report()

        if self.args.plot:
            self.plot()

        if self.args.file:
            self.merge()

comp = SeqComparator(args)
comp.run()

##Or export to some other thing

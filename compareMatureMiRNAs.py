#! /usr/bin/env python3

# Paper: https://academic.oup.com/bioinformatics/article/24/13/1530/238234
# DOCS: https://bioinf.eva.mpg.de/patman/patman-1.2.html
import argparse
import os
import shutil
import re
import subprocess
import pandas as pd
from tabulate import tabulate


class SeqComparator:
    def __init__(self, arguments):
        install_patman_path = f'{os.environ["HOME"]}/.Software/patman-1.2.2/patman'
        patman = shutil.which("patman")
        if patman is None:
            patman = install_patman_path
        if patman is not None:
            self.patman = patman
        else:
            print("Patman was not found!\nSource: https://bioinf.eva.mpg.de/patman/\nManual: https://bioinf.eva.mpg.de/patman/patman-1.2.html\nPaper: http://bioinformatics.oxfordjournals.org/cgi/reprint/24/13/1530 ")
        self.reportFile = "reports/Report.tsv"
        self.file = arguments.file
        self.args = arguments
        self.start= arguments.start
        self.end= arguments.end

    def append_star_sequences(self):
        print("Not implemented")

    def map(self):
        tsv = open(self.file, "r")
        sequences = [line.strip().split("\t")[0:2] for line in tsv.readlines()]
        for mismatch in range(self.start, self.end):
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
                patman = self.patman
                subprocess.run(f'{patman} -e {mismatch} -P {path}/{seq[0]}.fasta -D {path}/all.fasta -o {path}/result',
                               shell=True, check=True)
                result = open(f'{path}/result', "r")
                resultSize = os.stat(f'{path}/result')
                if resultSize.st_size == 0:
                    shutil.rmtree(path)
                else:
                    os.remove(f'{path}/all.fasta')
                    os.remove(f'{path}/{seq[0]}.fasta')
        self.merge()

    def merge(self):
        for mismatch in range(self.start, self.end):
            print(f"Processing {mismatch}/{self.end-1}")
            report_file = f'reports/Result-m{mismatch}.tsv'
            subprocess.run(f'printf "target\tquery\tstart\tend\tstrand\tedits\n"  > {report_file}; for seq in '
                           f'vs/{mismatch}/*; do cat $seq/result >> {report_file} ;done', shell=True, check=True)
        print("Finished merge")

    def report(self):
        # Deprecation Warning!
        # Merges all the reports together
        report = []
        print("Starting")
        for mismatch in range(self.start, self.end):
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

    def stats(self):
        for e in range(self.start, self.end):
            sub_report = f"reports/Result-m{e}.tsv"
            if os.path.exists(sub_report):
                df = pd.read_table(sub_report, sep="\t")
                seq_mapped = (len(df["query"].unique()))
                mapped = len(df.drop_duplicates())
                print(f'Edits={e} | {seq_mapped} sequences mapped to a total of {mapped} mappings.')
            else:
                print(f"No report found for edits={e}")
    def tab(self, data):
        print(tabulate(data, keys=data.keys))

    def group(self):
        # Group miRNAs into families for the most permissive edit report
        edits = self.end-1

        sub_report = f"reports/Result-m{edits}.tsv"
        df = pd.read_table(sub_report, sep="\t")

        df = df.astype({'edits': 'int'})
        df = df.astype({'start': 'int'})
        df = df.astype({'end': 'int'})

        df['queryName'] = df['query'].str.split('_').str[1]
        df['querySeq'] = df['query'].str.split('_').str[0]
        df['targetName'] = df['target'].str.split('_').str[1]
        df['targetSeq'] = df['target'].str.split('_').str[0]
        df['family'] = ""

        # novel-star count as novel
        novel_seq = df[df["queryName"].str.startswith("novel")]["querySeq"].unique()
        novel = len(novel_seq)
        family = 0
        print(f"Naming novel miRNAs with edit={edits}")
        # Iterate all sequences of novel miRNAs
        for querySequence in novel_seq:
            for mapping in df[df["querySeq"] == querySequence].itertuples():
                targetSequence = mapping.targetSeq
                if mapping.family == "":
                    df.loc[df["querySeq"] == querySequence, "family"] = family
                    df = self.add_name_to_query(df, targetSequence, family)
                    #print(family)
            if len(df["family"].unique()) > family + 1:
                family += 1

        num_families = len(df['family'].unique())
        print(f'That total of {num_families} families where generated out of {novel} novel sequences')
        df.to_csv(f"reports/Result-e{edits}-families.tsv", sep="\t")
        df.to_excel(f"reports/Result-e{edits}-families.xlsx")
        #count families with cons

    def add_name_to_query(self, dataf, query_sequence, counter):
        # Recursive function that names miRNAs in the <dataf> dataframe
        query_slice = dataf.loc[dataf["querySeq"] == query_sequence]
        for query_row in query_slice.itertuples():
            if query_row.family == "":
                dataf.loc[query_row.Index, "family"] = counter
                target_sequence = query_row.targetSeq
                dataf = self.add_name_to_query(dataf, target_sequence, counter)
        return dataf

    def plot(self):
        # Deprecation Warning!
        # Plots based on report file

        import matplotlib.pyplot as plt
        import pandas as pd
        df = pd.read_table(self.reportFile)
        print(df.head())
        df = df.astype({'len': 'int'})
        df = df.astype({'mismatch': 'int'})
        df = df.astype({'matches': 'int'})

        # Overview
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
        plt.savefig(f'figures/{title}.png')

        for mismatch in range(self.start, self.end):
            print(f'Mismatch:{mismatch}')
            subDf = df[df["mismatch"] == mismatch]
            # print(pd.crosstab(subDf['len'],subDf['mismatch']))
            plt.plot(pd.crosstab(subDf['len'], subDf['matches']))

            title = f'Mismatch{mismatch}-len vs mismatch.png'
            plt.title(title)

            plt.xlabel("Length of Sequence")
            plt.ylabel("Total number of matches with other mature")
            plt.savefig(f"figures/{title}.png")
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
            plt.savefig(f"figures/{title}.png")

            # for length in df["len"].unique():
            #    print(length)
            #    subDf=subDf[subDf["len"]==length]
            #    print(subDf.head())
            #    plt.plot(subDf['mismatch'])
            #    plt.show()


def main():
    parser = argparse.ArgumentParser(prog="compareMatureMiRNA",
                                     description="Compare all miRNAs against themselves")
    parser.add_argument("-f", "--file", help="TSV with expression")
    parser.add_argument("-r", "--report", help="Generate a report (Deprecated use)", action="store_true")
    parser.add_argument("-s", "--start", help="Edits range start", type=int, default=0)
    parser.add_argument("-e", "--end", help="Edits range end (<exclusive>, up to that number!)", type=int,
                        default=1)
    parser.add_argument("-p", "--plot", help="Generate plots (Deprecated use)", action="store_true")
    parser.add_argument("-d", "--stats", help="Show table of simple stats", action="store_true")
    parser.add_argument("-g", "--group", help="Group novel miRNAs into families of iso-forms (requires reports)", action="store_true")

    parser.add_argument("-m", "--merge", help="Merge mapping results", action="store_true")

    args = parser.parse_args()
    seq_comp = SeqComparator(args)

    if args.report:
        seq_comp.report()
    if args.plot:
        seq_comp.plot()
    if args.file:
        args.map()
    if args.merge:
        seq_comp.merge()
    if args.stats:
        seq_comp.stats()
    if args.group:
        seq_comp.group()


if __name__ == "__main__":
    main()


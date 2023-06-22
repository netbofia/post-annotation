#! /usr/bin/env python3

# Paper: https://academic.oup.com/bioinformatics/article/24/13/1530/238234
# DOCS: https://bioinf.eva.mpg.de/patman/patman-1.2.html
import argparse
import os
import shutil
import re
import subprocess
import sys
from pathlib import Path
import pandas as pd
from tabulate import tabulate
import glob2


# TODO make report,figure,tables folder??


class SeqComparator:
    def __init__(self, arguments):
        if arguments is None:
            print("Missing input Error - Constructor initiated without argparser.")
            sys.exit(1)
        else:
            install_patman_path = f'{os.environ["HOME"]}/.Software/patman-1.2.2/patman'
            patman = shutil.which("patman")
            if patman is None:
                patman = install_patman_path
            if patman is not None:
                self.patman = patman
            else:
                print(
                    "Patman was not found!\nSource: https://bioinf.eva.mpg.de/patman/\nManual: https://bioinf.eva.mpg.de/patman/patman-1.2.html\nPaper: http://bioinformatics.oxfordjournals.org/cgi/reprint/24/13/1530 ")
            self.reportFile = "reports/Report.tsv"
            self.file = arguments.file
            self.args = arguments
            self.start = arguments.start
            self.end = arguments.end
            self.mircat = arguments.mircat
            self.families = None

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
            print(f"Processing {mismatch}/{self.end - 1}")
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
        if not self.families:
            print("Not implemented yet!")

    def tab(self, dataframe):
        # https://www.geeksforgeeks.org/display-the-pandas-dataframe-in-table-style/
        print(tabulate(dataframe, headers='keys', tablefmt='psql'))

    def group(self):
        # Group miRNAs into families for the most permissive edit report
        edits = self.end - 1
        print(f'Running with edits={edits}')
        sub_report = f"reports/Result-m{edits}.tsv"
        if not os.path.exists(sub_report):
            print(f"Error {sub_report} missing. Run merge first")
            sys.exit(1)
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
                    # print(family)
            if len(df["family"].unique()) > family + 1:
                family += 1

        num_families = len(df['family'].unique())
        print(f'That total of {num_families} families where generated out of {novel} novel sequences')
        df.to_csv(f"reports/Result-e{edits}-families.tsv", sep="\t")
        df.to_excel(f"reports/Result-e{edits}-families.xlsx")
        print("Report exported to reports")

        self.identify_families_mapping_with_cons()

        # TODO  check if that are missing sequences !Attention these are only the miRNAs that mapped

    def join_families_with_expression_matrix(self, df_families):
        # Adds families to expression matrix
        if not os.path.exists(self.file):
            print("No file with the expression matrix")
            sys.exit(2)
        dfAll = pd.read_table(self.file, sep="\t")
        dfAll.index = dfAll['sequence']
        dfAll.loc[:, 'mapped'] = False
        dfAll.loc[:, 'MapsAsCons'] = False
        dfAll.loc[dfAll['sequence'].isin(df_families.querySeq.unique()), 'mapped'] = True
        dfAll.loc[dfAll['sequence'].isin(df_families[df_families['MapsAsCons']].querySeq.unique()), "MapsAsCons"] = True

        mapped = df_families.loc[
            (df_families.queryName.str.startswith("novel")), ['querySeq', "family", "MapsAsCons"]].drop_duplicates()
        # TODO are there miRNAS with Maps true and false?
        for row in mapped.itertuples():
            if row.MapsAsCons:
                # Sets family as conserved name
                # Might cause adding the last conservedTarget
                dfAll.loc[row.querySeq, "family"] = df_families.loc[(df_families['querySeq'] == row.querySeq) &
                                                                    df_families['targetName'].str.startswith(
                                                                        "mir"), "targetName"].unique()[0]
            else:
                dfAll.loc[row.querySeq, 'family'] = int(row.family)

        novel_that_are_conserved = df_families.loc[
            (df_families.targetName.str.startswith("novel")) & df_families.MapsAsCons, 'targetSeq'].unique()
        for false_novel in novel_that_are_conserved:
            cons_names = df_families.loc[(df_families.targetSeq == false_novel) & df_families.queryName.str.startswith(
                "mir"), "queryName"].unique()
            families = df_families.loc[df_families.querySeq == false_novel, "family"].unique()
            #print(f'Families:' + str(families) + ' Conserved names:' + " ".join(cons_names))
            dfAll.loc[dfAll["sequence"] == false_novel, "family"] = cons_names[0]
            if len(families) > 0:
                # Sets that family to the conserved one for all miRNAs that have that family
                for family in families:
                    dfAll.loc[dfAll["family"] == family, "family"] = cons_names[0]

        novel_counter = df_families[df_families.family != ""].family.astype("float").max().astype("int")
        for row in dfAll[dfAll['name'].str.startswith("novel")].itertuples():
            if pd.isna(row.family):
                novel_counter += 1
                dfAll.loc[row.sequence, 'family'] = novel_counter

        # From 31372 to 25077
        return dfAll

    def identify_families_mapping_with_cons(self):
        # Adds new column to file
        edits = self.end - 1
        families_report = f"reports/Result-e{edits}-families.tsv"
        if not os.path.exists(families_report):
            print(f"Missing file: {families_report}")
            sys.exit(2)
        df = pd.read_table(families_report, sep="\t")

        # Identifies miRNAs that can be considered conserved iso-forms
        df["MapsAsCons"] = (df["targetName"].str.startswith("mir")) & (df["queryName"].str.startswith("novel"))
        df.loc[
            (df["queryName"].str.startswith("mir")) & (df["targetName"].str.startswith("novel")), "MapsAsCons"] = True
        df.to_csv(families_report, sep="\t")

    def group_mature_star_into_same_family(self, retry=0):
        self.identify_families_mapping_with_cons()
        # Adds new column to file
        edits = self.end - 1
        families_report = f"reports/Result-e{edits}-families.tsv"

        if not os.path.exists(families_report):
            print(f"Missing file: {families_report}")
            sys.exit(2)
        df_families = pd.read_table(families_report, sep="\t")
        index_columns = [col for col in df_families.columns if col.startswith("Unnamed: 0")]
        df_families.drop(index_columns, axis=1, inplace=True)

        # Stop retry
        if df_families.keys().to_list().count("MapsAsCons") == 0 and retry < 2:
            self.identify_families_mapping_with_cons()
            retry += 1
            self.group_mature_star_into_same_family(retry)

        dfAll = self.join_families_with_expression_matrix(df_families)

        mircat_folder = self.mircat
        mircat_folder_path = Path(mircat_folder)
        mircat_files = None

        if mircat_folder_path.is_dir():
            # Remove fixed in production
            mircat_files = glob2.glob(mircat_folder + '/*_mirbase_noncons_output_filteredfixed.csv')
        else:
            print(f"Directory doesn't exist! {mircat_folder}")
            sys.exit(2)

        header = "Chromosome,Start,End,Orientation,Abundance,Sequence,sRNA length,# Genomic Hits,Hairpin Length," \
                 "Hairpin % G/C content,Minimum Free Energy,Adjusted MFE,miRNAstar,miRBaseID,pVal".split(",")
        df_list = [pd.read_csv(file, skiprows=1, header=None, sep=",", on_bad_lines="warn") for file in mircat_files]
        df = pd.concat(df_list)
        df.columns = header
        #####################################################3
        print(df.shape)
        # Should be 75246-9 headers =75237 checks out fine
        df = df.drop_duplicates()
        # 75237 --> 72506
        print(df.shape)

        starDf = pd.DataFrame()
        for starRow in df[df['miRNAstar'] != "NO"].itertuples():
            for miRNAstar in starRow.miRNAstar.strip().split(" "):
                starDf.loc[len(starDf), ["mature", "star", "starAbundance", "chromosome", "start", "end"]] = [
                    starRow.Sequence, miRNAstar.split("(")[0], miRNAstar.split(")")[0].split("(")[1],
                    starRow.Chromosome, starRow.Start, starRow.End]

        # 3336 rows
        # 1169 unique
        print(starDf.duplicated().value_counts())
        # Unique mature : star = 1988
        number_of_unique_stars = len(starDf.star.unique())
        print(f"Number of unique star sequences identified: {number_of_unique_stars}")

        unique_families = len(dfAll.family.unique())
        list_of_stars = starDf['star'].unique()



        ## Identify missing conserved-star
        dfAll['Star'] = dfAll.sequence.isin(list_of_stars)

        for star in list_of_stars:
            family_star = self.get_family_from_families_df(dfAll, star)
            for mature in starDf.loc[starDf['star'] == star, "mature"].unique():
                family_mature = self.get_family_from_families_df(dfAll, mature)
                if family_star != family_mature:
                    if pd.isna(family_star):
                        if pd.isna(family_mature):
                            print(star)
                        else:
                            dfAll.loc[dfAll['family'] == family_mature, "family"] = family_mature
                    else:
                        dfAll.loc[dfAll['family'] == family_mature, "family"] = family_star
        unique_families_reduced_by_precursor_based_pairing_of_ends = len(dfAll.family.unique())
        print(
            f"Reduced the number of families from {unique_families} to {unique_families_reduced_by_precursor_based_pairing_of_ends}")
        dfAll = self.remove_gaps_in_families(dfAll)
        dfAll.to_excel("tables/All_seqs-with-families.xlsx")

    def get_family_from_families_df(self, df, sequence):
        family = df[df['sequence'] == sequence].family.unique()
        if len(family) != 1:
            print(f"Error - Grouping generated multiple families ({len(family)}) for a miRNA, cannot proceed! "
                  f"Sequence:{sequence}")
            sys.exit(2)
        return family[0]

    def remove_gaps_in_families(self, df=None):
        if df is None:
            edits = self.end - 1
            exp_mat_path = f"tables/All_seqs-with-families.tsv"
            if not os.path.exists(exp_mat_path):
                print(f"Missing file: {exp_mat_path}")
                sys.exit(2)
            df = pd.read_table(exp_mat_path, sep="\t", dtype={"family": str})
            #df.fillna(-99999,inplace=True)
            df.family.astype(str)

        #Fill gaps in numbering by pull the higher numbers into the gaps.
        families = df.loc[(~df['family'].astype("str").str.startswith("mir")) & df["name"].str.startswith("novel"), 'family'].unique().astype(int)
        families.sort()
        families = [list(item) for item in zip(families, [-1]*len(families))]
        rm_idx = len(families)-1
        for idx, family in enumerate(families):
            if idx == 0:
                if family[0] != 0:
                    families[rm_idx][1] = 0
                    while families[rm_idx][1] != family[0]-1:
                        rm_idx -= 1
                        last_family_added = families[rm_idx+1][1]+1
                        families[rm_idx][1] = last_family_added+1
                    rm_idx -= 1
                else:
                    family[1] = 0
            else:
                #If is increment to previous family or last removed family
                previous_family = families[idx-1][1]
                last_removed_family = -999
                # Exception for first cycle
                if len(families) == rm_idx+1:
                    last_removed_family = families[rm_idx][1]
                else:
                    last_removed_family = families[rm_idx+1][1]
                if family[0]-1 == previous_family or family[0]-1 == last_removed_family:
                    if family[1] == -1:
                        family[1] = family[0]
                    else:
                        if family[1] > family[0]:
                            family[1] = family[0]

                else:
                    if family[1] != -1:
                        if family[1] > family[0]:
                            family[1] = family[0]
                    else:
                        families[rm_idx][1] = previous_family+1
                        while families[rm_idx][1] != family[0] - 1 and rm_idx > -1:
                            rm_idx -= 1
                            last_family_added = families[rm_idx + 1][1]
                            families[rm_idx][1] = last_family_added + 1
                        rm_idx -= 1
                        family[1] = family[0]
        for family in families:
            if family[0] != family[1]:
                #Lower the family number to fill the gaps
                df.loc[df['family'] == family[0], "family"] = family[1]
        df.to_csv("tables/all_seq-Renumbered.tsv", sep="\t")
        return df

    def generate_full_names(self, df=None):
        # First run use
        if df is None:
            exp_mat_path = f"tables/all_seq-Renumbered.tsv"
            if not os.path.exists(exp_mat_path):
                print(f"Missing file: {exp_mat_path}")
                sys.exit(2)
            df = pd.read_table(exp_mat_path, sep="\t", dtype={"family": str})
            df.drop("sequence.1", axis=1, inplace=True)
            df.drop("Unnamed: 0", axis=1, inplace=True)
            df.index = df.sequence
            df.drop("sequence", axis=1, inplace=True)
            df.dropna(how="all", inplace=True)
            df = df[df.index.notnull()]

        df["FullName"] = ""
        df.loc[df['name'].str.startswith("mir"), "FullName"] = \
            "qsu-miRmpc"+df.loc[df['name'].str.startswith("mir"), "name"].str.split("mir").str[1]

        # ## Adds in
        # TODO Information on orientation might help
        df.loc[(df['name'].str.startswith("mir")) & (df["Star"] == False), "FullName"] = df.loc[
            df['name'].str.startswith("mir"), "FullName"]+"-5p"
        df.loc[(df['name'].str.startswith("mir")) & (df["Star"] == True), "FullName"] = df.loc[
            df['name'].str.startswith("mir") & df["Star"], "FullName"] + "-3p"


        #novel-merged-in-conserved-group
        df.loc[df['name'].str.startswith("novel") & df["family"].astype(str).str.startswith("mir"), "FullName"] = \
            "qsu-miRmpc" + df.loc[df["family"].astype(str).str.startswith("mir"), "family"].str.split("mir").str[1]

        df.loc[(df['name'].str.startswith("novel")) & (df["Star"] == False), "FullName"] = \
            df.loc[df['family'].astype(str).str.startswith("mir"), "FullName"] + "-5p"
        df.loc[(df['name'].str.startswith('novel')) & (df["Star"] == True), "FullName"] = \
            df.loc[df['family'].astype(str).str.startswith("mir") & df["Star"], "FullName"] + "-3p"

        #print(df.loc[df['name'].str.startswith("novel"), 'FullName'].unique())

        # novel
        df.loc[df['name'].str.startswith("novel") & (~df["family"].astype(str).str.startswith("mir")), "FullName"] = \
            "qsu-miRmp" + df.loc[~df["family"].astype(str).str.startswith("mir"), "family"]

        df.loc[(df['name'].str.startswith("novel")) & (df["Star"] == False) & (~df["family"].astype(str).str.startswith("mir")), "FullName"] = \
            df.loc[~df['family'].astype(str).str.startswith("mir"), "FullName"] + "-5p"
        df.loc[(df['name'].str.startswith('novel')) & (df["Star"] == True) & (~df["family"].astype(str).str.startswith("mir")), "FullName"] = \
            df.loc[~df['family'].astype(str).str.startswith("mir") & df["Star"], "FullName"] + "-3p"

        df.to_csv("tables/all_seq-fullname.tsv", sep="\t")
        df.to_excel("tables/all_seq-fullname.xlsx")

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
    parser.add_argument("-n", "--mircat", help="path to miRCat directory")
    parser.add_argument("-p", "--plot", help="Generate plots (Deprecated use)", action="store_true")
    parser.add_argument("-d", "--stats", help="Show table of simple stats", action="store_true")
    parser.add_argument("-g", "--group", help="Group novel miRNAs into families of iso-forms (requires reports)",
                        action="store_true")
    parser.add_argument("-m", "--merge", help="Merge mapping results", action="store_true")
    parser.add_argument("--rm", help="Remove gaps in family sequential naming", action="store_true")
    parser.add_argument("--fullname", help="Adds a column with the full names", action="store_true")

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
    if args.mircat:
        seq_comp.group_mature_star_into_same_family()
    if args.rm:
        seq_comp.remove_gaps_in_families()
    if args.fullname:
        seq_comp.generate_full_names()


if __name__ == "__main__":
    main()

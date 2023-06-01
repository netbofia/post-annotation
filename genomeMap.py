#!/usr/bin/env python3

######################################



import argparse
import glob2
from pathlib import Path
import json
import re


parser = argparse.ArgumentParser(prog="genomeMapper", description="Generates a feature map from the miRCat hairpin files" )
parser.add_argument("-f", "--hairpinFolder")
args = parser.parse_args()
hairpinFolder = args.hairpinFolder
hairpinFolderPath = Path(hairpinFolder)
hairpinFiles = None
if hairpinFolderPath.is_dir():
    hairpinFiles = glob2.glob(hairpinFolder+'/*_miRNA_hairpins_filtered.txt')


def parse_hairpin_files(hairpin_files):
    #Reads the fasta files and parsed the header; sequence and hairpin, as named dictionary
    result = list()
    if len(hairpin_files) > 0:
        for fasta in hairpin_files:
            open_fasta = open(fasta, "r")
            while open_fasta.readable():
                current_line = open_fasta.readline()
                if len(current_line) == 0:
                    break
                if current_line.startswith(">"):
                    result += [{"header": current_line.strip().split("_", 1)[1].rsplit("/", 1), "sequence": open_fasta.readline().strip(), "structure": open_fasta.readline().strip()}]
    return result


lines = parse_hairpin_files(hairpinFiles)
##lines+=[ line.strip().split("_",1)[1].rsplit("/",1) for idx, line in enumerate(open(fasta,"r").readlines()) if re.match("^[().-<>]",line) ]


def locate_mature(mature_start, structure):
    match = re.match("[-<>]*", structure)
    start_c = int(mature_start) - match.start()
    end_c = int(start_c)+len(structure)
    return [start_c, end_c]

print(lines[1])
scaffolds = dict()

#
for line in lines:
    identifier = line["header"][0].split(" ", 1)[0]
    miCoordinates = line["header"][1].split("-")

    miStart = miCoordinates[0]
    miEnd = miCoordinates[1]
    coordinates = locate_mature(miStart, line["structure"])

    start = coordinates[0]
    end = coordinates[1]
    try:
        scaffolds[identifier][start]={"id": "", "start": start, "end": end, "type": "precursor"}
    except KeyError:
        scaffoldFeaturesSortedByStart=dict()
        scaffoldFeaturesSortedByStart[start]={"id": "", "start": coordinates[0], "end": coordinates[1], "type": "precursor"}
        scaffolds[identifier]=scaffoldFeaturesSortedByStart

#"size": 248956422,
#"bands": [
#            {"id": "p11.1", "start": 121700001, "end": 123400000, "type": "acen"},
#            {"id": "p11.2", "start": 120400001, "end": 121700000, "type": "gneg"},
#            {"id": "p12", "start": 117200001, "end": 120400000, "type": "gpos50"},
#            {"id": "p13.1", "start": 115500001, "end": 117200000, "type": "gneg"},
#
#for each scaffold

scaffoldName="NW_019827908.1"
#print(scaffolds["NW_019827908.1"]["1388709"])


scaffold=dict()
scaffold[scaffoldName]={"size":1,"bands":[]}

for idx,i in enumerate(scaffolds[scaffoldName]):
    scaffolds[scaffoldName][i]['id']="p"+str(idx)
    scaffold[scaffoldName]["bands"].append(scaffolds[scaffoldName][i])
scaffold[scaffoldName]["size"]=scaffold[scaffoldName]["bands"][-1]['end']
print(scaffold)

fw = open("/home/brunocosta/git/sRNA-Portal-workflow/public/javascripts/genomes/qsuber.json","w")
fw.write(json.dumps(scaffold,indent=4))
fw.flush()
fw.close()



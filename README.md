# Post-annotation
MiRPursuit Post Annotation tools

This should run after miRPursuit annotation pipeline is finished. 

Dependencies:
- python3
- pandas
- matplotlib.pyplot

Requires the mircat folder and the all-seq.tsv from count:
- counts/all-seq.tsv
- data/mircat

## Usage

```bash
python3 compareMatureMiRNA.py
```

usage: compareMatureMiRNA [-h] [-f FILE] [-r] [-s START] [-e END] [-p]
                          [-m]

Compare all miRNAs against themselves

options:
-  -h, --help              Show this help message and exit
-  -f FILE, --file FILE    TSV with expression [First step]
-  -r, --report            Generate a report
-  -s START, --start START
                           Edits range start
-  -e END, --end END       Edits range end (<exclusive>, up to that number!)
-  -p, --plot              Generate plots
-  -n PATH, --mircat PATH  Path to miRCat directory
-  -m, --merge             Merge mapping results
-  --rm                    Remove gaps in family sequential naming
-  --fullname              Adds a column with the full names


### Mapping
With -f it does the mapping of each miRNA against the others. The default range is 0 to 1 (edits=0) and does the merging afterward.

### Merge
You can do the merge step without running the mapping. This generates a Report of results for each edit value

### Group
Group miRNA together into one family

#### Group by precursor
If a miRNA is found on a precursor with a star sequence. The star sequence and the mature miRNA are bundled together into the same family

### Remove gaps in family naming
Due to further merging some families disappear and leave gaps in the sequencial numbering. This process gets the higher numbers into those gaps. 

### Full name
Adds a column with the full name of the novel miRNAs



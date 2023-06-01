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
-  -h, --help            show this help message and exit
-  -f FILE, --file FILE  TSV with expression [First step]
-  -r, --report          Generate a report
-  -s START, --start START
                        Edits range start
-  -e END, --end END     Edits range end (<exclusive>, up to that
                        number!)
-  -p, --plot            Generate plots
-  -m, --merge           Merge mapping results

### Mapping
With -f it does the mapping of each miRNA against the others. The default range is 0 to 1 (edits=0) and does the merging afterwards.

### Merge
You can do the merge step without running the mapping. This generates a Report of results for each edit value

### 



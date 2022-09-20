# Map Residue Features from UniProt

The downloadUniprotDomains.py script outputs a folder "Uniprot_Domains" with a subfolder named after the date. Output files are added to the subfolder.

UniProt annotations are filtered for human proteome reviewed enteries. For more downloadable annotation column options see <https://www.uniprot.org/help/uniprotkb_column_names>. To search a query string using a web browser, remove 'format', 'compress', and 'force'.

$ python downloadUniprotDomains.py  # To run script

### Functions:

```python
def queryUKB(SAVEFULL, SAVESUMMARY, FOLDERNAME)
    def parseRegion(row, splitKey)
    def parseSite(row, splitKey)

args for queryUKB():
1. full feature save filename contains unparsed data from uniprot api (full_human_features.txt) 
2. summary of feature positions filename contains parsed data and same line count as full feature download (features_summary.txt)
3. folder name for output file dir/ (Uniprot_Domains)

def compareUKBsequences(SAVESUMMARY, REFSEQ, QCoutput)

args for sequence comparison between queried data vs reference for cpdaa positions QC:
1. CCDSfromfasta_REFERENCE.csv
2. features_summary_QCd_sequences.csv

Run log:
2.12.22 line count of saved output files from run:
    20376 features_summary.txt
    20376 full_human_features.txt
    18417 features_summary_QCd_sequences.csv 
    
```

# Map Residue Features from UniProt

The downloadUniprotDomains.py script outputs a folder "Uniprot_Domains" with a subfolder named after the date. Output files are added to the subfolder.

UniProt annotations are filtered for human proteome reviewed enteries. For more downloadable annotation column options see <https://www.uniprot.org/help/uniprotkb_column_names>. To search a query string using a web browser, remove 'format', 'compress', and 'force'.

```bash
python downloadUniprotDomains.py  # To run script
```

### Functions:

```python
def queryUKB(SAVEFULL, SAVESUMMARY, FOLDERNAME)
    def parseRegion(row, splitKey)
    def parseSite(row, splitKey)
```

args for queryUKB():
1. full feature save filename contains unparsed data from uniprot api (full_human_features.txt) 
2. summary of feature positions filename contains parsed data and same line count as full feature download (features_summary.txt)
3. folder name for output file dir/ (Uniprot_Domains)

```python
def compareUKBsequences(SAVESUMMARY, REFSEQ, QCoutput)
```

args for sequence comparison between queried data vs reference for cpdaa positions QC:
1. CCDSfromfasta_REFERENCE.csv
2. features_summary_QCd_sequences.csv

---

Run log:
2.12.22 line count of saved output files from run:
    20376 features_summary.txt
    20376 full_human_features.txt
    18417 features_summary_QCd_sequences.csv 
    
notes from query link:
    2021_04 version of ukb was downloaded on feb. 6 using browser share button link copy and edit for this script.

----Only enteries w/ 3d structure, n=7,376------
https://www.uniprot.org/uniprot/?query=reviewed%3Ayes%20organism%3A%22Homo%20sapiens%20(Human)%20%5B9606%5D%22%20proteome%3Aup000005640%20keyword%3A%223D-structure%20%5BKW-0002%5D%22&columns=

------All reviewed human proteome enteries, n=20,360------------
https://www.uniprot.org/uniprot/?query=&fil=proteome%3AUP000005640%20AND%20reviewed%3Ayes%20AND%20organism%3A%22Homo%20sapiens%20(Human)%20%5B9606%5D%22&columns=

------------------------------------------------------------------------------
download link failed w/ line..
"&fil=reviewed:yes%20organism:%22Homo%20sapiens%20(Human)%20[9606]%22" \

CONVERTING URL to string that works with script:
    %3A is :
    %5B is [
    %5D is ]
    %2C is ,
    %22 is "

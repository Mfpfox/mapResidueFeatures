# Map Residue Features from UniProt

Project has two main parts:
1. downloadUniprotDomains.py
2. isposition_atSite_inDomain.py

## [1] downloadUniprotDomains.py
* script outputs a folder "Uniprot_Domains" with a subfolder named after the date, all output files are added to this subfolder.
* UniProt annotations are filtered for human proteome reviewed enteries. For more downloadable annotation column options see <https://www.uniprot.org/help/uniprotkb_column_names>. To search a query string using a web browser, remove 'format', 'compress', and 'force'.

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



```
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
```



##  [2] isposition_atSite_inDomain.py

### Functions:
```python
addSiteAnnotations(infile, head, outfile, siteCol, window)
addRegionAnnotations(infile, head, outfile, siteCol, window)
ispositionAtSiteInDomain()
```


```
--------------updated 7/21/22--------------
cleaned comments updated descriptions,
added () to proximal boundary if statements in both functions,
added low_memory=False to import of position level reference

------------------input files------------------
features_summary_QCd_IdenticalSeq_18231proteins.csv
    functions use input file because it only includes annotations for proteins with identical sequence between reference of project and 2/12/2022 uniprotDomainDownload
    * made in isposition_atSite.ipynb

REF_posID_level_18827cpdaa_4535UKBIDs.csv
    (all CpDAA positions)

codonMaxMeanScored_mergedDescribeProt_1231_MendelCpD.csv v1 run
    ($$$ all aa positions in mendelCpD)

codonMaxMeanScored_mergedDescribeProt_3853_Mendelian v1 run

--------------output files------------------
[input name]_SITES_ANNOTATED_window_[windowSize].csv

[input name]_REGIONS_ANNOTATED_window_[windowSize].csv

    ** Notes on output: each run produces two unneeded files that have extension "SITESmerged.csv" or "REGIONSmerged.csv" , deleting these output files


```

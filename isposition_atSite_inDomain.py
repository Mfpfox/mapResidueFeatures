# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
sys.path.append("/Users/mariapalafox/Desktop/TOOLBOXPY")
from all_funx import *

"""
--------------updated 7/21/22--------------
cleaned comments updated descriptions,
added () to proximal boundary if statements in both functions,
added low_memory=False to import of position level reference

------------------functions------------------
addSiteAnnotations(infile, head, outfile, siteCol, window)
addRegionAnnotations(infile, head, outfile, siteCol, window)
ispositionAtSiteInDomain()

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

"""

def addSiteAnnotations(infile, head, outfile, siteCol, window):
    f = open(infile)
    df = list(csv.reader(f))
    with open(outfile, 'w') as out:
            csvWriter = csv.writer(out)
            csvWriter.writerow(head)
    # loop over rows skip header
    for row in df[1:]:
        newrow = []
        newrow.append(row[0]) # posID
        pos = int(row[1])
        newrow.append(pos)
        annoRow = row[2:] # select site annotation columns, skip posID and pos
        for i, v in enumerate(annoRow):
            flounder = False # FLAG EXACT
            flounder2 = False # FLAG PROXIMAL
            siteType = siteCol[i]
            # nan annotation value case
            if (len(v) < 1) | (v == ''):
                newrow.append(v) 
                newrow.append(flounder) # exact
                newrow.append(flounder2) # proximal
            # not nan annotation value case
            if len(v) >= 1: 
                if siteType != 'DisulfideBond':
                    vsplit = v.split(";")
                    # EXACT
                    for exac in vsplit:
                        try:
                            posCheck = int(exac.strip(" "))
                        except:
                            print("error posCheck:", posCheck)
                            posCheck = -100
                        # check exact match
                        if pos == posCheck:
                            flounder = True
                            break
                    # PROXIMAL
                    for prox in vsplit:
                        try:
                            posCheck = int(prox.strip(" "))
                        except:
                            print("error posCheck:", posCheck)
                            posCheck = -100
                        # check proximal match
                        if (pos >= posCheck-6 and pos < posCheck) | (pos <= posCheck+6 and pos > posCheck):
                            flounder2 = True
                            break
                    newrow.append(v)
                    newrow.append(flounder)
                    newrow.append(flounder2)
                # disulfide case
                if siteType == 'DisulfideBond':
                    vsplit = v.split(";")
                    for exac in vsplit:
                        posCheck = exac.strip(" ")
                        posCheck = posCheck.split("..")
                        for disul in posCheck:
                            try:
                                di = int(disul)
                            except: 
                                print("error disul:", disul)
                                di = -100
                            # check exact match
                            if pos == di:
                                flounder = True
                                break
                    for prox in vsplit:
                        posCheck = prox.strip(" ")
                        posCheck = posCheck.split("..")
                        for pp in posCheck:
                            try:
                                pc = int(pp)
                            except:
                                print("error pp:", pp)
                                pc = -100
                            # check proximal match
                            if (pos >= (pc - window) and pos < pc) | (pos <= (pc + window) and pos > pc):
                                flounder2 = True
                                break
                    newrow.append(v)
                    newrow.append(flounder)
                    newrow.append(flounder2)
        with open(outfile, 'a') as out:
            csvWriter = csv.writer(out)
            csvWriter.writerow(newrow)


def addRegionAnnotations(infile, head, outfile, siteCol, window):
    """
    -------what code does--------------
    exact overlap example for posID C5, annotations range 3-5 or range 2-6 (inclusive range boundaries)
    proximal overlap example for posID C5, range 2-4, but not 2-5 (mutally exclusive from exact overlap boolean column)
    proximal defined as:
        # -[window integer] from lower bound of range and 
        # +[window integer] for upper bound of range
    case where exact is TRUE and proximal is TRUE example: 
        # ranges 3-5;10-12. position C5 is exact for range 3-5 and proximal for range 10-12 (10 is lower bound, and positions for proximal as TRUE include 9,8,7,6,5 positions)
    """
    f = open(infile)
    df = list(csv.reader(f))
    with open(outfile, 'w') as out:
            csvWriter = csv.writer(out)
            csvWriter.writerow(head)
    for row in df[1:]:
        newrow = []
        newrow.append(row[0]) # posID
        pos = int(row[1]) # aa position to check in domain range
        newrow.append(pos)
        annoRow = row[2:] # select region columns, skip posID and pos cols
        for i, v in enumerate(annoRow):
            flounder = False # FLAG EXACT
            flounder2 = False # FLAG PROXIMAL
            siteType = siteCol[i]
            # nan annotation value case, ADD EMPTY STRING TO FILTER
            if (len(v) < 1) | (v == ''):
                newrow.append(v) 
                newrow.append(flounder) # exact 
                newrow.append(flounder2) # proximal
            # not nan annotation value case
            if len(v) >= 1:
                vsplit = v.split(";")
                # EXACT
                for exac in vsplit:
                    posCheck = exac.strip(" ")
                    posCheck = posCheck.split("-")
                    if len(posCheck) == 2:
                        b1 = int(posCheck[0])
                        b2 = int(posCheck[1])
                        r = [*range(b1, b2+1, 1)] # include 22,24 is 22,23,24 range
                        if pos in r:
                            # check exact match 
                            flounder = True
                            break      
                    if len(posCheck) == 1:
                        try:
                            b1 = int(posCheck[0])
                        except:
                            print("error posCheck:", posCheck)
                            b1=-100
                        if pos == b1:
                            # check exact match 
                            flounder = True
                            break
                # PROXIMAL
                for prox in vsplit:
                    posCheck = prox.strip(" ")
                    posCheck = posCheck.split("-")
                    try:
                        b1 = int(posCheck[0])
                        b2 = int(posCheck[1])
                    except:
                        print("error posCheck:", posCheck)
                        b1=-100
                        b2=-100
                    # check proximal match
                    if (pos >= (b1 - window) and pos < b1) | (pos <= (b2 + window) and pos > b2):
                        flounder2 = True
                        break
                newrow.append(v)
                newrow.append(flounder)
                newrow.append(flounder2)
        with open(outfile, 'a') as out:
            csvWriter = csv.writer(out)
            csvWriter.writerow(newrow)

def ispositionAtSiteInDomain():
    """
    -------what code does--------------
    imports and merges on position index ref ukbIDs, so only looking at proteins with annotations from ukb
    adds an integer pos column from split posID
    splits merged df into site and region dfs, drops columns other than posID, pos, site|regoin annotions
    *makes new header column with exact and proximal colnames added for site and region dfs
    ---------variables---------
    """
    # import QC ukb annotations ref (made in isposition_atSite_inDomain.ipynb)
    # only includes annotations for proteins with identical sequence
    UKBINFILE = "features_summary_QCd_IdenticalSeq_18231proteins.csv"

    # import posID level ref to merge with annotations
    dirr = "/Users/mariapalafox/Desktop/BRIDGE/disorder/dbNSFPmapping/MAXMEAN_SCORING/"
    # dir for REF 18k: "/Users/mariapalafox/Desktop/BRIDGE/DISORDER/ChemoproteomicData/"

    POSINFILE = "codonMaxMeanScored_mergedDescribeProt_3853_Mendelian.csv"
    # POSINFILE = "codonMaxMeanScored_mergedDescribeProt_1231_MendelCpD.csv"
    # POSINFILE = "REF_posID_level_18827cpdaa_4535UKBIDs.csv"

     # '.csv' must be included
    saveSite = POSINFILE.replace(".csv", "_SITESmerged.csv")
    saveRegion = POSINFILE.replace(".csv", "_REGIONSmerged.csv")
    # outfile naming for site and region
    outSite = saveSite.replace("merged.csv", "_ANNOTATED_window_3.csv")
    outRegion = saveRegion.replace("merged.csv", "_ANNOTATED_window_3.csv")
    # change number in header string to match PROXIMAL size set
    window = 3
    proxHEADER = "_Proximal3"
    # name of column with UKBIDs for both dataframes
    IDname = "UKBID"
    posIDname = "posID"
    #----------------------------------------
    positions = pd.read_csv(dirr + POSINFILE, low_memory=False)
    ukb = pd.read_csv(UKBINFILE)
    print("position level reference head:")
    print(positions.head())
    print(positions.shape)
    uniqueCount(positions, posIDname)
    uniqueCount2(positions, IDname)
    print()
    print("UKB identical sequence QC head:")
    print(ukb.head())
    print(ukb.shape)
    uniqueCount(ukb, IDname)
    print()
    # merge ukb annotations on pos level df
    merged1 = pd.merge(positions, ukb, on=[IDname], how="inner")
    # add position integer col from split posID
    merged1= getAAposition(merged1, posIDname, 'pos')
    print("merged df head:")
    print(merged1.head())
    print(merged1.shape)
    uniqueCount(merged1, posIDname)
    uniqueCount(merged1, IDname)
    print()
    # make header cols for merged df's specific to site and region annos.:
    basicHead = [posIDname, 'pos']

    # sites
    sitesCol= ['SubstrateBindingSite', 'ActiveSite', 'MetalBindSite',
           'ModifiedResidue', 'Lipidation', 'Glycosylation', 'Mutagenesis',
           'DisulfideBond'] 
    sitesDF = merged1[basicHead + sitesCol].copy()
    sitesDF.to_csv(saveSite, index=False) # save site and region df's

    sitesAddCol = [] # make exact and proximal header additions SITES
    for coli in sitesCol:
        exx = coli + '_Exact'
        prx = coli + proxHEADER
        sitesAddCol.append(coli)
        sitesAddCol.append(exx)
        sitesAddCol.append(prx)
    # regions
    regionCol = [
            'NucleotideBindRegion', 'DNAbindRegion', 'ZincFinger',
           'Motif', 'DisorderedRegion', 'Domain', 'BetaStrand', 'Helix', 
            'Turn']
    regionDF = merged1[basicHead + regionCol].copy()
    regionDF.to_csv(saveRegion, index=False) # save site and region df's

    regAddCol = [] # make exact and proximal header additions REGIONS
    for coli in regionCol:
        exx = coli + '_Exact'
        prx = coli + proxHEADER # variable
        regAddCol.append(coli)
        regAddCol.append(exx)
        regAddCol.append(prx)

    # new headers
    newheaderSites = basicHead + sitesAddCol
    newheaderRegion = basicHead + regAddCol

    print("calling addSiteAnnotations()...")
    addSiteAnnotations(saveSite, newheaderSites, outSite, sitesCol, window)

    print("calling addRegionAnnotations()...")
    addRegionAnnotations(saveRegion, newheaderRegion, outRegion, regionCol, window)

ispositionAtSiteInDomain()
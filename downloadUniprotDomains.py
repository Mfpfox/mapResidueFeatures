from requests import get
from time import strftime
from os import path, makedirs
from csv import reader
import argparse
from all_funx import *

# arguments added here
# parser = argparse.ArgumentParser(description=" ")
# parser.add_argument("queryresults", help=" ")
# parser.add_argument("summaryresults", help=" ")
# parser.add_argument("cpdaa", help="Residue-level chemoproteomic .csv filename")
# args = parser.parse_args()
# FULLFILE = args.queryresults
# SUMFILE = args.summaryresults
# cpdaa = args.cpdaa

FULLFILE = "full_human_features.txt" # ukb api dump
SUMFILE = "features_summary.txt" # parsed and formatted version of dump
cpdaa = "QC/REF_posID_level_18827cpdaa_4535UKBIDs.csv" # cpdaa positions
referenceseq = "QC/CCDSfromfasta_REFERENCE.csv" # sequences that cpdaa and missense positions are based on
qcout = "features_summary_QCd_sequences.csv"

# get directory and set, returns the directory component of os.path
FOLDERNAME = "Uniprot_Domains"
folder = path.dirname(path.abspath(__file__))
marktime = strftime("%Y%m%d") # subfolder name= today's date
output_folder = folder + "/" + FOLDERNAME + "/" + marktime + "/"
if not path.exists(output_folder):
    makedirs(output_folder) # updates directories, if directory name doesnt exist, creates dir

# INPUTS to queryUKB()
SAVEFULL = output_folder + FULLFILE
SAVESUMMARY = output_folder + SUMFILE

# INPUTS to Sequence QC
REFSEQ = folder + "/" + referenceseq
CpDAA = folder + "/" + cpdaa
QCoutput = output_folder + "/" + qcout

def parseSite(row, splitKey):
    # called by queryUKB, for every protein row returns pos of single sites ex. P12335 1; 3; 6; 7
    row = row.split(splitKey)
    coords = []
    for item in row:
        # removed two whitespaces in partition for case of no note/evidence- replaced split("..") w strip() 
        coords.append(item.partition(";")[0].strip()) 
        coords2 = []
        for item in coords:
            if len(item) > 1:
                # checks if there is isoform id preface or if range for other type of annotation
                if "-" not in item and " " not in item and ".." not in item:
                    coords2.append(item)
            else:
                continue
    return(coords2)

def parseRegion(row, splitKey):
    # called by queryUKB,for every protein row returns list of feature position ranges ex. P12335 1-3; 6-10
    coords = []
    row = row.split(splitKey)
    for item in row:
        if ".." in item:
            coords.append(item.partition(";  ")[0].split(".."))  # checked if its a range 1st
    coords2 = []
    for item in coords:
        if len(item) > 1:
            # checks if there is isoform id preface or if range for other type of annotation
            if "-" not in item[0] and " " not in item[0]:
                coords2.append(item)
    ranges = [x + "-" + y for x, y in coords2]
    return(ranges)

def queryUKB(SAVEFULL, SAVESUMMARY, FOLDERNAME):
    # online query to ukb api
    download_link = "https://www.uniprot.org/uniprot" \
            "?query=reviewed:yes" \
            "&fil=organism:%22Homo%20sapiens%20(Human)%20[9606]%22" \
            "&format=tab" \
            "&compress=no" \
            "&force=yes" \
            "&columns=id," \
            "entry%20name," \
            "protein%20names," \
            "genes(PREFERRED)," \
            "length," \
            "sequence," \
            "comment(SUBUNIT%20STRUCTURE)," \
            "3d," \
            "database(PDB)," \
            "interactor," \
            "comment(POST-TRANSLATIONAL%20MODIFICATION)," \
            "feature(BINDING%20SITE)," \
            "feature(ACTIVE%20SITE)," \
            "feature(METAL%20BINDING)," \
            "feature(MODIFIED%20RESIDUE)," \
            "feature(LIPIDATION)," \
            "feature(GLYCOSYLATION)," \
            "feature(MUTAGENESIS)," \
            "feature(NP%20BIND)," \
            "feature(DNA%20BINDING)," \
            "feature(ZINC%20FINGER)," \
            "feature(MOTIF)," \
            "feature(REGION)," \
            "feature(DOMAIN%20EXTENT)," \
            "feature(DISULFIDE%20BOND)," \
            "feature(BETA%20STRAND)," \
            "feature(HELIX)," \
            "feature(TURN)" \
           "&sort=score"
    
    # saving full ukb query results
    download_database = get(download_link).text
    with open(SAVEFULL, "w") as f:
        f.write(download_database)

    # read full ukb file and create feature summary file
    with open(SAVEFULL, "r") as qdata, \
        open(SAVESUMMARY, "w") as output:
        qdata_csv = reader(qdata, delimiter = "\t")
        next(qdata) # skip header
        # write header to summary file
        output.write("UKBID" + "\t" + \
            "Protein" + "\t" + \
            "Gene" + "\t" + \
            "Length" + "\t" + \
            "Sequence" + "\t" + \
            "SubunitStructure" + "\t" + \
            "3D" + "\t" + \
            "PDB" + "\t" + \
            "InteractsWith" + "\t" + \
            "PTM" + "\t" + \
            "SubstrateBindingSite" + "\t" + \
            "ActiveSite" + "\t" + \
            "MetalBindSite" + "\t" + \
            "ModifiedResidue" + "\t" + \
            "Lipidation" + "\t" + \
            "Glycosylation" + "\t" + \
            "Mutagenesis" + "\t" + \
            "NucleotideBindRegion" + "\t" + \
            "DNAbindRegion" + "\t" + \
            "ZincFinger" + "\t" + \
            "Motif" + "\t" + \
            "DisorderedRegion" + "\t" + \
            "Domain" + "\t" + \
            "DisulfideBond" + "\t" + \
            "BetaStrand" + "\t" + \
            "Helix" + "\t" + \
            "Turn" + "\n")
        
        # save row values as variables to functions to use
        for line in qdata_csv:
            entry = line[0]             #output.write
            protein_names = line[2]      #output.write
            gene_names = line[3]        #output.write
            length = line[4]            #output.write
            sequence = line[5]          #output.write
            subunit = line[6]         #output.write
            d3 = line[7]              #output.write
            pdb = line[8]             #output.write
            interactors = line[9]     #output.write
            ptmcc = line[10]          #output.write
            # all below columns required parsing
            binding = line[11]           #output.write.parsed
            active = line[12]           #output.write.parsed
            metal = line[13]               #output.write.parsed
            modresidue = line[14]                #output.write.parsed
            lipid = line[15]   #output.write.parsed
            glycos = line[16]           #output.write.parsed
            mutagen = line[17]              #output.write.parsed
            nucleotide = line[18]             #output.write.parsed
            dna = line[19]            #output.write.parsed
            zinc = line[20]           #output.write.parsed
            motif = line[21]           #output.write.parsed
            region = line[22]   #output.write.parsed
            domain = line[23]   #output.write.parsed
            bonds = line[24]   #output.write.parsed
            beta = line[25]             #output.write.parsed
            helix = line[26]            #output.write.parsed
            turn = line[27]             #output.write.parsed
        
            # PARSING SITE ANNOTATIONS-----------------------
            active_coords2 = parseSite(active, "ACT_SITE ")
            metal_coords2 = parseSite(metal, "METAL ")
            mod_coords2 = parseSite(modresidue, "MOD_RES ")
            lipid_coords2 = parseSite(lipid, "LIPID ")
            gly_coords2 = parseSite(glycos, "CARBOHYD ")
            mut_coords2 = parseSite(mutagen, "MUTAGEN ")
            
            # PARSING MULTI-POSITION ANNOTATIONS---------------
            npbind_ranges = parseRegion(nucleotide, "NP_BIND ")
            dnabind_ranges = parseRegion(dna, "DNA_BIND ")
            zinc_ranges = parseRegion(zinc, "ZN_FING ")
            motif_ranges = parseRegion(motif, "MOTIF ")
            domain_ranges = parseRegion(domain, "DOMAIN ")
            strand_ranges = parseRegion(beta, "STRAND ")
            helix_ranges = parseRegion(helix, "HELIX ")
            turn_ranges = parseRegion(turn, "TURN ")
            
            #-----special case 1: binding site-substrate-------
            binding = binding.split("BINDING ")
            """
            **what parsing code does..**
            * splits string on certain all caps word like BINDING, ZN_FING, STRAND, DOMAIN, etc.
            * ignores positions prefaced by something other than all caps word used to split string and positions prefaced by 
                isoform id, ex P1234-3:2..4 ignored and TURN 10..14 ignored for BINDING split()
            * only parses positions that are describing 'Substrate' binding site
            * parses single position only and ignores 2 position annotation (checks for ..)
            * outputs single positions as list
            """
            bind_coords=[]
            for item in binding:
                if "Substrate" in item:
                    # removed two white spaces in partition in case there is no note or evidence
                    # replaced split("..") with strip() to remove whitespaces
                    bind_coords.append(item.partition(";")[0].strip())
            bind_coords2 = []
            for item in bind_coords:
                if len(item) > 1:
                    # checks if there is isoform id preface or if range for other type of annotation
                    if "-" not in item and " " not in item and ".." not in item:
                        bind_coords2.append(item)
                else:
                    continue
               
            #-----special case 2: disulfide bond----------
            bond_coords2 = []
            bonds = bonds.replace(" DISULFID ", "")
            bonds = bonds.replace("DISULFID ", "")
            bonds=bonds.split(";")
            for item in bonds:
                if len(item) > 1:
                    # 2 positions
                    if '..' in item:
                        bond_coords = [item.split("..")[0], item.split("..")[1]]
                        # checks if there is isoform id preface or if range for other type of annotation
                        if "-" not in bond_coords[0] and " " not in bond_coords[0]:
                            bond_coords2.append(bond_coords)
                    # 1 position
                    if '..' not in item:
                        # checks if there is isoform id preface or if range for other type of annotation
                        if "-" not in item and " " not in item and "=" not in item:
                            v = []
                            v.append(item)
                            v.append(item)
                            bond_coords2.append(v)
            bridges = [x + ".." + y for x, y in bond_coords2]
            
            #-----special case 3: region-disordered--------------
            region = region.split("REGION ")
            disorder_coords=[]
            for item in region:
                if "Disordered" in item:
                    disorder_coords.append(item.partition(";  ")[0].split(".."))
            disorder_coords2 = []
            for item in disorder_coords:
                if len(item) > 1:
                    # checks if there is isoform id preface or if range for other type of annotation
                    if "-" not in item[0] and " " not in item[0]:
                        disorder_coords2.append(item)
                elif len(item) == 1:
                    # case of single position, make into range format
                    # checks if there is isoform id preface or if range for other type of annotation
                    if "-" not in item[0] and " " not in item[0]:
                        v = []
                        v.append(item[0])
                        v.append(item[0])
                        disorder_coords2.append(v)
                else:
                    continue
            disorder_ranges = [x + "-" + y for x, y in disorder_coords2]
           
            # write final parsed line to output
            output.write(entry + "\t" + protein_names + "\t" + \
                gene_names + "\t" + length + "\t" + sequence + "\t" + \
                subunit + "\t" + d3 + "\t" + pdb + "\t" + \
                interactors + "\t" + ptmcc + "\t" + \
                "; ".join(bind_coords2) + "\t" + \
                "; ".join(active_coords2) + "\t" + \
                "; ".join(metal_coords2) + "\t" + \
                "; ".join(mod_coords2) + "\t" + \
                "; ".join(lipid_coords2) + "\t" + \
                "; ".join(gly_coords2) + "\t" + \
                "; ".join(mut_coords2) + "\t" + \
                "; ".join(npbind_ranges) + "\t" + \
                "; ".join(dnabind_ranges) + "\t" + \
                "; ".join(zinc_ranges) + "\t" + \
                "; ".join(motif_ranges) + "\t" + \
                "; ".join(disorder_ranges) + "\t" + \
                "; ".join(domain_ranges) + "\t" + \
                "; ".join(bridges) + "\t" + \
                "; ".join(strand_ranges) + "\t" + \
                "; ".join(helix_ranges) + "\t" + \
                "; ".join(turn_ranges) + "\n")


def compareUKBsequences(summaryFile, refFile, qcname):
    dfsummary = read_pd_txt(summaryFile) # Sequence, UKBID
    ref = pd.read_csv(refFile) # proSequence, UKBID
    ids = list(set(ref.UKBID))
    dfsummary = addcolumnconditionalDropFalse(ids, dfsummary, 'UKBID') # filter summary df
    pepdict = dict(zip(ref.UKBID, ref.proSequence)) # make dict from ref IDs and sequences
    # new columns
    otherLen = []
    identical = []
    BminusA = []
    diffcol = []
    posdiffcol = []
    diffposID = []
    # loop over summary df
    for index, row in dfsummary.iterrows():
        pepSumr = row['Sequence'] # save summarydf pep var
        ukbID = row['UKBID'] # save ukb id
        pepRef = pepdict[ukbID] # save pep in dictionary matching ukb id
        pepSumr = str(pepSumr)
        pepRef = str(pepRef)
        lenSumr = len(pepSumr)
        lenRef = len(pepRef)
        # 100% identical case
        if lenRef == lenSumr:
            if pepRef == pepSumr:
                otherLen.append(lenSumr)
                identical.append(True)
                BminusA.append(0)
                diffcol.append(np.nan)
                posdiffcol.append(np.nan)
                diffposID.append(np.nan)
            # same length but not identical case
            else:
                diffaa_list = [pepRef[i] for i in range(len(pepSumr)) if pepSumr[i] != pepRef[i]]
                diffpos_list = [i+1 for i in range(len(pepSumr)) if pepSumr[i] != pepRef[i]]
                otherLen.append(lenSumr)
                identical.append(False)
                BminusA.append(0)
                diffcol.append(diffaa_list)
                posdiffcol.append(diffpos_list)
                posids = []
                for ii,v in enumerate(diffaa_list):
                    p = ukbID + '_' + v + str(diffpos_list[ii])
                    posids.append(p)
                posids = ';'.join(posids)
                diffposID.append(posids)
        # REF longer than Summaryfrom querySeq case
        if lenRef > lenSumr:
            diffaa_list = [pepRef[i] for i in range(len(pepSumr)) if pepSumr[i] != pepRef[i]]
            diffpos_list = [i+1 for i in range(len(pepSumr)) if pepSumr[i] !=
                               pepRef[i]]
            otherLen.append(lenRef)
            identical.append(False)
            BminusA.append((lenSumr-lenRef))
            diffcol.append(diffaa_list)
            posdiffcol.append(diffpos_list)
            posids = []
            for ii,v in enumerate(diffaa_list):
                p = ukbID + '_' + v + str(diffpos_list[ii])
                posids.append(p)
            posids = ';'.join(posids)
            diffposID.append(posids)
        # REF shorter than Summaryfrom querySeq case
        if lenRef < lenSumr:
            diffaa_list = [pepSumr[i] for i in range(len(pepRef)) if pepRef[i] != pepSumr[i]]
            diffpos_list = [i+1 for i in range(len(pepRef)) if pepRef[i] != pepSumr[i]]
            otherLen.append(lenRef)
            identical.append(False)
            BminusA.append((lenSumr-lenRef))
            diffcol.append(diffaa_list)
            posdiffcol.append(diffpos_list)
            posids = []
            for ii,v in enumerate(diffaa_list):
                p = ukbID + '_' + v + str(diffpos_list[ii])
                posids.append(p)
            posids = ';'.join(posids)
            diffposID.append(posids)
    #dfsummary.loc[:, 'REF.length'] = otherLen
    dfsummary.loc[:, 'identicalToRefSequence'] = identical # boolean
    dfsummary.loc[:, 'query.minus.REF.difference'] = BminusA
    # aa diff based on shorter seq
    #dfsummary.loc[:, 'residue.difference'] = diffcol
    # position diff based on shorter seq
    #dfsummary.loc[:, 'position.difference'] = posdiffcol
    dfsummary.loc[:, 'QC.posID.difference'] = diffposID
    print()
    checkColumnValues(dfsummary, "identicalToRefSequence")
    dfsummary.to_csv(qcname, index=False)

if __name__ == '__main__':
    queryUKB(SAVEFULL, SAVESUMMARY, FOLDERNAME)
    compareUKBsequences(SAVESUMMARY, REFSEQ, QCoutput)

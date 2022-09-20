"""
feb.4 2022 $python downloadUniprotTopological.py 
code derived from SNAP-IT Update_Annotations.py

created folder "Uniprot_Topological"
created sub-folder "$ todays date"
output files in sub-folder:
    1. full_human_topo_domain.txt
    2. cysteine_summary.txt

updates:
2/6/22 fixed natural variants in extra cellular domain range counts by adding try and except statement

"""
from requests import get
from time import strftime
from os import path, makedirs
from csv import reader

# Set directories, returns the directory component of a pathname os.path
folder = path.dirname(path.abspath(__file__))

# Update with today's date
# time.strftime: format time, %m is month, %Y is year %d is day %M is minute
today = strftime("%Y%m%d")

# Create updated directories, if directory name doesnt exist, make it
updated_annotation_folder = folder + "/Uniprot_Topological/" + today + "/"
if not path.exists(updated_annotation_folder):
    makedirs(updated_annotation_folder)

"""
This code currently filters for human proteins that have their topological domains annotated. To see full list of downloadable annotation columns, browse the 'Column' button in the UniProtKb query <https://www.uniprot.org/help/uniprotkb_column_names>.

To search this query on uniprot.org, remove the commands for 'format', 'compress', and 'force'.

url copied from share button in browser, the output of this url is 3901 row table displayed on website, and original script output file full_human_topo_domain.txt is 3904 rows.

https://www.uniprot.org/uniprot/?query=reviewed%3Ayes%20organism%3A%22Homo%20sapiens%20(Human)%20%5B9606%5D%22%20annotation%3A(type%3Atopo_dom)&columns=id%2Centry%20name%2Cprotein%20names%2Cgenes(PREFERRED)%2Clength%2Cfeature(NATURAL%20VARIANT)%2Csequence%2Ckeywords%2Cgo%2Cfeature(TOPOLOGICAL%20DOMAIN)%2Ccomment(SUBCELLULAR%20LOCATION)%2Cfeature(TRANSMEMBRANE)%2Cfeature(INTRAMEMBRANE)%2Cfeature(DISULFIDE%20BOND)&sort=score

-------------------------------------------------------------------------------
TESTING if script works with url link from share button in browser, it doesnt

CONVERTING URL to string that works with script:
    %3A is :
    %5B is [
    %5D is ]
    %2C is ,

FAILED w/
#  "&fil=reviewed:yes%20organism:%22Homo%20sapiens%20(Human)%20[9606]%22" \
"""
download_link = "https://www.uniprot.org/uniprot" \
       "?query=%22topological%20domain%22" \
       "&fil=organism:%22Homo%20sapiens%20(Human)%20[9606]%22" \
       "&format=tab" \
       "&compress=no" \
       "&force=yes" \
       "&columns=id," \
       "entry%20name," \
       "protein%20names," \
       "genes(PREFERRED)," \
       "length," \
       "feature(NATURAL%20VARIANT)," \
       "sequence," \
       "keywords," \
       "go," \
       "feature(TOPOLOGICAL%20DOMAIN)," \
       "comment(SUBCELLULAR%20LOCATION)," \
       "feature(TRANSMEMBRANE)," \
       "feature(INTRAMEMBRANE)," \
       "feature(DISULFIDE%20BOND)" \
       "&sort=score"

file = updated_annotation_folder + "full_human_topo_domain.txt"
output_file = updated_annotation_folder + "cysteine_summary.txt"

download_database = get(download_link).text
with open(file, "w") as f:
    f.write(download_database)

with open(file, "r") as topo, \
    open(output_file, "w") as output:
    topo_csv = reader(topo, delimiter = "\t")
    next(topo)

    output.write("UniProtKb" + "\t" + "Gene" + "\t" + "Protein" + "\t" + "Extracellular_Domains" + "\t" + "Natural_C_Variants" + "\t" + "Acquired_C" + "\t" + "Acquired_C_from_R" + "\t" + "Lost_C" + "\t" + "Extracellular_Sequences" + "\t" + "Extracellular_C" + "\t" + "Extracellular_Disulfide_Bonds" + "\n")

    for line in topo_csv:
        entry = line[0]
        entry_name = line[1]
        protein_names = line[2]
        gene_names = line[3]
        length = line[4]
        natural_variant = line[5]
        sequence = line[6]
        keywords = line[7]
        gene_ontology = line[8]
        topological_domain = line[9]
        subcellular_location = line[10]
        transmembrane = line[11]
        intramembrane = line[12]
        disulfide_bond = line[13]

        # Get extracullar domain range from topological domain column
        extracell_coords = []
        domains = topological_domain.split("TOPO_DOM ")
        for item in domains:
            if "Extracellular" in item:
                extracell_coords.append(item.partition(";  ")[0].split(".."))
        extracell_coords2 = []
        for item in extracell_coords:
            if len(item) > 1:
                extracell_coords2.append(item)
            elif len(item) == 1:
                v = []
                v.append(item[0])
                v.append(item[0])
                extracell_coords2.append(v)
            else:
                continue
        extracell_ranges = [x + "-" + y for x, y in extracell_coords2]

        # Check for acquired cysteine within extracellular range from natural variant column
        ##(ignored) VARIANT 634..635;
        ##(ignored) VARIANT Q9NZC2-2:200;
        ##(parsed) VARIANT 634;  /note="C -> CHELC
        variant_coords = []
        acquired_C = 0
        acquired_C_from_R = 0
        lost_C = 0
        
        variants = natural_variant.split("VARIANT ")
        for item in variants:
            if "-> C" in item and ".." not in item:
                print("gain of C item: ", item)
                for v in extracell_coords2:
                    a = v[0]
                    b = v[1]
                    try:
                        # isoform ukb ids filtered out with try statement
                        varpos = int(item.split(";")[0])
                        if int(a) <= varpos <= int(b):
                            variant_coords.append([item.split(";")[0], item.split("/")[1].replace("; ","").strip().replace('\"',""), item.split("/")[2].replace("; ","").strip().replace('\"',"")])
                            acquired_C += 1
                            if "R ->" in item:
                                acquired_C_from_R += 1
                    except ValueError:
                        continue
            if "C ->" in item and ".." not in item:
                print("loss of C item: ", item)
                for v in extracell_coords2:
                    a = v[0]
                    b = v[1]
                    try:
                        # isoform ukb ids filtered out with try statement
                        varpos = int(item.split(";")[0])
                        if int(a) <= varpos <= int(b):
                            variant_coords.append([item.split(";")[0], item.split("/")[1].replace("; ","").strip().replace('\"',""), item.split("/")[2].replace("; ","").strip().replace('\"',"")])
                            lost_C += 1
                    except ValueError:
                        continue
        variant_events = [x + y.replace("note="," ").replace(" -> ","/") + z.replace("evidence="," ") for x, y, z in variant_coords] # to simplify row values, only include x and y values (z is notes)

        # Substring extracellular peptide sequence and count number of extracellular cysteines
        # correct error with position index, must subtract 1 from pos to get actual aa at the aapos
        extracell_seqs = []
        count_C = 0
        for item in extracell_coords2:
            try:
                extracell_seqs.append(sequence[int(item[0])-1 : int(item[1])])
                count_C += sequence[int(item[0])-1 : int(item[1])].count("C")
            except ValueError:
                pass

        # Check for disulfide bonds within extracellular range
        # only returns bridges that are within range
        bond_coords = []
        bonds = disulfide_bond.replace(" DISULFID ", "")
        bonds = bonds.replace("DISULFID ", "")
        bonds = bonds.split(";")
        for item in bonds:
            if '..' in item:
                coords = [item.split("..")[0], item.split("..")[1]]
                for v in extracell_coords2:
                    a = v[0]
                    b = v[1]
                    try:
                        if int(a) <= int(coords[0]) <= int(b) or int(a) <= int(coords[1]) <= int(b):
                            bond_coords.append(coords)
                    except ValueError:
                        pass
        bridges = [x + ".." + y for x, y in bond_coords]

        output.write(entry + "\t" + entry_name + "\t" + gene_names + "\t" + ", ".join(extracell_ranges) + "\t" + ", ".join(variant_events) + "\t" + str(acquired_C) + "\t" + str(acquired_C_from_R) + "\t" + str(lost_C) + "\t" + ", ".join(extracell_seqs)  + "\t" + str(count_C) + "\t" + ", ".join(bridges) + "\n")

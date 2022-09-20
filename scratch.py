















# redundancy  for function testing

        active = active.split("ACT_SITE ")
        active_coords=[]
        for item in active:
            # removed two white spaces in partition in case there is no note or evidence
            # replaced split("..") with strip() to remove whitespaces
            active_coords.append(item.partition(";")[0].strip()) 
        active_coords2 = []
        for item in active_coords:
            if len(item) > 1:
                # checks if there is isoform id preface or if range for other type of annotation
                if "-" not in item and " " not in item and ".." not in item:
                    active_coords2.append(item)
            else:
                continue

        metal = metal.split("METAL ")
        metal_coords=[]
        for item in metal:
            # removed two white spaces in partition in case there is no note or evidence
            # replaced split("..") with strip() to remove whitespaces
            metal_coords.append(item.partition(";")[0].strip()) 
        metal_coords2 = []
        for item in metal_coords:
            if len(item) > 1:
                # checks if there is isoform id preface or if range for other type of annotation
                if "-" not in item and " " not in item and ".." not in item:
                    metal_coords2.append(item)
            else:
                continue


        modresidue = modresidue.split("MOD_RES ")
        mod_coords=[]
        for item in modresidue:
            # removed two white spaces in partition in case there is no note or evidence
            # replaced split("..") with strip() to remove whitespaces
            mod_coords.append(item.partition(";")[0].strip()) 
        mod_coords2 = []
        for item in mod_coords:
            if len(item) > 1:
                # checks if there is isoform id preface or if range for other type of annotation
                if "-" not in item and " " not in item and ".." not in item:
                    mod_coords2.append(item)
            else:
                continue

        lipid = lipid.split("LIPID ")
        lipid_coords=[]
        for item in lipid:
            # removed two white spaces in partition in case there is no note or evidence
            # replaced split("..") with strip() to remove whitespaces
            lipid_coords.append(item.partition(";")[0].strip()) 
        lipid_coords2 = []
        for item in lipid_coords:
            if len(item) > 1:
                # checks if there is isoform id preface or if range for other type of annotation
                if "-" not in item and " " not in item and ".." not in item:
                    lipid_coords2.append(item)
            else:
                continue

        glycos = glycos.split("CARBOHYD ")
        gly_coords=[]
        for item in glycos:
            # removed two white spaces in partition in case there is no note or evidence
            # replaced split("..") with strip() to remove whitespaces
            gly_coords.append(item.partition(";")[0].strip()) 
        gly_coords2 = []
        for item in gly_coords:
            if len(item) > 1:
                # checks if there is isoform id preface or if range for other type of annotation
                if "-" not in item and " " not in item and ".." not in item:
                    gly_coords2.append(item)
            else:
                continue

        mutagen = mutagen.split("MUTAGEN ")
        mut_coords=[]
        for item in mutagen:
            # removed two white spaces in partition in case there is no note or evidence
            # replaced split("..") with strip() to remove whitespaces
            mut_coords.append(item.partition(";")[0].strip()) 
        mut_coords2 = []
        for item in mut_coords:
            if len(item) > 1:
                # checks if there is isoform id preface or if range for other type of annotation
                if "-" not in item and " " not in item and ".." not in item:
                    mut_coords2.append(item)
            else:
                continue
#########################################################################################################################################################


        nucleotide = nucleotide.split("NP_BIND ")
        nuc_coords=[]
        for item in nucleotide:
            # check if its a range
            if ".." in item:
                # ignores non NP_BIND #..# in item that begins with NP_BIND by spliting on NP_BIND
                nuc_coords.append(item.partition(";  ")[0].split(".."))
        nuc_coords2 = []
        for item in nuc_coords:
            if len(item) > 1:
                # checks if there is isoform id preface or if range for other type of annotation
                if "-" not in item[0] and " " not in item[0]:
                    nuc_coords2.append(item)
        npbind_ranges = [x + "-" + y for x, y in nuc_coords2]

        dna = dna.split("DNA_BIND ")
        dna_coords=[]
        for item in dna:
            # check if its a range
            if ".." in item:
                # ignores non DNA_BIND prefaced ranges in item by spliting on DNA_BIND
                dna_coords.append(item.partition(";  ")[0].split(".."))

        dna_coords2 = []
        for item in dna_coords:
            if len(item) > 1:
               # checks if there is isoform id preface or if range for other type of annotation
                if "-" not in item[0] and " " not in item[0]:
                    dna_coords2.append(item)
        dnabind_ranges = [x + "-" + y for x, y in dna_coords2]


        zinc = zinc.split("ZN_FING ")
        z_coords=[]
        for item in zinc:
            # check if its a range
            if ".." in item:
                # ignores non ZN_FING prefaced ranges in item by spliting on ZN_FING
                z_coords.append(item.partition(";  ")[0].split(".."))
        z_coords2 = []
        for item in z_coords:
            if len(item) > 1:
                # checks if there is isoform id preface or if range for other type of annotation
                if "-" not in item[0] and " " not in item[0]:
                    z_coords2.append(item)
        zinc_ranges = [x + "-" + y for x, y in z_coords2]


        motif = motif.split("MOTIF ")
        motif_coords=[]
        for item in motif:
            # check if its a range
            if ".." in item:
                # ignores non MOTIF prefaced ranges in item by spliting on MOTIF
                motif_coords.append(item.partition(";  ")[0].split(".."))
        motif_coords2 = []
        for item in motif_coords:
            if len(item) > 1:
                # checks if there is isoform id preface or if range for other type of annotation
                if "-" not in item[0] and " " not in item[0]:
                    # ignore ranges for non-canonical sequence since im only checking pos and seq identify of canonical 
                    motif_coords2.append(item)
        motif_ranges = [x + "-" + y for x, y in motif_coords2]


        domain = domain.split("DOMAIN ")
        domain_coords=[]
        for item in domain:
            # check if its a range
            if ".." in item:
                # ignores non DOMAIN prefaced ranges in item by spliting on DOMAIN
                domain_coords.append(item.partition(";  ")[0].split(".."))
        domain_coords2 = []
        for item in domain_coords:
            if len(item) > 1:
                # checks if there is isoform id preface or if range for other type of annotation
                if "-" not in item[0] and " " not in item[0]:
                    domain_coords2.append(item)
        domain_ranges = [x + "-" + y for x, y in domain_coords2]

        strand_coords = []
        beta = beta.split("STRAND ")
        for item in beta:
            # check if its a range
            if ".." in item:
                # ignores non STRAND #..# in item by spliting on STRAND
                strand_coords.append(item.partition(";  ")[0].split(".."))
        strand_coords2 = []
        for item in strand_coords:
            if len(item) > 1:
                # checks if there is isoform id preface or if range for other type of annotation
                if "-" not in item[0] and " " not in item[0]:
                    strand_coords2.append(item)
        strand_ranges = [x + "-" + y for x, y in strand_coords2]

        helix_coords = []
        helix = helix.split("HELIX ")
        for item in helix:
            # check if its a range
            if ".." in item:
                # ignores non HELIX #..# in item by spliting on HELIX
                helix_coords.append(item.partition(";  ")[0].split(".."))
        helix_coords2 = []
        for item in helix_coords:
            if len(item) > 1:
                # checks if there is isoform id preface or if range for other type of annotation
                if "-" not in item[0] and " " not in item[0]:
                    helix_coords2.append(item)
        helix_ranges = [x + "-" + y for x, y in helix_coords2]

        turn_coords = []
        turn = turn.split("TURN ")
        for item in turn:
            # check if its a range
            if ".." in item:
                # ignores non TURN #..# in item by spliting on TURN
                turn_coords.append(item.partition(";  ")[0].split(".."))
        turn_coords2 = []
        for item in turn_coords:
            if len(item) > 1:
                # checks if there is isoform id preface or if range for other type of annotation
                if "-" not in item[0] and " " not in item[0]:
                    turn_coords2.append(item)
        turn_ranges = [x + "-" + y for x, y in turn_coords2]

import pandas as pd
import gumpy, numpy

count = 0
snp_count = []
def parse_who_catalog(filename):
    '''Parses the WHO TB catalog

    Args:
        filename (str): Path to the WHO catalog
    Returns:
        pd.dataframe: Dataframe containing the mutations
    '''
    df = pd.read_excel(filename, sheet_name="Genome_indices")
    return df

def rev_comp_snp(reference, gene, pos, ref, alt, masks):
    '''Convert a mutation into the appropriate number of SNPs in GARC, converting to amino acids as required

    Args:
        reference (gumpy.Genome): Reference genome object
        gene (str): Gene name
        pos (int): Genome index
        ref (str): Reference base(s)
        alt (str): Mutant base(s)
        masks (dict): Dictionary used for caching the masks required for rebuilding genes
    Returns:
        list(str): List of SNP mutations in GARC
    '''
    mutations = []
    ref_seq = reference.nucleotide_sequence.copy()
    for (index, (r, a)) in enumerate(zip(ref, alt)):
        if r != a:
            if (pos + index) - reference.genes_lookup[gene]["start"] <= 0:
                #Past the end of the gene so just return
                return mutations
            if reference.genes_lookup[gene]["end"] - (pos + index) < 0 or reference.genes[gene].codes_protein == False:
                #Promoter or non-coding so return the difference in nucleotides
                r,a  = gumpy.Gene._complement([r, a])
                mutations.append(gene + "@" + r + str(reference.genes_lookup[gene]["end"] - (pos + index)) + a)
            else:
                ref_seq[pos + index - 1] = a
    if reference.genes[gene].codes_protein:
        stacked_mask, mask = masks[gene]

        g = gumpy.Gene(   name=gene,\
                    nucleotide_sequence=ref_seq[mask],\
                    index=reference.nucleotide_index[mask],\
                    nucleotide_number=reference.stacked_nucleotide_number[stacked_mask],\
                    is_cds=reference.stacked_is_cds[stacked_mask],\
                    is_promoter=reference.stacked_is_promoter[stacked_mask],\
                    is_indel=reference.is_indel[mask],
                    indel_length=reference.indel_length[mask],
                    codes_protein=reference.genes_lookup[gene]['codes_protein'],\
                    reverse_complement=reference.genes_lookup[gene]['reverse_complement'],\
                    feature_type=reference.genes_lookup[gene]['type'])
        
        aa_mut = [(i+1, ref_aa, alt_aa) for (i, (ref_aa, alt_aa)) in enumerate(zip(reference.genes[gene].codons, g.codons)) if ref_aa != alt_aa]
        for (pos_, r, a) in aa_mut:
            r = g.codon_to_amino_acid[r]
            a = g.codon_to_amino_acid[a]
            if r == a:
                mutations.append(gene + "@" + str(pos_) + "=")
            else:
                mutations.append(gene + "@" + r + str(pos_) + a)
    return mutations

def rev_comp(reference, gene, pos, ref, alt, masks):
    '''Handle reverse complement changes required

    Args:
        reference (gumpy.Genome): Reference genome object
        gene (str): Name of the gene
        pos (int): Genome index
        ref (str): Reference bases
        alt (str): Mutant bases
        masks (dict): Dictionary used for caching the masks for rebuilding genes
    Returns:
        list(str): List of mutations in GARC
    '''
    global count
    global snp_count
    if len(ref) == len(alt):
        #SNP so convert to AA too
        return rev_comp_snp(reference, gene, pos, ref, alt, masks)
    elif len(ref) > len(alt):
        #Del
        snp = rev_comp_snp(reference, gene, pos, ref, alt, masks)
        count += 1
        snp_count.append(snp)
        dels = [gene + "@" + str((reference.genes_lookup[gene]["end"] - (pos + (len(ref) - len(alt))))) + "_del_" + "".join(gumpy.Gene._complement(list(ref[len(alt):])))]
        return dels + snp
    elif len(ref) < len(alt):
        #Ins
        snp = rev_comp_snp(reference, gene, pos, ref, alt, masks)
        count += 1
        snp_count.append(snp)
        ins = [gene + "@" + str((reference.genes_lookup[gene]["end"] - (pos + (len(alt) - len(ref))))) + "_ins_" + "".join(gumpy.Gene._complement(list(alt[len(ref):])))]
        return ins + snp

def snps(reference, gene, pos, ref, alt, masks):
    '''Convert a mutation into the appropriate number of SNPs in GARC, converting to amino acids as required

    Args:
        reference (gumpy.Genome): Reference genome object
        gene (str): Gene name
        pos (int): Genome index
        ref (str): Reference base(s)
        alt (str): Mutant base(s)
        masks (dict): Dictionary used for caching the masks required for rebuilding genes
    Returns:
        list(str): List of SNP mutations in GARC
    '''
    mutations = []
    ref_seq = reference.nucleotide_sequence.copy()
    for (index, (r, a)) in enumerate(zip(ref, alt)):
        if r != a:
            if reference.genes_lookup[gene]["end"] - (pos + index) <= 0:
                #Past the end of the gene so just return
                return mutations
            if (pos + index) - reference.genes_lookup[gene]["start"] < 0 or reference.genes[gene].codes_protein == False:
                #Promoter or non-coding so return the difference in nucleotides
                mutations.append(gene + "@" + r + str((pos + index) - reference.genes_lookup[gene]["start"]) + a)
            else:
                ref_seq[pos + index - 1] = a
    if reference.genes[gene].codes_protein:
        stacked_mask, mask = masks[gene]

        g = gumpy.Gene(   name=gene,\
                    nucleotide_sequence=ref_seq[mask],\
                    index=reference.nucleotide_index[mask],\
                    nucleotide_number=reference.stacked_nucleotide_number[stacked_mask],\
                    is_cds=reference.stacked_is_cds[stacked_mask],\
                    is_promoter=reference.stacked_is_promoter[stacked_mask],\
                    is_indel=reference.is_indel[mask],
                    indel_length=reference.indel_length[mask],
                    codes_protein=reference.genes_lookup[gene]['codes_protein'],\
                    reverse_complement=reference.genes_lookup[gene]['reverse_complement'],\
                    feature_type=reference.genes_lookup[gene]['type'])
        
        aa_mut = [(i+1, ref_aa, alt_aa) for (i, (ref_aa, alt_aa)) in enumerate(zip(reference.genes[gene].codons, g.codons)) if ref_aa != alt_aa]
        for (pos_, r, a) in aa_mut:
            r = g.codon_to_amino_acid[r]
            a = g.codon_to_amino_acid[a]
            if r == a:
                mutations.append(gene + "@" + str(pos_) + "=")
            else:
                mutations.append(gene + "@" + r + str(pos_) + a)
    return mutations

class filterList(list):
    def __eq__(self, other):
        if len(self) != len(other):
            # print("Different lengths")
            return False
        check = True
        for (x, y) in zip(self, other):
            if x is not None and y is not None:
                #Ignore None values
                check = check and x == y
        return check


def del_calls(ref, alt, added=0):
    '''Deal with del calls. Attempts to identify dels mid-sequence. If a repeated base is deleted (aaa->aa), it is assumed
        that the first base is deleted.
    
    Args:
        ref (list): Reference bases
        alt (list): Alternative bases
        added (int, optional): Number of blanks added to check for mid-seq dels
    Returns:
        list(str): List of mutations in GARC
    '''
    print(ref, alt, len(ref), len(alt), added)
    input("?")
    if len(ref) == len(alt):
        return ref, alt
    for index in range(len(ref)):
        if index < len(alt) and ref[index] != alt[index]:
            if len(alt) >= len(ref):
                #Definitely not a mid-seq del
                return [], []
            #The bases differ, so check if this was a mid-seq del
            alt.insert(index, None)
    if len(ref) != len(alt) and added < 3:
        ref, alt = del_calls(ref, alt, added=added+1)
    return ref, alt



def to_garc(reference, gene, pos, ref, alt, masks):
    '''Convert to GARC

    Args:
        reference (gumpy.Genome): Reference genome object
        pos (int): Genome index
        ref (str): Reference base(s)
        alt (str): Mutant base(s)
    Returns:
        list(str): List of mutations in GARC
    '''
    global count
    global snp_count
    #Reverse complement genes need some alterations
    if reference.genes[gene].reverse_complement:
        return rev_comp(reference, gene, pos, ref, alt, masks)
    #Just convert to GARC
    else:
        if len(ref) == len(alt):
            return snps(reference, gene, pos, ref, alt, masks)
        elif len(ref) > len(alt):
            '''
            Indels are weirder than I first thought, they don't always indicate a del/ins at the end of the seq
                e.g. agctctagtg -> agtctagta has a `c` being deleted mid-seq (seq[2])
            With this kind of mid-seq indel, it is not possible to determine which position an indel occurs at, especially as there are often SNPs too:
                e.g:
                    tccggtctg -> a is ambiguous as to which values are delted/SNPs so the positions reported could be wrong
                    accg -> a is ambiguous as to if this is del(ccg) or del some other 3 bases and a SNP
            Also, GARC does not have syntax to suport both insertions and deltions simaltaneously

            '''
            #Del
            snp = snps(reference, gene, pos, ref, alt, masks)
            count += 1
            snp_count.append(snp)
            dels = del_calls(filterList(ref), filterList(alt))
            dels = [gene + "@" + str(((pos + (len(ref) - len(alt))) - reference.genes_lookup[gene]["start"])) + "_del_" + ref[len(alt):]]
            return snp + dels
        elif len(ref) < len(alt):
            snp = snps(reference, gene, pos, ref, alt, masks)  
            count += 1
            snp_count.append(snp)     
            ins = [gene + "@" + str(((pos + (len(alt) - len(ref))) - reference.genes_lookup[gene]["start"])) + "_ins_" + alt[len(ref):]]     
            return ins + snp


def get_masks(reference, gene):
    '''Find the numpy masks for the arrays within the reference genome for the specified gene.
    The masks are used to speed up the instanciation of new Gene objects for SNP finding. Finding the masks
        takes some time, so this is cached so the mask only has to be found once per gene rather than once per row
    
    Args:
        reference (gumpy.Genome): Reference genome object
        gene (str): Gene name
    Returns:
        (numpy.array, numpy.array): Tuple of 2 numpy arrays. First denotes the stacked mask for mutli-dimensional
                                    attributes. Second denotes the mask for 1D attributes
    '''
    #The mask for all stacked arrays (N-dim)
    stacked_mask = reference.stacked_gene_name == gene
    #The mask for singular arrays (1-dim) by collapsing stacked mask to 1-dim
    mask = numpy.any(stacked_mask, axis=0)
    return stacked_mask, mask

if __name__ == "__main__":
    #Load the reference genome
    reference = gumpy.Genome.load("reference.json.gz")

    #Load the catalogue
    data = parse_who_catalog("WHO-UCN-GTB-PCI-2021.7-eng.xlsx")

    #The catalogue is grouped by gene, so we can store gene masks until they are no longer required
    masks = {}

    #Setup details for drugs
    drug_columns = [
        'RIF_Conf_Grade',
        'INH_Conf_Grade',
        'EMB_Conf_Grade',
        'PZA_Conf_Grade',
        'LEV_Conf_Grade',
        'MXF_Conf_Grade',
        'BDQ_Conf_Grade',
        'LZD_Conf_Grade',
        'CFZ_Conf_Grade',
        'DLM_Conf_Grade',
        'AMI_Conf_Grade',
        'STM_Conf_Grade',
        'ETH_Conf_Grade',
        'KAN_Conf_Grade',
        'CAP_Conf_Grade']
    drugs = {
        drug.split("_")[0]: {
            "R": set(),
            "U": set(),
            "S": set(),
            "F": set()} for drug in drug_columns}

    #Iterate over the catalogue
    for (index, row) in data.iterrows():
        garc = []
        #Pull out gene name, pos, ref and alt
        gene = row["gene_name"]
        if masks.get(gene) is None:
            #Cache the masks
            masks = {gene: get_masks(reference, gene)}
        pos = str(row["final_annotation.Position"])
        ref = row["final_annotation.ReferenceNucleotide"]
        alt = row["final_annotation.AlternativeNucleotide"]

        #Check for multiple positions defined within pos
        if len(pos.split(",")) > 1:
            #There is more than 1 position defined
            
            #Case where the number of positions matches the number of reference bases
            if len(pos.split(",")) == len(ref) == len(alt):
                pos = pos.split(",")
                for (p, r, a) in zip(pos, ref, alt):
                    garc += to_garc(reference, gene, int(p), r, a, masks)
            #Case where the number of positions is less than the number of reference bases
            elif len(pos.split(",")) < len(ref):
                index = 0
                pos = pos.split(",")
                if len(ref) != len(alt):
                    raise Exception("Different ref and alt lengths")
                for (r, a) in zip(ref, alt):
                    if r != a:
                        garc += to_garc(reference, gene, int(pos[index]), r, a, masks)
                        index += 1
            else:
                raise Exception("Weird format: "+ref+", "+alt)
        else:
            garc += to_garc(reference, gene, int(pos), ref, alt, masks)
        # if reference.genes[gene].reverse_complement:
        # print(index, garc, gene, pos, ref, alt)
        for drug in drug_columns:
            col = row[drug]
            drug = drug.split("_")[0]
            category = None
            if pd.isnull(col):
                continue
            if "1)" in col:
                # Resistance
                category = "R"
            elif "3)" in col:
                # Uncertain
                category = "U"
            elif "2)" in col or "4)" in col or "5)" in col or col == "Synonymous":
                # Not resistant
                category = "S"
            else:
                category = "F"

            for mutation in garc:
                drugs[drug][category].add(mutation)

    with open("output.csv", "w") as f:
        header = "GENBANK_REFERENCE,CATALOGUE_NAME,CATALOGUE_VERSION,CATALOGUE_GRAMMAR,PREDICTION_VALUES,DRUG,MUTATION,PREDICTION,SOURCE,EVIDENCE,OTHER\n"
        common_all = "NC_000962.3,WHO-UCN-GTB-PCI-2021.7,1.0,GARC1,RSUF,"
        f.write(header)
        resist = 0
        for drug in drugs.keys():
            print(drug)
            common = common_all + drug + ","
            for category in drugs[drug].keys():
                for mutation in sorted(drugs[drug][category]):
                    f.write(common + mutation + "," + category + ",{},{},{}\n")
            print({key: len(drugs[drug][key]) for key in drugs[drug].keys()})
            print({key: len(drugs[drug]["R"].intersection(drugs[drug][key])) 
            for key in drugs[drug].keys() 
            if key != "R" and len(drugs[drug]["R"].intersection(drugs[drug][key])) > 0 })
            synon_resist = len([m for m in drugs[drug]["R"] if "=" in m])
            print(synon_resist)
            print()
    # print(garc)
    print("Total indel: ", count)
    print("Total with SNPs: ", len([i for i in snp_count if len(i)>0]))
    print("Max SNPs: ", max([len(i) for i in snp_count]))
    print("Mean SNPs: ", sum([len(i) for i in snp_count])/len([i for i in snp_count if len(i)>0]))
    print("Median SNPs: ", sorted([len(i) for i in snp_count if len(i) > 0])[len([i for i in snp_count if len(i)>0])//2])
    print(sorted([len(i) for i in snp_count if len(i)>0]))

import functools
import re
import sys

import numpy
import pandas as pd

import gumpy


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
    offset = 0
    for (index, (r, a)) in enumerate(zip(ref, alt)):
        if r is None or a is None:
            offset += 1
            continue
        if r is not None and a is not None and r != a:
            if (pos + index - offset) - reference.genes_lookup[gene]["start"] <= 0:
                #Past the end of the gene so just return
                return mutations
            if reference.genes_lookup[gene]["end"] - (pos + index - offset) < 0 or reference.genes[gene].codes_protein == False:
                #Promoter or non-coding so return the difference in nucleotides
                r,a  = gumpy.Gene._complement([r, a])
                mutations.append(gene + "@" + r + str(reference.genes_lookup[gene]["end"] - (pos + index - offset)) + a)
            else:
                ref_seq[pos + index - offset - 1] = a
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
    offset = 0
    for (index, (r, a)) in enumerate(zip(ref, alt)):
        if r is None or a is None:
            offset += 1
            continue
        if r is not None and a is not None and r != a:
            if reference.genes_lookup[gene]["end"] - (pos + index - offset) <= 0:
                #Past the end of the gene so just return
                return mutations
            if (pos + index - offset) - reference.genes_lookup[gene]["start"] < 0 or reference.genes[gene].codes_protein == False:
                #Promoter or non-coding so return the difference in nucleotides
                mutations.append(gene + "@" + r + str((pos + index - offset) - reference.genes_lookup[gene]["start"]) + a)
            else:
                ref_seq[pos + index - offset - 1] = a
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


def del_calls(reference, gene, pos, ref, alt, masks, rev_comp=False):
    '''Deal with del calls. Attempts to identify dels mid-sequence. 
        If a repeated base is deleted (aaa->aa), it is assumed that the first base is deleted.
    
    Args:
        reference (gumpy.Genome): Reference genome object
        gene (str): Gene name
        pos (int): Genome position
        ref (list): Reference bases
        alt (list): Alternative bases
        masks (dict): Dictionary of gene_name->(stacked_mask, mask)
        rev_comp (bool, optional): Flag to determine if reverse complement. Defaults to False
    Returns:
        list(str): List of mutations in GARC
    '''
    #Del has len(alt) < len(ref)
    del_len = len(ref) - len(alt)
    current = None
    current_snps = 999
    start = 0
    ref1 = list(ref)
    #Iterate through the positions at which the ins could occur, checking which has the lowest overall SNPs
    for x in range(len(alt)+1):
        alt1 = [alt[i] for i in range(x)]+[None for i in range(del_len)]+[alt[i] for i in range(x, len(alt))]
        if snp_number(ref1, alt1) < current_snps:
            current = alt1
            current_snps = snp_number(ref1, alt1)
            start = x
    #Position with the best SNPs is the best position for the ins
    seq = [ref[i] for i in range(len(current)) if current[i] is None]
    if rev_comp:
        p = reference.genes_lookup[gene]["end"] - (pos + start)
        r = ''.join(gumpy.Gene._complement(seq))
        snp = rev_comp_snp(reference, gene, pos, ref, current, masks)
        if p > reference.genes_lookup[gene]["end"]:
            #Del happened past the 3' end of the gene so ignore it
            return snp
    else:
        p = pos - reference.genes_lookup[gene]["start"] + start
        r = ''.join(seq)
        snp = snps(reference, gene, pos, ref, current, masks)
        if p > reference.genes_lookup[gene]["end"]:
            #If the del happened past the 3' end of the gene, ignore it
            return snp
    return snp + [gene + "@" + str(p) + "_del_" + r]

def snp_number(ref, alt):       
    '''Helper function to find the SNP distance between two arrays, ignoring None values

    Args:
        ref (list): List of bases
        alt (list): List of bases

    Returns:
        int: SNP distance ignoring None values
    '''    
    snps = 0
    for (a, b) in zip(ref, alt):
        if a is not None and b is not None and a != b:
            snps += 1
    return snps

def ins_calls(reference, gene, pos, ref, alt, masks, rev_comp=False):
    '''Deal with ins calls. Attempts to detect mid-sequence insertions.
        If a repeated base has an insertion, it is assumed that the insertion occured at first base (`aa`->`aaa` infers ins @ seq[0])

    Args:
        reference (gumpy.Genome): Reference Genome object
        gene (str): Gene name
        pos (int): Genome index
        ref (str): Reference bases
        alt (str): Alternative bases
        masks (dict): Dictionary of gene_name->(stacked_mask, mask)
        rev_comp (bool, optional): Flag to show if the gene is reverse complement. Defaults to False
    Returns:
        list(str): List of mutations in GARC
    '''    
    #Ins has len(ref) < len(alt)
    ins_len = len(alt) - len(ref)
    current = None
    #Arbitrarily high SNPs so it can only decrease
    current_snps = 999
    start = 0
    alt1 = list(alt)
    #Iterate through the positions at which the ins could occur, checking which has the lowest overall SNPs
    for x in range(len(ref)+1):
        ref1 = [ref[i] for i in range(x)]+[None for i in range(ins_len)]+[ref[i] for i in range(x, len(ref))]
        if snp_number(ref1, alt1) < current_snps:
            current = ref1
            current_snps = snp_number(ref1, alt1)
            start = x
    #Position with the best SNPs is the best position for the ins
    seq = [alt[i] for i in range(len(current)) if current[i] is None]
    if rev_comp:
        p = reference.genes_lookup[gene]["end"] - (pos + start)
        a = ''.join(gumpy.Gene._complement(seq))
        snp = rev_comp_snp(reference, gene, pos, current, alt, masks)
        if p > reference.genes_lookup[gene]["start"]:
            #Past the 3' end so ignore
            return snp
    else:
        p = pos - reference.genes_lookup[gene]["start"] + start
        a = ''.join(seq)
        snp = snps(reference, gene, pos, current, alt, masks)
        if p > reference.genes_lookup[gene]["end"]:
            #Past the 3' end so ignore
            return snp
    return snp + [gene + "@" + str(p) + "_ins_" + a]


def to_garc(reference, gene, pos, ref, alt, masks):
    '''Convert to GARC

    Args:
        reference (gumpy.Genome): Reference genome object
        pos (int): Genome index
        ref (str): Reference base(s)
        alt (str): Mutant base(s)
        masks (dict): Dictionary of gene_name->(stacked_mask, mask)
    Returns:
        list(str): List of mutations in GARC
    '''
    rev_comp = reference.genes[gene].reverse_complement
    if len(ref) == len(alt):
        if rev_comp:
            return rev_comp_snp(reference, gene, pos, ref, alt, masks)
        else:
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
        In order to make some sense of this, the indel is selected based on where in the sequence it causes the least SNPs. If there are
        repeating sequences as detailed above, the first one is selected. This is not a perfect solution, but allows some form of standardisation
        with the least possible mutations from a single row.
        '''
        #Del
        return del_calls(reference, gene, pos, ref, alt, masks, rev_comp=rev_comp)
    elif len(ref) < len(alt):
        #Ins
        return ins_calls(reference, gene, pos, ref, alt, masks, rev_comp=rev_comp)


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

def generalise(mutation):
    '''Takes a mutation in GARC and produces a list of generalised mutations

    Args:
        mutation (str): Mutation
    Returns:
        list(str): List of generalised versions of mutations
    '''
    if re.compile(
            r"([a-zA-Z0-9]+)@((-?\d+=)|([acgtA-Z!])(-?\d+)([actgA-Z!]))").fullmatch(mutation):
        # SNP so no generalisation required
        return [mutation]
    else:
        if "ins" in mutation:
            # Of the form gene@pos_ins_bases
            gene_name = mutation.split("@")[0]
            pos, ins, bases = mutation.split("@")[1].split("_")
            if len(bases) % 3 != 0:
                # Frame shift
                fs = [gene_name + "@" + pos + "_fs"]
            else:
                fs = []
            items =  fs + [
                mutation,
                gene_name + "@" + pos + "_indel",
                gene_name + "@" + pos + "_ins",
                gene_name + "@" + pos + "_ins_" + str(len(bases)),
            ]
        elif "del" in mutation:
            # Of the form gene@pos_del_bases
            gene_name = mutation.split("@")[0]
            pos, del_, bases = mutation.split("@")[1].split("_")
            if len(bases) % 3 != 0:
                # Frame shift
                fs = [gene_name + "@" + pos + "_fs"]
            else:
                fs = []
            items =  fs + [
                mutation,
                gene_name + "@" + pos + "_indel",
                gene_name + "@" + pos + "_del",
                gene_name + "@" + pos + "_del_" + str(len(bases)),
            ]
        # items = [i for i in items if i not in seen]
        return items

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
            # "F": set()
            } for drug in drug_columns}
    #Iterate over the catalogue
    for (index, row) in data.iterrows():
        garc = []
        #Pull out gene name, pos, ref and alt
        gene = row["gene_name"]
        if masks.get(gene) is None:
            #Cache the masks
            masks = {gene: get_masks(reference, gene)}
        pos = str(row["final_annotation.Position"])#Cast to a str for str.split(',')
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
        if True in [m in 
                    {'rpoB@1296_ins_ttc', 'pncA@317_del_t', 'pncA@394_del_gctggtgta', 'pncA@391_ins_gg', 
                    'pncA@-4_del_g', 'pncA@456_ins_c', 'pncA@389_del_tgta', 'eis@c-13t', 'gid@102_del_g', 
                    'gid@351_del_g', 'ethA@110_del_a'} 
                    for m in garc]:
            print(ref)
            print(alt)
            print(garc)
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
                print(drug)
                print(category)
            print()
        if "generalise" in sys.argv:
            #Produce more general cases such as gene@pos_indel if a flag is given
            garc = functools.reduce(lambda x, y: x + y, [generalise(g) for g in garc], [])
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
            # else:
                # category = "F"

            for mutation in garc:
                drugs[drug][category].add(mutation)


            
            

    with open("output.csv", "w") as f:
        header = "GENBANK_REFERENCE,CATALOGUE_NAME,CATALOGUE_VERSION,CATALOGUE_GRAMMAR,PREDICTION_VALUES,DRUG,MUTATION,PREDICTION,SOURCE,EVIDENCE,OTHER\n"
        common_all = "NC_000962.3,WHO-UCN-GTB-PCI-2021.7,1.0,GARC1,RUS,"
        f.write(header)
        resist = 0
        for drug in drugs.keys():
            print(drug)
            common = common_all + drug + ","
            for category in drugs[drug].keys():
                #As there are some mutations which are R and another category,
                #just assume that they belong to R so piezo can parse the catalogue.
                m = drugs[drug][category]
                for c in drugs[drug].keys():
                    if c != category:
                        m = m.difference(drugs[drug][c])
                mutations = sorted(list(m))
                # if category != "R":
                #     mutations = sorted(list(drugs[drug][category].difference(drugs[drug]["R"])))
                # else:
                #     mutations = sorted(list(drugs[drug][category]))
                for mutation in mutations:
                    f.write(common + mutation + "," + category + ",{},{},{}\n")
            # print({key: len(drugs[drug][key]) for key in drugs[drug].keys()})
            print({key: drugs[drug]["R"].intersection(drugs[drug][key])
            for key in drugs[drug].keys() 
            if key != "R" and len(drugs[drug]["R"].intersection(drugs[drug][key])) > 0 })
            synon_resist = len([m for m in drugs[drug]["R"] if "=" in m])
            # print(synon_resist)
            print()

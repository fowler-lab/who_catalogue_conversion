import re, copy

import pandas as pd
from Bio import SeqIO

import gumpy

'''Version of translate which uses the raw postion, reference and alt values to bypass the problems identified with HGVS usage.
Requires gumpy and BioPython
'''
count = 0
count2 = 0

def parse_who_catalog(filename):
    '''Parses the WHO TB catalog

    Args:
        filename (str): Path to the WHO catalog
    Returns:
        pd.dataframe: Dataframe containing the mutations
    '''
    df = pd.read_excel(filename, sheet_name="Genome_indices")
    return df


def get_snps(gene, pos, end, ref, alt):
    '''Get the SNPs from given pos, ref and alt values

    Args:
        gene (str): Gene name
        pos (int): Gene index of the item
        end (int): Maximum gene index
        ref (str): Reference nucleotide(s)
        alt (str): Mutated nucleotide(s)

    Returns:
        list(str): List of SNP mutations in GARC
    '''
    if pos > end:
        return []
    if len(ref) == 1 and len(alt) == 1:
        return [gene + "@" + ref.lower() + str(pos) + alt.lower()]
    else:
        mutations = []
        for (i, (r, a)) in enumerate(zip(ref, alt)):
            if r != a:
                # There is a SNP here
                if pos + i > end:
                    return mutations
                mutations.append(gene + "@" + r.lower() +
                                 str(pos + i) + a.lower())
    return mutations


def to_garc(gene, pos, end, ref, alt):
    '''Converts a given mutation to GARC

    Args:
        gene (str): Gene name
        pos (int): Gene index
        end (int): Maximum index for this gene. Required as catalogue specifies bases past the 3' end
        ref (str): Reference nucleotide(s)
        alt (str): Mutant nucleotide(s)

    Returns:
        list(str): List of the mutations in GARC
    '''
    if len(ref) == len(alt):
        # SNP(s)
        return get_snps(gene, pos, end, ref, alt)
    elif len(ref) > len(alt):
        # Del (and SNPs)
        # Check for SNPs within the shared bases
        snps = get_snps(gene, pos, end, ref, alt)
        # Get the dels
        dels = [gene + "@" + str(pos + len(alt)) + "_del_" + ref[len(alt):]]
        return snps + dels
    elif len(ref) < len(alt):
        # Ins (and SNPs)
        # Check for SNPs within shared bases
        snps = get_snps(gene, pos, end, ref, alt)
        # Get the ins
        ins = [gene + "@" + str(pos + len(alt)) + "_ins_" + alt[len(ref):]]
        return snps + ins

def to_aa(reference, mutation):
    '''Converts mutation from nucleotide to amino acid as appropriate

    Args:
        reference (gumpy.Genome): Reference genome object
        mutations (str): Mutation in GARC
    '''
    if re.compile(r"([a-zA-Z0-9]+)@([acgt])(\d+)([acgt])").fullmatch(mutation):
        # SNP format
        gene, ref, pos, alt = re.compile(
            r"([a-zA-Z0-9]+)@([acgt])(\d+)([acgt])").findall(mutation)[0]
    else:
        # Just return frame shifts and indels as they make weird changes
        return mutation
    gene_name = gene
    gene = reference.genes[gene]
    if gene.reverse_complement:
        global count
        count += 1

        # if gene.nucleotide_sequence[gene.is_cds][int(pos)] != ref:
        #     print("##",mutation, gene.name, pos, gene.nucleotide_sequence[gene.is_cds][int(pos)], ref, alt, gene.reverse_complement)
        #Undo the reverse complementing
        gene.nucleotide_sequence=gene._complement(gene.nucleotide_sequence[::-1])
        gene.index=gene.index[::-1]
        gene.nucleotide_number=gene.nucleotide_number[::-1]
        gene.is_cds=gene.is_cds[::-1]
        gene.is_promoter=gene.is_promoter[::-1]
        gene.is_indel=gene.is_indel[::-1]
        gene.indel_length=gene.indel_length[::-1]
        # pos = int(pos) - 1
    #Rebuild a gene object with the alt value
    if gene.nucleotide_sequence[gene.is_cds][int(pos)] != ref:
        print(mutation, gene.name, pos, gene.nucleotide_sequence[gene.is_cds][int(pos)], ref, gene.reverse_complement)
        # if reference.nucleotide_sequence[reference.genes_lookup[gene_name]["start"] + int(pos) -1] != ref:
        #     print(reference.nucleotide_sequence[reference.genes_lookup[gene_name]["start"] + int(pos) - 1], int(pos), reference.genes_lookup[gene_name]["start"])
        global count2
        count2 += 1
        print()
    gene.nucleotide_sequence[gene.is_cds][int(pos)] = alt
    g = gumpy.Gene(
        name=gene.name,
        nucleotide_sequence = gene.nucleotide_sequence,
        index = gene.index,
        nucleotide_number = gene.nucleotide_number,
        is_cds = gene.is_cds,
        is_promoter = gene.is_promoter,
        is_indel = gene.is_indel,
        indel_length = gene.indel_length,
        reverse_complement = gene.reverse_complement,
        codes_protein = gene.codes_protein,
        feature_type = gene.feature_type
    )
    aa_index = 3

    # # if gene.codes_protein:
    # #     # Get reference AA
    # #     if gene.reverse_complement:
    # #         triplets = gene.triplet_number
    # #         gene.nucleotide_sequence = gene.nucleotide_sequence[::-1]
    # #         is_cds = gene.is_cds[::-1]
    # #         ref = gene._complement([ref])[0]
    # #         alt = gene._complement([alt])[0]
    # #         pos = (int(pos) - len(gene.is_cds)) * -1
    # #     else:
    # #         triplets = gene.triplet_number
    # #         is_cds = gene.is_cds


    # #     aa_index = triplets[int(pos)] - 1

    #     if gene.nucleotide_sequence[is_cds][int(pos)] != ref:
    #         print(gene_name, gene.nucleotide_sequence[is_cds][int(pos)], ref, pos, aa_index)
    #     gene.nucleotide_sequence[is_cds][int(pos)] = alt
    #     ref = gene.amino_acid_sequence[aa_index]
    #     gene._translate_sequence()

    #     alt = gene.amino_acid_sequence[aa_index]
        # print(mutation, aa_index, ref, alt)
        # print()
    if ref != alt:
        return gene_name + "@" + ref + str(aa_index + 1) + alt
    else:
        # Synonymous mutation
        return gene_name + "@" + str(aa_index + 1) + "="

if __name__ == "__main__":
    # Load the reference genome for gathering gene indices
    genbank = SeqIO.read("NC_000962.3.gbk", "genbank")
    # reference_genome = gumpy.Genome("NC_000962.3.gbk", multithreaded=True)
    # reference_genome.save("reference.json.gz", compression_level=1)
    reference_genome = gumpy.Genome.load("reference.json.gz")
    gene_starts = {}
    gene_ends = {}
    codes_protein = {}
    for feature in genbank.features:
        if 'gene' in feature.qualifiers.keys():
            gene_name = feature.qualifiers['gene'][0]
        elif 'locus_tag' in feature.qualifiers.keys():
            gene_name = feature.qualifiers['locus_tag'][0]
        else:
            continue
        gene_starts[gene_name] = int(feature.location.start) + 1
        gene_ends[gene_name] = int(feature.location.end)
        codes_protein[gene_name] = feature.type != "rRNA"

    data = parse_who_catalog("WHO-UCN-GTB-PCI-2021.7-eng.xlsx")
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
            "?": set()} for drug in drug_columns}
    incorrect = 0
    for (index, row) in data.iterrows():
        garc = []
        # Get pos, ref, alt, gene_name
        gene = row["gene_name"]
        if "," in str(row["final_annotation.Position"]):
            # Mutliple indcies are included
            if len(str(row["final_annotation.ReferenceNucleotide"])) == len(str(
                    row["final_annotation.AlternativeNucleotide"])) == str(row["final_annotation.Position"]).count(","):
                # There are the same number of indices as there are bases
                # So each index corresponds to each base
                for (pos, ref, alt) in zip(row["final_annotation.Position"].split(","), list(
                        row["final_annotation.ReferenceNucleotide"]), list(row["final_annotation.AlternativeNucleotide"])):
                    garc += to_garc(gene, int(pos) -
                                    gene_starts[gene], gene_ends[gene] -
                                    gene_starts[gene], ref, alt)
            else:
                # Weird edge case where the number of nucleotides is different than the number of indices
                # This means that we have to determine which bases are
                # referenced by the indices based on which are different
                index = 0
                for (ref, alt) in zip(list(row["final_annotation.ReferenceNucleotide"]), list(
                        row["final_annotation.AlternativeNucleotide"])):
                    if ref != alt:
                        pos = int(
                            row["final_annotation.Position"].split(",")[index])
                        garc += to_garc(gene, pos -
                                        gene_starts[gene], gene_ends[gene] -
                                        gene_starts[gene], ref, alt)
                        index += 1
        else:
            # Single index, so continue as expected
            pos = int(row["final_annotation.Position"]) - gene_starts[gene]
            end = gene_ends[gene] - gene_starts[gene]
            ref = row["final_annotation.ReferenceNucleotide"]
            alt = row["final_annotation.AlternativeNucleotide"]
            garc += to_garc(gene, pos, end, ref, alt)

        # Convert to AA as appropriate if gene codes protein
        if codes_protein[gene]:
            aa_garc = []
            for mutation in garc:
                if "-" in mutation:
                    # Promoter so skip
                    aa_garc.append(mutation)
                else:
                    mutation = to_aa(reference_genome, mutation)
                    aa_garc.append(mutation)
            # Remove duplicates
            garc = sorted(list(set([f"{i}" for i in aa_garc])))
        else:
            # Remove duplicates
            garc = sorted(list(set([f"{i}" for i in garc])))
        # print(garc)

        for drug in drug_columns:
            col = row[drug]
            drug = drug.split("_")[0]
            category = None
            if pd.isnull(col):
                continue
            if "1)" in col or "2)" in col:
                # Resistance
                category = "R"
            elif "3)" in col:
                # Uncertain
                category = "U"
            elif "4)" in col or "5)" in col:
                # Not resistant
                category = "S"
            else:
                category = "?"

            # Check for indels first
            if len([m for m in garc if "ins" in m or "del" in m]) > 0:
                for m in [m for m in garc if "ins" in m or "del" in m]:
                    drugs[drug][category].add(m)
            else:
                # No indels, so just add all
                for mutation in garc:
                    drugs[drug][category].add(mutation)

    with open("output.csv", "w") as f:
        header = "GENBANK_REFERENCE,CATALOGUE_NAME,CATALOGUE_VERSION,CATALOGUE_GRAMMAR,PREDICTION_VALUES,DRUG,MUTATION,PREDICTION,SOURCE,EVIDENCE,OTHER\n"
        common_all = "NC_000962.3,WHO-UCN-GTB-PCI-2021.7,1.0,GARC1,RSU?,"
        f.write(header)
        resist = 0
        for drug in drugs.keys():
            print(drug)
            common = common_all + drug + ","
            for category in drugs[drug].keys():
                for mutation in sorted(drugs[drug][category]):
                    f.write(common + mutation + "," + category + ",{},{},{}\n")
            print({key: len(drugs[drug][key]) for key in drugs[drug].keys()})
            synon_resist = len([m for m in drugs[drug]["R"] if "=" in m])
            print(synon_resist)
            resist += synon_resist
            print()
    # print(resist)
    print(count)
    print(count2)
    # print(garc)
    # print(len([i for i in garc if "fs" in i]))

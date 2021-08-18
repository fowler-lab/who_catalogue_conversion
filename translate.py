import functools
import re
from collections import defaultdict

import pandas as pd

'''
Problems:
    Use of x_ydelABCinsBBC:
        This is an undocumented use of the grammar (it stipulates use of EITHER del or ins)
        Most of these items could be denoted as SNPs -
        Some of the sequences denoted in this manner are 50+ bases long, some for only 4 SNPs
        Some items are included which cross the 3' end which seems unhelpful?

    Protein mutations are not all correct, especially for long x_ydelABCinsDBC
        e.g use of 'p.=' despite having 3 different amino acids within the indel
    Alternatively there are several protein mutations denoted not as 'p.=' (synonymous protein mutation),
        but are infact synonymous e.g 'p.Gln431Gln'

'''


class Mutation(object):
    def __init__(self, gene, pos, end, ref, alt, ins, del_):
        self.gene = gene
        self.pos = pos
        self.end = end
        self.ref = ref
        self.alt = alt
        self.ins = ins
        self.del_ = del_

    def __to_garc(self, gene, pos, ref, alt, ins, del_):
        '''Convert a single mutation to GARC

        Args:
            gene (str): Gene name
            pos (int): Position index
            ref (str): Reference nucleotide
            alt (str): Alternative nucleotide
            ins (str): Insertion nucleotide sequence
            del_ (str): Deletion nucleotide sequence

        Returns:
            str: Mutation in GARC
        '''
        gene = gene + "@"
        # Case for SNP
        if ins is None and del_ is None:
            return gene + ref.lower() + str(pos) + alt.lower()
        # Case for frame shift
        if (ins is not None and len(ins) % 3 != 0) or (
                del_ is not None and len(del_) % 3 != 0):
            return gene + str(pos) + "_fs"
        # Case for ins
        if ins is not None and del_ is None:
            return gene + str(pos) + "_ins_" + ins.lower()
        # Case for del
        if ins is None and del_ is not None:
            return gene + str(pos) + "_del_" + del_.lower()

    def to_garc(self):
        '''Converts this to GARC. As mutliple positions can be specifed, returns a list of GARC mutations

        Returns:
            [str]: List of strings of GARC mutations
        '''
        if self.end is None:
            # Single mutation
            if self.ins is not None and self.del_ is not None:
                # Indel so use both an ins and del
                return sorted(list(set([
                    self.__to_garc(
                        self.gene,
                        self.pos,
                        self.ref,
                        self.alt,
                        self.ins,
                        None),
                    self.__to_garc(
                        self.gene,
                        self.pos,
                        self.ref,
                        self.alt,
                        None,
                        self.del_),
                ])))
            return [self.__to_garc(self.gene, self.pos,
                                   self.ref, self.alt, self.ins, self.del_)]
        else:
            # Deal with ranges of bases...
            # Case for *x
            if "*" in self.end:
                # Adjust to ignore mutations past the 3' end
                ignore_length = int(self.end.replace("*", ""))
                if self.ins is not None and self.del_ is not None:
                    # Indel so check if the ins would be in gene range
                    if len(self.del_) != len(self.ins):
                        # There is actually an indel here
                        if len(self.del_) > len(self.ins):
                            # There is more deleted than inserted, so check for
                            # changes required for ins
                            if len(self.ins) < ignore_length:
                                # Insertion is smaller than ignore length so
                                # leave ins alone
                                self.del_ = self.del_[:-ignore_length]
                            else:
                                self.ins = self.ins[:-ignore_length]
                                self.del_ = self.del_[:-ignore_length]
                        else:
                            self.ins = self.ins[:-ignore_length]
                            self.del_ = self.del_[:-ignore_length]
                    else:
                        self.ins = self.ins[:-ignore_length]
                        self.del_ = self.del_[:-ignore_length]
                elif self.ins is not None:
                    self.ins = self.ins[:-ignore_length]
                elif self.del_ is not None:
                    self.del_ = self.del_[:-ignore_length]
            # If an indel, replace with ins and del
            if self.ins is not None and self.del_ is not None:
                # Check if an indel is actually an indel as some indels are
                # actually SNPs
                if len(self.ins) == len(self.del_):
                    garc = []
                    # Same length ins/del so can actually be given as SNPs
                    for (index, (i, d)) in enumerate(zip(self.ins, self.del_)):
                        if i != d:
                            # SNP
                            garc.append(
                                self.__to_garc(
                                    self.gene, int(
                                        self.pos) + index, d, i, None, None))
                    return sorted(list(set(garc)))
                elif len(self.ins) != len(self.del_) and (len(self.ins) - len(self.del_)) % 3 != 0:
                    # Frame shifting mutation
                    return [self.gene + "@" + str(self.pos) + "_fs"]
                return sorted(list(set([
                    self.__to_garc(
                        self.gene,
                        self.pos,
                        self.ref,
                        self.alt,
                        self.ins,
                        None),
                    self.__to_garc(
                        self.gene,
                        self.pos,
                        self.ref,
                        self.alt,
                        None,
                        self.del_),
                ])))
            # Otherwise, just treat it as a usual ins/del
            else:
                return [self.__to_garc(
                    self.gene, self.pos, self.ref, self.alt, self.ins, self.del_)]


def parse_who_catalog(filename):
    '''Parses the WHO TB catalog

    Args:
        filename (str): Path to the WHO catalog
    Returns:
        pd.dataframe: Dataframe containing the mutations
    '''
    df = pd.read_excel(filename, sheet_name="Genome_indices")
    return df


def parse_mutation(mutation, gene):
    '''Parse a mutation to a data structure

    Args:
        mutation (str): Mutation in HGVS
        gene (str): Name of the gene
    Returns:
        Mutation: A Mutation object to use as an itermediary store
    '''
    # Ensure the row is of the correct format
    pattern = re.compile(r"""
                    ([cgmnopr])\. #Mutation type is denoted using a single letter - group 0
                    ([\*\-\+]?\d+)_?([\*\-\+]?\d+)? #Position (range)? - group 1, 2
                    ( #Picking a mutation type - group 3
                        ([ACTG]>[ACTG]) #SNP - group 4
                        |(ins[ACTG]+) #Ins - group 5
                        |(del[ACTG]+ins[ACTG]+) #Indel - group 6
                        |(dup[ACTG]+) #Duplicate - group 7
                        |(del[ACTG]+) #Del - group 8
                    )""", re.VERBOSE)
    if pattern.fullmatch(mutation) is None:
        raise Exception("Not in correct format: " + mutation)
    else:
        # The mutation is valid so use this pattern to pull out groups
        groups = pattern.findall(mutation)[0]
        mut_type = groups[0]
        pos = groups[1]
        if groups[2] != "":
            # There is a secondary position, so a range was given
            end = groups[2]
        else:
            end = None

        # Only one of these should exist
        if groups[4] != "":
            # SNP
            ref, alt = groups[4].split(">")
            ins = del_ = None
        elif groups[5] != "":
            # Ins
            ins = groups[5].replace("ins", "")
            del_ = ref = alt = None
        elif groups[6] != "":
            # Indel
            del_, ins = groups[6].replace("del", "").split("ins")
            ref = alt = None
        elif groups[7] != "":
            # Dup which converts to an ins
            ins = groups[7].replace("dup", "")
            del_ = ref = alt = None
        elif groups[8] != "":
            # Del
            del_ = groups[8].replace("del", "")
            ins = ref = alt = None
        else:
            raise Exception("No groups detected: " + mutation)
        return Mutation(gene, pos, end, ref, alt, ins, del_)


if __name__ == "__main__":
    # Read the xslx
    data = parse_who_catalog("WHO-UCN-GTB-PCI-2021.7-eng.xlsx")
    # Pull out the mutations in HGVS
    # mutations = data["gene_name", "final_annotation.TentativeHGVSNucleotidicAnnotation"]
    fs = 0
    fs_genes = defaultdict(int)
    for (index, row) in data.iterrows():
        gene = row["gene_name"]
        mutation = row["final_annotation.TentativeHGVSNucleotidicAnnotation"]
        if not pd.isnull(mutation):
            for m in mutation.split(","):
                # if "*" in m:
                #     print()
                #     print(gene+"@"+m)
                m = parse_mutation(m, gene)
                garc = m.to_garc()
                if functools.reduce(lambda x, y: x or y, [
                                    "fs" in i for i in garc], False):
                    print(mutation, garc)
                    fs += 1
                    fs_genes[gene] += 1
    print(fs)
    print(fs_genes)

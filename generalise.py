import pandas as pd
from collections import defaultdict, Counter
import re
'''
Due to how few mutations there are at several points, especially when considering just
R and S predictions, there is flimsy/no basis for creating general rules such as `rpoB@142_del_* R`

This would make it more usable, but when for most indels conferring either R or S, there is at most 2 mutations,
this seems a little presumptive...
Naively converting all indels to be more general (i.e ins_ac->(ins_2, ins, indel, fs)) would naively produce such rules, 
but there would likely be collisons such as when `13_ins_t R` and `13_del_t S` would both produce `13_indel` and `13_fs`,
which would have to be removed to work with piezo - and the resulting `13_ins_1` is only based on a single `ins_t` so `ins_g` 
may not have the same effect

Similarly for SNP mutations, there are no cases where every SNP mutation has been detailed with the same prediction, 
and often there are some mutations at a single gene index where certain amino acids infer R, and some infer S - 
preventing a neat `rpoB@A23* S` style rule
'''


if __name__ == "__main__":
    #Load the output csv
    data = pd.read_csv("output.csv")
    snps = defaultdict(set)
    refs = defaultdict(Counter)
    '''GENBANK_REFERENCE,CATALOGUE_NAME,CATALOGUE_VERSION,CATALOGUE_GRAMMAR,PREDICTION_VALUES,DRUG,MUTATION,PREDICTION,SOURCE,EVIDENCE,OTHER'''
    for (index, row) in data.iterrows():
        if index == 0:
            continue
        mutation = row["MUTATION"]
        drug = row["DRUG"]
        pred = row["PREDICTION"]
        #Get the SNPs as indels are already generalised
        snp = re.compile(r"([a-zA-Z0-9_]+)@([acgtzxA-Z])(-?[0-9]+)([acgtzxA-Z])")
        synon = re.compile(r"([a-zA-Z0-9_]+)@([0-9]+)=")
        indel = re.compile(r"([a-zA-Z0-9_]+)@(-?[0-9]+)_(ins|del)_([acgtzx]+)")
        if snp.fullmatch(mutation):
            #The mutation is a SNP
            gene, ref, pos, alt = snp.fullmatch(mutation).groups()
        elif synon.fullmatch(mutation):
            gene, pos = synon.fullmatch(mutation).groups()
            ref = "="
            alt = "="
        elif indel.fullmatch(mutation):
            gene, pos, type_, bases = indel.fullmatch(mutation).groups()
            if type_ == "ins":
                ref = "ins"
                alt = "ins_"+bases
            elif type_ == "del":
                ref = bases
                alt = "del_"+bases
        # if synon.fullmatch(mutation) or snp.fullmatch(mutation) or indel.fullmatch(mutation):
        if indel.fullmatch(mutation):
            refs[(gene, int(pos))].update([ref])
            # if refs.get((gene, int(pos))) is not None and refs.get((gene, int(pos)))!= ref:
            #     print("???", gene, pos, refs.get((gene, int(pos))), ref)
            # else:
            #     refs[(gene, int(pos))] = ref
            snps[( drug, pred)].add((gene, int(pos), alt))
    by_drug = defaultdict(dict)
    for key in snps.keys():
        by_drug[key[0]][key[1]] = snps[key]
        # print(key, len(snps[key]))
        # print()
    for drug in by_drug.keys():
        print(drug, by_drug[drug].keys())
        interesting = {pos: (gene, pos, alt) for (gene, pos, alt) in by_drug[drug].get('R', set()).union(by_drug[drug].get('S', set())).union(by_drug[drug].get('U', set()))}
        for pos in sorted(list(interesting.keys())):
            r = sorted(list([i for i in by_drug[drug].get('R', []) if i[0]==interesting[pos][0] and i[1]==pos] ))
            s = sorted(list([i for i in by_drug[drug].get('S', [])if i[0]==interesting[pos][0] and i[1]==pos] ))
            u = sorted(list([i for i in by_drug[drug].get('U', [])if i[0]==interesting[pos][0] and i[1]==pos] ))
            if sum([1 for i in [r, s, u] if len(i) > 0]) > 1:
                print(drug, interesting[pos])
                print("R", r)
                print("S", s)
                print("U", u)
                print()
            
        # for (gene, pos, alt) in sorted(list(set([x for x in by_drug[drug].get('R', set())]))):
        #     print(drug, gene, list(refs[(gene, pos)].keys())[0], pos, alt, "R")
        # for (gene, pos, alt) in sorted(list(set([x for x in by_drug[drug].get('S', set())]))):
        #     print(drug, gene, list(refs[(gene, pos)].keys())[0], pos, alt, "S")
        # print(drug, len(by_drug[drug].get('R', set()).intersection(by_drug[drug].get('S', set())) ))
        # for prediction in by_drug[drug].keys():
        #     print(drug, prediction, len(by_drug[drug][prediction]))

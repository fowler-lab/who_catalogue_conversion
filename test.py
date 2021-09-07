import pandas as pd
from functools import reduce
from collections import Counter

'''Testing some aspects of the catalogue...
'''
def parse_who_catalog(filename):
    '''Parses the WHO TB catalog

    Args:
        filename (str): Path to the WHO catalog
    Returns:
        pd.dataframe: Dataframe containing the mutations
    '''
    df = pd.read_excel(filename, sheet_name="Genome_indices")
    return df

if __name__ == "__main__":
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
    
    indelsR = []
    indelsS = []
    fsR = []
    fsS = []
    
    for (index, row) in data.iterrows():
        gene = row["gene_name"]
        pos = str(row["final_annotation.Position"])
        ref = row["final_annotation.ReferenceNucleotide"]
        alt = row["final_annotation.AlternativeNucleotide"]
        confs = ["" if pd.isnull(row[drug]) else row[drug] for drug in drug_columns]
        R = reduce(lambda x, y: x or y, ["1)" in c for c in confs], False)
        if R:
            #There was either assoc R or not assoc R
            #Check for indels
            if len(ref) != len(alt):
                indelsR.append([gene, pos, ref, alt, confs])
                if (len(ref) - len(alt)) % 3 != 0:
                    #Frame shift
                    fsR.append([gene, pos, ref, alt, confs])
        S = reduce(lambda x, y: x or y, ["5)" in c for c in confs], False)
        if S:
            #There was either assoc R or not assoc R
            #Check for indels
            if len(ref) != len(alt):
                indelsS.append([gene, pos, ref, alt, confs])
                if (len(ref) - len(alt)) % 3 != 0:
                    #Frame shift
                    fsS.append([gene, pos, ref, alt, confs])
        if len(ref) != len(alt):
            if "aaa" in ref and "aaa" not in alt and ("aa" in alt):
                print(ref, alt, R, S)
            if "ccc" in ref and "ccc" not in alt and ("cc" in alt):
                print(ref, alt, R, S)
            if "ggg" in ref and "ggg" not in alt and ("gg" in alt):
                print(ref, alt, R, S)
            if "ttt" in ref and "ttt" not in alt and ("tt" in alt):
                print(ref, alt, R, S)

    print(indelsR)
    print(len(indelsR))
    print(sorted(list(set([i[0] for i in indelsR]))))
    print(Counter([i[0] for i in indelsR]))
    print()
    print(indelsS)
    print(len(indelsS))
    print(sorted(list(set([i[0] for i in indelsS]))))
    print(Counter([i[0] for i in indelsS]))
    print()
    # for i in [i for i in indelsR+indelsS if i[0] == "rpoB"]:
    #     print(i[:4], i[4])
    # print(set([i[0] for i in indelsR]).intersection(set([i[0] for i in indelsS])))
    for i in indelsR+indelsS:
        print(i[0], i[1], i[2], i[3], sep="\t")

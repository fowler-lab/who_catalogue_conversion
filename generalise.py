import re
import copy
'''
Utilises the specific mutations generated with `parse.py` to produce a more generalised list of mutations
This is more helpful for piezo
'''


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
            return fs + [
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
            return fs + [
                mutation,
                gene_name + "@" + pos + "_indel",
                gene_name + "@" + pos + "_del",
                gene_name + "@" + pos + "_del_" + str(len(bases)),
            ]


if __name__ == "__main__":
    # Read the specific csv
    with open("output.csv", "r") as f:
        data = [line.split(",") for line in f]
    new = [data[0]]
    for (i, row) in enumerate(data[1::]):
        for m in generalise(row[6]):
            row[6] = m
            new.append(copy.deepcopy(row))

    with open("generalised.csv", "w") as f:
        for row in new:
            f.write(",".join(row))

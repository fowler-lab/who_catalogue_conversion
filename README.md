# Conversion of the WHO TB catalogue to GARC
Convert mutations detailed within the [WHO TB catalogue](https://www.who.int/publications/i/item/9789240028173 "WHO TB catalogue") to [GARC](https://fowlerlab.org/?p=5642 "GARC") for the purpose of compatability with [piezo](https://github.com/oxfordmmm/piezo "piezo") for use as drug resistance prediction.

## Indels
Within the WHO catalogue, the format of indel calls is somewhat ambiguous due to inconsistent positions of the indels, and that indels can be described in several ways when mixed with SNPs (which they often are within this catalogue). Furthermore, it is not uncommon for the WHO catalogue to detail long sequences of the same length for both `ref` and `alt`, which only consist of SNPs. It is also possible for the same mutation to be described through several different combinations of SNPs and indels.
To standardise this, an approach for deciphering the indels had to be developed:

Using a sliding window, the indel is placed in the position which causes the least SNPs. E.g. `actg`->`acgtg` has an insertion at `seq[2]`

If there are repeating sections which are ambiguous, the first item is chosen as the indel position. E.g. `aaa`->`aaaa` has an insertion at `seq[0]`

As long as there is a single indel within the sequence, this works. More than 1 indel in a sequence has not been observed

## Requirements
Requires [gumpy](https://github.com/oxfordmmm/gumpy "gumpy") for finding amino acid changes. Everything should be installable through pip:
```
pip install -r requirements.txt
```

## Parsing
Due to the unstable nature of the grammar used within the catalogue, it is simpler to parse the mutations and build from the ground up rather than translation. General rules are also added to ensure that any mutation input can be predicted for. These include rules such as `*= S` which means that any synonymous mutation infers susceptibility (i.e any non-change in amino acids should remain susceptable), and `*? U` which means that any nonsynonymous mutation infers an unknown outcome. This means that any mutation passed to `piezo` using this catalogue can produce a prediction.
```
python parse.py > skipped.tsv
```
This produces `output.csv` which is of a similar format to the catalogues ingested by piezo. This lacks the generalised rules which occur in some of the mutations denoted within [existing piezo catalogues](https://github.com/oxfordmmm/tuberculosis_amr_catalogues "existing piezo catalogues").
`skipped.tsv` is also produced, which details some rows of the catalogue which were skipped for various issues (see below).

## Issues
There are some rows of the WHO TB catalogue which do not cleanly produce a singluar mutation for a singular row of the catalogue. There are cases where there are several SNPs within a single row, or a mixture of SNPs and an indel (as detailed above). This becomes an issue due to the fact that these mutations should therefore be coupled together - all mutations detailed in a row should be present for the resistance prediction.
This makes it impossible to translate these rows to GARC as there is not support for chaining mutations like this, i.e something along the lines of `whiB6@a-77g & whiB6@-74_del_g` is not valid GARC syntax, and to extend GARC to handle such detail would involve serious work with all tools which utilise it.
Therefore, rows which produce more than 1 mutation must be skipped within the output.

Furthermore, there are some rows which detail mutations past the 3' end of the gene, and so are not helpful and are also ommitted.


Moreover, there exist some rows which produce mutations which are placed in more than 1 resistance category for a drug. These may be distinct rows, but as the mutations are the same with different results, they also have to be ommitted. 

E.g. gid@351_del_g is formed by `cgc -> cg @ 4407850`, but can also be formed by `gc -> c @ 4407851` - both giving deletion of the same `c` base in different rows. The former is accociated with susceptibility to streptomycin, whereas the latter is associated with resistance to streptomycin. This is due to the `interim` values used within the catalogue not concretely inferring resistance, so susceptibility must be selected instead.

A full breakdown of skipped rows can be found in `skipped.tsv`, with the general format for rows being `label   gene    pos ref alt` with some exceptions such as for mutations in > 1 category, and added mutation data for rows with multiple mutations. 

Labels:
* `Multiple mutations per row:` Is the result of skipping due to > 1 mutation in this row.
* `Cut off` Is the result of skipping due to changes detailed past the 3' end of the gene. There should only be 2 cases of this.
* `Exists in >1 resistance category:` Is the result of skipping due to the mutation belonging to > 1 resistance category for a given drug.

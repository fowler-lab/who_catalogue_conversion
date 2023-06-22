# Conversion of the WHO TB catalogue to GARC
Convert mutations detailed within the [WHO TB catalogue](https://www.who.int/publications/i/item/9789240028173 "WHO TB catalogue") to [GARC](https://fowlerlab.org/?p=5642 "GARC") for the purpose of compatability with [piezo](https://github.com/oxfordmmm/piezo "piezo") for use as drug resistance prediction.

A copy of the parsed catalogue is provided within this repo - see `WHO-UCN-GTB-PCI-2021.7.GARC.csv`. It is also stored alongside a number of tuberculosis resistance catalogues in another GitHub repo [tuberculosis-amr-catalogues](https://github.com/oxfordmmm/tuberculosis_amr_catalogues).

## Requirements
Requires [gumpy](https://github.com/oxfordmmm/gumpy "gumpy") for finding amino acid changes. Everything should be installable through pip:
```
pip install -r requirements.txt
```

## Running
First time running will take a considerable amount of time (30+ minutes) due to consistent rebuilding of `gumpy.Gene` objects for mutations. After this, values will be cached to disk using pickle. Due to the security implications of the pickle module, **DO NOT SEND/RECIVE PICKLES**. They are blacklisted in the `.gitignore` for this reason.
```
python parse.py
```

### Development
During future development, you may wish to force it to re-parse the data rather than using pickles. To do this, just add any argument to the script:
```
python parse.py <anything>
```

## Parsing notes
Due to some issues in the way the WHO catalogue is built, there are a few issues which drove design decisions. The issues found are detailed below, along with the solutions currently in place.
### Indels
Within the WHO catalogue, the format of indel calls is somewhat ambiguous due to inconsistent positions of the indels, and that indels can be described in several ways when mixed with SNPs (which they often are within this catalogue). Furthermore, it is not uncommon for the WHO catalogue to detail long sequences of the same length for both `ref` and `alt`, which only consist of SNPs. It is also possible for the same mutation to be described through several different combinations of SNPs and indels.
To standardise this, an approach for deciphering the indels had to be developed:

Using a sliding window, the indel is placed in the position which causes the least SNPs. E.g. `actg`->`acgtg` has an insertion at `seq[2]`

If there are repeating sections which are ambiguous, the first item is chosen as the indel position. E.g. `aaa`->`aaaa` has an insertion at `seq[0]`

As long as there is a single indel within the sequence, this works. More than 1 indel in a sequence has not been observed

### Multiple mutations per row
This mainly arrises due to handling of the `Indels` detailed above. This results in a split of an `indel`, and any number of `snps` from a single row of the catalogue. However, there are also cases where a single row of the catalogue has no `indel`, but rather a series of `snps` - still producing multiple mutations per row.

As these arrise from allelic variant calling, and so imply that all mutations are required for the drug predictions. To accomodate this, an extension to GARC had to be made to allow for this chaining together of mutations. This is a very simple rule which states that any number of mutations can be joined together with `&` to produce a new valid mutation.
E.g:
    ['pncA@R140S', 'pncA@453_del_gttcggta'] --> 'pncA@R140S&pncA@453_del_gttcggta'
For the sake of consistency, it is recommended that such mutations are sorted before joining, but this is also enforced within `piezo`, so is not required.

Alterations had to be made within `piezo` and `gnomon` to accomodate such mutations, but should now be working.

### Mutations past the gene end
There are a few cases where the catalogue details mutations which cross the gene end. This produces an issue as there is no other mutation at the genome level, and as one of the key assumptions made when building the catalogue was that such mutations (non-coding and non-promoter) do not confer resistance, we have opted to censor such mutations at the gene end.
This means that the parsing attempts to parse the mutations in the same manner as usual, rejecting any which fall past the gene end.

These cut off mutations are very rare, and only result in a single SNP being accepted as a valid mutation within the entire catalogue.

**Note that this will produce some STDOUT lines where this happens when parsing.** These are intended for debugging, and also show what mutations are accepted from a row in these cases.

### Mutations existing in >1 resistance category
This occurs at several points within the catalogue, and are cases in which there are rows which parse to the same mutation which have different values for resistance for the same drug.
For example, ['pncA@-5_del_g', 'pncA@317_del_t', 'pncA@386_del_atgt']  are both R and S for  PZA

The current solution is to ignore this during parsing. It is possible that this will be resolved in a future version of the catalogue. This creates a catalogue which may contain >1 row for a single mutation's effect on a drug. Alterations have been made to piezo to deal with such a catalogue, and produce a prediction based on the most significant predicted value. 
For example, `pncA@-5_del_g` has predictions of both `R` and `S` for `PZA`. Using a prioritisation of `R > U > S`, the final prediction will be `R`

## Expert rules
Separate from the catalogue, there are a set of `expert rules` which support it. The below table maps row numbers (inclusive) of `expertRules.csv` to the associated expert rule. There are further rules based on literature and similar within the report, but these are all already included within the parsed catalogue.

| Row start | Row end | Expert rule |
| --------- | ------- | ----------- |
| 2 | 3 | Borderline rpoB mutations --> RIF R |
| 4 | 103 | Non-synonymous RRDR mutations --> RIF RI | 
| 104 | 105 | Premature stop codon or indel in katG --> INH RI |
| 106 | 107 | Premature stop codon or indel in pncA --> PZA RI |
| 108 | 109 | Premature stop codon or indel in gid --> STM RI |
| 110 | 111 | Premature stop codon or indel in ethA --> ETH RI |
| 112 | 116 | According to previous WHO guidance not already in the catalogue |
| 117 | 117 | Alternate form of row above (pncA@-126_del_c) |
| 116 | 119 | From literature, not already in the catalogue |
| - | - | Any mutation within gyrA/B which confers R --> LFX R |

### Not included
Some rules cannot be included due to limitations of GARC:

| Rule | Drug | Prediction |
| ---- | ---- | ---------- |
| RIF resistance & (pncA@*? or pncA@*_indel) | PZA | R |

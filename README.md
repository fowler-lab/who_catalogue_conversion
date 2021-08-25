# Conversion of the WHO TB catalogue to GARC
Convert mutations detailed within the [WHO TB catalogue](https://www.who.int/publications/i/item/9789240028173 "WHO TB catalogue") to [GARC](https://fowlerlab.org/?p=5642 "GARC") for the purpose of compatability with [piezo](https://github.com/oxfordmmm/piezo "piezo") for use as drug resistance prediction.

## Indels
Within the WHO catalogue, the format of indel calls is somewhat ambiguous due to inconsistent positions of the indels, and that indels can be described in several ways when mixed with SNPs (which they often are within this catalogue). Furthermore, it is not uncommon for the WHO catalogue to detail long sequences of the same length for both `ref` and `alt`, which only consist of SNPs. It is also possible for the same mutation to be described through several different combinations of SNPs and indels.
To standardise this, an approach for deciphering the indels had to be developed:

Using a sliding window, the indel is placed in the position which causes the least SNPs. E.g. `actg`->`acgtg` has an insertion at `seq[2]`

If there are repeating sections which are ambiguous, the first item is chosen as the indel position. E.g. `aaa`->`aaaa` has an insertion at `seq[0]`

As long as there is a single indel within the sequence, this works. More than 1 indel in a sequence has not been observed

## Requirements
Requires [gumpy](https://github.com/oxfordmmm/gumpy "gumpy") for finding amino acid changes.

## Parsing
Due to the unstable nature of the grammar used within the catalogue, it is simpler to parse the mutations and build from the ground up rather than translation.
```
python parse.py
```
This produces `output.csv` which is of a similar format to the catalogues ingested by piezo, but can be generalised (`python generalise.py`) to denote less specific mutations e.g `Rv2752c@1579_del_gccta` could become `Rv2752@1579_indel`.
This is closer to the style mutations are denoted within [existing piezo catalogues](https://github.com/oxfordmmm/tuberculosis_amr_catalogues "existing piezo catalogues").

import argparse
import os
import re

import pandas as pd


def get_args():
    parser = argparse.ArgumentParser(description="Extract taxonomic subset.")
    parser.add_argument("Inputdir", help="directory containing the whole dataset")
    parser.add_argument("Outputdir",help="Path to output directory")
    return parser.parse_args()

def taxon_matching(rank_taxo, taxonomy_splitted):
    for taxa in taxonomy_splitted:
        search = re.findall(f".*{rank_taxo}", taxa)
        if search:
            return "[]".join(search)
    return ""

def clean_taxonomy(path,out_path):
    suffix_taxa = {"phylum": "mycota", "subphylum": "mycotina", "class": "mycetes", "subclass": "mycetidae",
                   "superorder": "anae", "order": "ales", "suborder": "ineae", "infraorder": "aria",
                   "superfamily": "acea", "family": "aceae", "subfamily": "oideae"}
    header = list(suffix_taxa.keys())
    header.append("genus")
    header.insert(0, "Species")
    df = pd.read_csv(path, sep="\t")
    species = df["#Species"]
    result = []
    taxonomy = pd.DataFrame(columns=[header])
    for index, value in enumerate(df["Taxonomy"]):
        result.clear()
        taxonomy_splitted = value.split("; ")
        result.append(species[index])
        for rank, suffix in suffix_taxa.items():
            result.append(taxon_matching(suffix, taxonomy_splitted))
        result.append(taxonomy_splitted[-1])
        taxonomy.loc[index] = result
    out_path = os.path.join(out_path, "taxonomy_clean.csv")
    taxonomy.to_csv(out_path, sep = ",", decimal='.', index=False)
    return taxonomy


if __name__ == '__main__':
    args = get_args()
    taxo_clean = clean_taxonomy(args.Inputdir,args.Outputdir)



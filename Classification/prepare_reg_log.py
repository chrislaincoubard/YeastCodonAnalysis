import pandas as pd
import os
import json
import argparse


def get_args():
    parser = argparse.ArgumentParser(prog="prepare log reg", description="data preparation for logisitic regression")
    parser.add_argument("Workingdir", help="directory containing the whole dataset")
    parser.add_argument("Outputdir",
                        help="Path to output directory (will be created if it doesn't exist)", default=os.getcwd())
    parser.add_argument("Taxopath",
                        help="Path to taxonomy file")
    parser.add_argument("-f","--family", help="families of interest", nargs = "*")
    parser.add_argument("-s", "--species",help="list of species to keep", nargs = "*")
    parser.add_argument("-gc", "--genecodetype", help="Choose which genecode type you want to include (default = W)",
                        default="W", nargs="*")
    return parser.parse_args()



def clean_file_name(file_path):
    clean_file = file_path.split("\\")[-1]
    clean_name = clean_file.split("_")[:-1]
    clean_name = " ".join(clean_name)
    return clean_name


def prepare_df(file_path, df_taxo):
    df = pd.read_csv(file_path, sep="\t")
    df = df.drop(df.columns[[2, 3, 4, 5, 6]], axis=1)
    name = clean_file_name(file_path)
    index = df_taxo.index[df_taxo['Species'] == name].tolist()
    fam = df_taxo.at[index[0], "family"]
    if fam == "Debaryomycetaceae":
        df["family"] = 0
    if fam == "Metschnikowiaceae":
        df["family"] = 1
    if fam == "Pichiaceae":
        df["family"] = 2
    if fam == "Saccharomycetaceae":
        df["family"] = 3
    return df


def extract_subset(data_taxo,families, species_to_extract =""):
    list_df = []
    for index, name_taxon in enumerate(data_taxo["family"]):
        species_name = data_taxo.at[index, "Species"]
        species_name = species_name.replace(" ","_")
        if name_taxon in families:
            for value in args.genecodetype:
                filepath = f"{args.Workingdir}\\{species_name}_{value}.csv"
                list_df.append(prepare_df(filepath, data_taxo))
                print(f"Done for this file : {species_name}_{value}.csv")
    return list_df


if __name__ == '__main__':
    args = get_args()
    try:
        os.mkdir(args.Outputdir)
    except OSError:
        pass
    df_taxo = pd.read_csv(args.Taxopath, sep="\t")
    list_df = extract_subset(df_taxo, args.family, args.species)
    df = pd.concat(list_df)
    path = os.path.join(args.Outputdir, "subset_data_reg.csv")
    print("converting to csv")
    df.to_csv(path, sep=",", decimal=".", index=False)
    print("Done")
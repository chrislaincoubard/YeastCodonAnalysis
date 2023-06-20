import pandas as pd
import os
import argparse
import numpy as np

def get_args():
    parser = argparse.ArgumentParser(prog="extract_subset", description="Extract taxonomic subset.")
    parser.add_argument("Workingdir", help="directory containing the whole dataset")
    parser.add_argument("Outputdir",
                        help="Path to output directory (will be created if it doesn't exist)", default=os.getcwd())
    parser.add_argument("Taxopath",
                        help="Path to taxonomy file")
    return parser.parse_args()



def clean_file_name(file_path):
    clean_file = file_path.split("\\")[-1]
    clean_name = clean_file.split("_")[:-1]
    return " ".join(clean_name)



def prepare_data(df):

    clean_name = clean_file_name(file)
    df2 = df.drop(df.columns[[0, 1, 2, 3, 4, 5, 6]], axis=1)
    colnames = list(df2)

    df3 = pd.DataFrame(columns=(["Species"] + colnames))

    index = list(df_taxo["Species"]).index(clean_name)

    liste = [list(np.quantile(df[x], [0.5])) for x in list(df2)]
    liste.insert(0, [clean_file_name(file)])
    print(liste)
    for index, value in enumerate(list(df3)):
        df3[value] = liste[index]
    return df3



if __name__ == '__main__':
    args = get_args()
    try:
        os.mkdir(args.Outputdir)
    except OSError as error:
        pass
    df_taxo = pd.read_csv(args.Taxopath, sep="\t")
    median_A, median_W, median_C = [], [], []
    for file in os.listdir(args.Workingdir):
        if file.endswith(".csv"):
            f = os.path.join(args.Workingdir, file)
            df = pd.read_csv(f, sep="\t")
            if file.endswith("W.csv"):
                median_W.append(prepare_data(df))
            if file.endswith("C.csv"):
                median_C.append(prepare_data(df))
            if file.endswith("A.csv"):
                median_A.append(prepare_data(df))
    df_to_save_A = pd.concat(median_A)
    df_to_save_W = pd.concat(median_W)
    df_to_save_C = pd.concat(median_C)
    df_to_save_A.to_csv(os.path.join(args.Outputdir,"Median_metrics_A.csv"), index = False)
    df_to_save_C.to_csv(os.path.join(args.Outputdir,"Median_metrics_C.csv"), index = False)
    df_to_save_W.to_csv(os.path.join(args.Outputdir,"Median_metrics_W.csv"), index = False)

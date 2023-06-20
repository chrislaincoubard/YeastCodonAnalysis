import argparse
import os
import random

import numpy as np
import pandas as pd
import sklearn.metrics as mt
from imblearn.over_sampling import ADASYN
from sklearn.linear_model import SGDClassifier
from sklearn.model_selection import train_test_split
from sklearn.utils import resample


def get_args():
    parser = argparse.ArgumentParser(description="Perform classification on the data set with multiple conditions")
    parser.add_argument("Inputdir", help="directory containing the whole dataset")
    parser.add_argument("Outputdir",
                        help="Path to output directory ")
    return parser.parse_args()

def test_model(model, x_train, y_train, x_test, y_test, df, name):
    y_test2 = np.array(y_test)
    metrics = []
    fit = model.fit(x_train, y_train)
    y_pred = model.predict(x_test)
    cnf = mt.confusion_matrix(y_test, y_pred)
    print("Real CNF \n")
    print(cnf)
    cnf_dataframe = pd.DataFrame(cnf, columns=["Debaryomycetaceae","Metschnikowiaceae","Pichiaceae","Saccharomycetaceae"])
    cnf_dataframe.to_csv(f"{args.Outputdir}\\cnf_{name}_{j}.csv", sep = ",", index = False)
    print("\n")
    print("CNF Dataframe \n")
    print(cnf_dataframe)
    coef_dataframe = pd.concat([pd.DataFrame(x_train.columns), pd.DataFrame(np.transpose(fit.coef_))], axis = 1)
    coef_dataframe.to_csv(f"{args.Outputdir}\\coef_{name}_{j}.csv", sep=",", index = False)
    metrics.append(mt.balanced_accuracy_score(y_test, y_pred))
    metrics.append(mt.matthews_corrcoef(y_test, y_pred))
    metrics.extend(mt.recall_score(y_test, y_pred, average=None))
    metrics = [round(float(x), 3) for x in metrics]
    metrics.insert(0, name)
    df.loc[len(df)] = metrics
    y_test1 = np.array(y_test)
    if (y_test2 == y_test1).all() :
        print("yipee")
    else :
        print("nope")
    return 0


def subset_species(data, list_species, p=0.75):
    df_to_return = pd.DataFrame(columns=data.columns.values.tolist())
    for species in list_species:
        rows, columns = np.where(data == species)
        index_to_keep = random.sample(list(rows), int(round(p * len(rows), 0)))
        df = data.iloc[index_to_keep]
        df_to_return = pd.concat([df_to_return, df])
    return df_to_return

def upsample(data, label):
    data['family'] = label.values
    list_class = []
    Upsample = max(data["family"].value_counts())
    for value in set(label.values) :
        fam = data[data["family"] == value]
        list_class.append(resample(fam, n_samples=Upsample, replace=True))
    X = pd.concat(list_class)
    Y = X["family"]
    X = X.drop(["family"], axis=1)
    return X, Y

def downsample(data, label):
    data['family'] = label.values
    list_class = []
    Downsample = min(data["family"].value_counts())
    for value in set(label.values):
        fam = data[data["family"] == value]
        list_class.append(resample(fam, n_samples=Downsample, replace=True))
    X = pd.concat(list_class)
    Y = X["family"]
    X = X.drop(["family"], axis=1)
    return X, Y


def prepare_data_train(data):
    data = data.dropna()
    print(list(data))
    if "#CDSid" in list(data):
        Y_label = data[["family"]]
        X_values = data.drop(["#CDSid", "Species", "family"], axis=1)
        data_partition = train_test_split(X_values, Y_label.values.ravel(),stratify=Y_label.values.ravel(), test_size=0.2)
        train = [data_partition[0], data_partition[2]]
    else :
        Y_label = data[["family"]]
        X_values = data.drop(["family"], axis=1)
        data_partition = train_test_split(X_values, Y_label.values.ravel(),stratify=Y_label.values.ravel(), test_size=0.2)
        train = [data_partition[0], data_partition[2]]
    return train


def prepare_data_test(data):
    data = data.dropna()
    Y_label = data[["family"]]
    if "#CDSid" in list(data):
        X_values = data.drop(["#CDSid", "Species", "family"], axis = 1)
        data_partition = train_test_split(X_values, Y_label.values.ravel(),stratify=Y_label.values.ravel(), test_size=0.2)
        test = [data_partition[1], data_partition[3]]
    else :
        X_values = data.drop(["family"], axis=1)
        data_partition = train_test_split(X_values, Y_label.values.ravel(),stratify=Y_label.values.ravel(), test_size=0.2)
        test = [data_partition[1], data_partition[3]]
    return test

if __name__ == '__main__':
    args = get_args()
    dir = args.Inputdir
    for file in os.listdir(dir):
        df_tosave = pd.DataFrame(columns=["type", "MCC", "Bacc", "r0", "r1", "r2", "r3"])
        if file.endswith(".csv"):
            print(f"Processing of : {file}")
            path  = os.path.join(dir, file)
            data = pd.read_csv(path, sep=',')
            data = data.dropna()

            df_noMC = data.drop(
                ["GCcds", "CpG", "SCUO", "CBI", "FPC", "FAC", "BOC", "BIC", "Gravy"], axis=1)
            Y = df_noMC["family"]
            X = df_noMC.drop(["#CDSid", "Species", "family"], axis=1)
            species = set(data["Species"].tolist())


            x_train_noMC, y_train_noMC = prepare_data_train(df_noMC)

            X_upsampling_noMC, Y_upsampling_noMC = upsample(x_train_noMC.copy(), pd.Series(y_train_noMC))
            X_downsampling_noMC, Y_downsampling_noMC = downsample(x_train_noMC.copy(), pd.Series(y_train_noMC))
            x_ada, y_ada = ADASYN().fit_resample(X, Y)
            data_partition = train_test_split(x_ada, y_ada.values.ravel(), stratify=y_ada.values.ravel(),test_size=0.2)
            x_train_ada, y_train_ada = (data_partition[0], data_partition[2])
            print(list(x_train_ada))


            elnet_SGD = SGDClassifier(max_iter=1000, tol=1e-3, loss="log_loss", penalty='elasticnet')
            lasso_SGD = SGDClassifier(max_iter=1000, tol=1e-3, loss="log_loss", penalty='l1')
            ridge_SGD = SGDClassifier(max_iter=1000, tol=1e-3, loss="log_loss", penalty='l2')
            classic_SGD = SGDClassifier(max_iter=1000, tol=1e-3, loss="log_loss", penalty=None)
            print("Models creaton : Done")
            if file.endswith(".csv"):
                for j in [1, 2, 3, 4, 5]:

                    print(f"Processing run {j} out of 5.")
                    x_test, y_test = prepare_data_test(df_noMC)
                    x_test2 = np.array(x_test)

                    test_model(classic_SGD, x_train_noMC, y_train_noMC, x_test,y_test, df_tosave, f"NoPenalty")
                    x_test1 = np.array(x_test)
                    if (x_test2 == x_test1).all():
                        print("im here")
                    else :
                        print("stop me")
                    test_model(elnet_SGD, x_train_noMC, y_train_noMC, x_test, y_test, df_tosave, f"ElasticNet")
                    test_model(lasso_SGD, x_train_noMC, y_train_noMC, x_test, y_test, df_tosave, f"Lasso")
                    test_model(ridge_SGD, x_train_noMC, y_train_noMC, x_test, y_test, df_tosave, f"Ridge")
                    test_model(lasso_SGD, X_upsampling_noMC, Y_upsampling_noMC, x_test, y_test, df_tosave, f"Oversamplng")
                    test_model(lasso_SGD, X_downsampling_noMC, Y_downsampling_noMC, x_test, y_test, df_tosave, f"Undersampling")
                    test_model(lasso_SGD, x_ada, y_ada, x_test, y_test, df_tosave, f"ADASYN")
                print(f"Done for {file}")
            else :
                for j in [1,2,3,4,5]:
                    print(f"Processing run {j} out of 5.")
                    x_test, y_test = prepare_data_test(df_noMC)

                    test_model(lasso_SGD, x_train_noMC, y_train_noMC, x_test, y_test, df_tosave, f"Lasso")
                print(f"Done for {file}")

        df_tosave.to_csv(f"{args.Outputdir}\\{file}_final_bench.csv", sep=",")
    print("Execution done")

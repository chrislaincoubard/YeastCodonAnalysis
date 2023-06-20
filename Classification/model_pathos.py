import os

import pandas as pd
import sklearn.metrics as mt
from imblearn.over_sampling import ADASYN
from sklearn.linear_model import SGDClassifier
from sklearn.model_selection import train_test_split
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="Extract taxonomic subset.")
    parser.add_argument("File", help="File with dataset")
    parser.add_argument("Outputdir",
                        help="Path to output directory")
    return parser.parse_args()

if __name__ == '__main__':
    args = get_args()
    data_ft_and_pathos = pd.read_csv(args.File, sep='\t')
    Y = data_ft_and_pathos["Species"]
    X = data_ft_and_pathos.drop(data_ft_and_pathos.columns[[26]], axis=1)
    X_ada, Y_ada = ADASYN().fit_resample(X, Y)
    X_train_ft_pathos, X_test_ft_pathos, Y_train_ft_pathos, Y_test_ft_pathos = train_test_split(X, Y.ravel(), test_size=0.2)


    def test_model(model, x_train, y_train, x_test, y_test, df, name):
        metrics = []
        fit = model.fit(x_train, y_train)
        y_pred = model.predict(x_test)
        cnf = mt.confusion_matrix(y_test, y_pred)
        cnf_dataframe = pd.DataFrame(cnf, columns=["Pathos", 'FT (non patho)'])
        print(cnf_dataframe)
        cnf_dataframe.to_csv(f"cnf_{name}_{i}.csv", sep=",", index=False)
        metrics.append(mt.balanced_accuracy_score(y_test, y_pred))
        metrics.append(mt.matthews_corrcoef(y_test, y_pred))
        metrics.extend(mt.recall_score(y_test, y_pred, average=None))
        metrics = [round(float(x), 3) for x in metrics]
        metrics.insert(0, name)
        df.loc[len(df)] = metrics
        return 0


    df_tosave = pd.DataFrame(columns=["type", 'MCC', 'Bacc', 'r0', 'r1'])
    data_partition = train_test_split(X_ada, Y_ada.values.ravel(), stratify=Y_ada.values.ravel(), test_size=0.2)
    X_train_ada, Y_train_ada = (data_partition[0], data_partition[2])

    elnet_SGD = SGDClassifier(max_iter=1000, tol=1e-3, loss="log_loss", penalty='elasticnet')
    lasso_SGD = SGDClassifier(max_iter=1000, tol=1e-3, loss="log_loss", penalty='l1')
    ridge_SGD = SGDClassifier(max_iter=1000, tol=1e-3, loss="log_loss", penalty='l2')
    classic_SGD = SGDClassifier(max_iter=1000, tol=1e-3, loss="log_loss", penalty=None)

    for i in range(1, 11):
        print(f"Processing run {i} out of 10.")
        partition = train_test_split(X_ada, Y_ada.values.ravel(), test_size=0.2)
        X_test_ada, Y_test_ada = (partition[1], partition[3])
        test_model(lasso_SGD, X_train_ada, Y_train_ada, X_test_ada, Y_test_ada, df_tosave, f"Adasyn")
        print(f"Done for lasso --> {i}")

    print("Saving to csv : ....")
    path = os.path.join(args.Outputdir, "ft_pathos_only_result.csv")
    df_tosave.to_csv(path, sep=",", index=False)
    print("Done")

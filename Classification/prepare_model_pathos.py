import pandas as pd
import os

os.chdir("C:/Users/chris/OneDrive/Documents/Master_bio_info_Biblio/Stage_M2/code/data/pathos_subset")

X_pathos = pd.read_csv("les_pathos.csv",sep=",", encoding="utf-8")

X = pd.read_csv("subset_pathos_all.csv", sep="\t",encoding="utf-7" )
X = X.dropna(axis = 0)
# X_pathos = X_pathos.dropna(axis=0)
X_FT = pd.read_csv("only_ft.csv", sep = "\t")
X_FT = X_FT.dropna()
Species = X["Species"]
# CDS_list_patho = X_pathos["#CDSid"]
CDS_list_FT = X_FT["#CDSid"]
CDS_ALL = X["#CDSid"]
Species_class = []
CDS_class = []
X = X.drop(X.columns[[0,1,2,3,4,5,6]], axis=1)



# X_pathos = X_pathos.drop(X_pathos.columns[[0,1,2,3,4,5,6,7]], axis = 1)
X_FT = X_FT.drop(X_FT.columns[[0,1,2,3,4,5,6,7]], axis = 1)



print(len(Species_class))
X["Species"] = Species_class
# X.to_csv("only_pathos_with_class.csv", sep=',', index=False)





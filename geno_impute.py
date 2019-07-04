import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt

from sklearn.impute import SimpleImputer


def geno_imputation(geno_data_file, pheno_data_file, merge_imputed_file):
	# GET DATA FROM FILES
	A = pd.read_csv(geno_data_file, delim_whitespace=True)
	del A["FID"] # for 3K data FID is the same as IID. For other datasets, we may concatenate FID and IID to make a unique ID
	del A["PAT"] # These are not variable
	del A["MAT"] #
	del A["SEX"] #
	del A["PHENOTYPE"] # This is all-missing. The phenotype will be provided from another file

	P = pd.read_csv(pheno_data_file)
	# Adjust the phenotype data so that IDs can match. SNP-Seek IDs have a space, should be converted to "_"
	newID = np.array([ x.replace(" ", "_") for x in P["IRIS ID"] ])
	P["IRIS ID"] = newID
	#---------------------
	# Some samples may have been removed from genotype data.

	"""
	pheno = P.iloc[:,5]
	plt.hist(pheno, bins=30, color="#339933")
	plt.show()
	"""
	#Merging genotype and phenotype
	M = pd.merge(P, A, left_on = "IRIS ID", right_on ="IID" )
	[*M][:10]
	Y = M.iloc[:,5]
	X = M.iloc[:, 7:]

	#Imputation
	# use mean imputation
	imp = SimpleImputer(missing_values=np.nan, strategy='mean', copy=True)
	A_geno = np.array(A.iloc[:,1:])
	imp.fit(A_geno) # this comutes the means of columns, to be used instead of missing values
	A_geno_imputed = imp.transform(A_geno)

	# Combine the imputed data with the IDs
	A_imputed = pd.DataFrame(A_geno_imputed, columns=[*A][1:] )
	A_imputed["IID"] = A["IID"]
	A_imputed2 = A_imputed.iloc[:, [A_imputed.shape[1]-1, *list(range(A_imputed.shape[1]-1)) ] ]
	imputed_file = geno_data_file + "_imputed.csv"
	A_imputed2.to_csv(imputed_file)

	M = pd.merge(P, A_imputed2, left_on = "IRIS ID", right_on ="IID" )
	[*M][:10]
	M.to_csv(merge_imputed_file, index=False)

	Y = M.iloc[:,5]
	X = M.iloc[:, 7:]

	return X, Y
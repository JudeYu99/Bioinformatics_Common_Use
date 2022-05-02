#! /usr/bin

# @Author: Yu Zhu
# @Email: 1830416012@stu.suda.edu.cn
# @Address: Department of Bioinformatics, Medical College, Soochow University

## NOTE: Python 3.7.0 code

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn import preprocessing
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

###########################################################################
#=======================---------------------------=======================#
#                                   PCA                                   #
#=======================---------------------------=======================#
###########################################################################

def myPCA(filename):
    
    data = pd.read_csv(filename)

    #########################
    #
    # Perform PCA on the data
    #
    #########################
    scaled_data = StandardScaler().fit_transform(data)
    pca = PCA()
    pca.fit(scaled_data) 
    pca_data = pca.transform(scaled_data) 

    #########################
    #
    # Draw a scree plot and a PCA plot
    #
    #########################
    per_var = np.round(pca.explained_variance_ratio_* 100, decimals = 1)
    labels = ['PC' + str(x) for x in range(1, len(per_var) + 1)]

    plt.bar(x = range(1, len(per_var) + 1), height = per_var)
    plt.ylabel('Percentage of Explained Variance')
    plt.xlabel('Principal Component')
    plt.title('Scree Plot')
    plt.savefig(filename + "_Scree.png")
    plt.show()

    pca_df = pd.DataFrame(pca_data, columns = labels)

    plt.scatter(pca_df.PC1, pca_df.PC2)
    plt.title('PCA Graph')
    plt.xlabel('PC1 - {0}%'.format(per_var[0]))
    plt.ylabel('PC2 - {0}%'.format(per_var[1]))

    for item in pca_df.index:
        plt.annotate(item, (pca_df.PC1.loc[item], pca_df.PC2.loc[item]))

    plt.savefig(filename + "_PCA.png")
    plt.show()



###########################################################################
#=======================---------------------------=======================#
#                                Analysis                                 #
#=======================---------------------------=======================#
###########################################################################

def get_PCs(per_var):
    PC_sums = 0
    PC_counts = 0
    for i in per_var:
        PC_sums += i
        PC_counts += 1
        if PC_sums >= 80:
            break
    return PC_counts

def ZSCORE(Company_num, scaled_data, pca, per_var, PCs):
    Zscore = 0
    for z in range(PCs):
        for j in range(PCs):
            No_y = 0
            for i in range(len(scaled_data[Company_num])):
                No_y +=  scaled_data[Company_num][i] * pca.components_[j][i]

            Zscore += (per_var/100)[z] * No_y
    #print(Zscore)
    return Zscore

def PCA_all_score(filename):
    data = pd.read_csv(filename)
    scaled_data = StandardScaler().fit_transform(data)
 
    pca = PCA()
    pca.fit(scaled_data)
    pca_data = pca.transform(scaled_data) 
    
    per_var = np.round(pca.explained_variance_ratio_* 100, decimals = 1)
    
    PCs = get_PCs(per_var)
    
    result = []
    for k in range(scaled_data.shape[0]):
        result.append(ZSCORE(k, scaled_data, pca, per_var, PCs))
        
    return result

  

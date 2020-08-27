from mpl_toolkits.mplot3d import Axes3D
import csv
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from collections import Counter


def del_list_inplace(l, id_to_del):
    for i in sorted(id_to_del, reverse=True):
        del l[i]
    return l

def shufflingBackwords(clustData, seqIDs):

    fn = 'data/sars_cov_2_fixed.fasta'
    seqIDs_initial_order = []
    
    # load seqIDs with the initial order
    with open(fn, 'r') as fh:
        for line in fh:
            if line[0] == '>':
                seqIDs_initial_order.append(findSequenceId(line))
    
    reorderingVector = []
    for id in seqIDs:
        reorderingVector.append(seqIDs_initial_order.index(id))
    
    seqIDs = [x for _,x in sorted(zip(reorderingVector,seqIDs))]
    clustData = [x for _,x in sorted(zip(reorderingVector,clustData))]

    return clustData, seqIDs

def findSequenceId(myline):
    
    index = myline.find('|')
    seqID = myline[1:index-1]

    return seqID


def loadClusteringData():

    # Path from data File and path for clustering
    dataFile = 'Output/sars_cov_2_fixed/output.csv'
    seqIndicesFile = 'Output/sars_cov_2_fixed/seqIndices.csv'
    seqIDsFile = 'Input/sars_cov_2_fixed_sequencesIDs.csv'
    clusteringDataPath = 'ClusteringData/clustData.csv'
    clusteringDataPathNoHeaers = 'ClusteringData/clustData_no_headres.csv'
    
    with open(dataFile, newline='') as csvfile:
        cpamreader = csv.reader(csvfile)
        inputData = list(cpamreader)
    
    with open(seqIndicesFile, newline='') as csvfile:
        cpamreader = csv.reader(csvfile)
        seqIndices = list(cpamreader)
    
    with open(seqIDsFile, newline='') as csvfile:
        cpamreader = csv.reader(csvfile)
        seqIDslist = list(cpamreader)
        seqIDs = seqIDslist[0]
    
    to_be_deleted = []
    for i in range(len(inputData)):
        if inputData[i][2] == str(0) or  inputData[i][2] == str(281):
            to_be_deleted.append(i)

    inputData = del_list_inplace(inputData, to_be_deleted)
    seqIndices =  del_list_inplace(seqIndices, to_be_deleted)
    kmersList = [inputData[i][0] for i in range(len(inputData))]
    
    # initialization
    clustData = [[0 for _ in range(len(kmersList))] for _ in range(len(seqIDs))]
    
    for i in range(len(seqIndices)):
        for j in range(len(seqIDs)):
            if str(j) in seqIndices[i]:
                clustData[j][i] = 1

    clustData, seqIDs = shufflingBackwords(clustData, seqIDs)      
    dfObj = pd.DataFrame(clustData, columns = kmersList, index= seqIDs)
    dfObj.to_csv(clusteringDataPath, index=seqIDs, columns=kmersList)

    with open(clusteringDataPathNoHeaers, 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        for i in range(len(clustData)):
            # Add contents of list as last row in the csv file
            csv_writer.writerow(clustData[i])
    
    return clustData, dfObj

def clusteringKmeans(clustData, numClusters):

    # PCA algorithm
    pca = PCA(n_components=3)
    principalComponents = pca.fit_transform(clustData)
    
    # Kmeans algorithm
    kmeans = KMeans(n_clusters=numClusters).fit(clustData)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('Principal Component 1', fontsize = 15)
    ax.set_ylabel('Principal Component 2', fontsize = 15)
    ax.set_zlabel('Principal Component 3', fontsize = 15)
    ax.set_title('Clustering SARS-Cov-2 genome sequences', fontsize = 20)
    
    targets = kmeans.labels_

    colors = np.array([x for x in 'bgrcmykbgrcmykbgrcmykbgrcmyk'])
    colors = np.hstack([colors] * 20)
    center_colors = colors[:numClusters]

    for i in range(len(principalComponents)):
        ax.scatter(principalComponents[i][0],principalComponents[i][1],principalComponents[i][2], c = center_colors[targets[i]], label = str(targets[i])) 

    ax.legend()
    ax.grid()
    plt.show()

    return principalComponents, targets

def classificationWithMetadata(principalComponents, clustData, targets):
    
    metaDatafile = 'ClusteringDataMixed\kmer_analysis_and_meta_data__fixed_merged.csv'
    metadata = pd.read_csv(metaDatafile)
    date = metadata['Collection_Date']
    loc = metadata['Geo_Location']
    source = metadata['Isolation_Source']
    seqIDs = metadata['Accession']

    c = list(zip(seqIDs, date, loc, source, targets))
    c = [x for _,x in sorted(zip(targets,c))]

    cols=['Accession', 'Collection_Date', 'Geo_Location', 'Isolation_Source', 'Cluster']
    dfObj = pd.DataFrame(c, columns = cols)
    dfObj.to_csv('metanalysis.csv',index = False,  columns = cols)    

    return


if __name__ == "__main__":
    
    # have to do sth
    clustData, dfObj = loadClusteringData()
    principalComponents, targets = clusteringKmeans(clustData,5)
    classificationWithMetadata(principalComponents, clustData, targets)
    




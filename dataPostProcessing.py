import csv
import pandas as pd
# import numpy as np
# from collections import Counter

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
                seqIDs_initial_order.append(line)
    
    reorderingVector = []
    for id in seqIDs:
        reorderingVector.append(seqIDs_initial_order.index(id[0]))
    
    seqIDs = [x for _,x in sorted(zip(reorderingVector,seqIDs))]
    clustData = [x for _,x in sorted(zip(reorderingVector,clustData))]

    return clustData, seqIDs


def loadClusteringData():

    # Path from data File and path for clustering
    dataFile = 'Output/sars_cov_2_fixed/output.csv'
    seqIndicesFile = 'Output/sars_cov_2_fixed/seqIndices.csv'
    timesPerSeqFile = 'Output/sars_cov_2_fixed/timesPerSeq.csv'
    seqIDsFile = 'Input/sars_cov_2_fixed_sequencesIDs.csv'
    clusteringDataPath = 'ClusteringData/clustData.csv'
    clusteringDataPathNoHeaers = 'ClusteringData/clustData_no_headres.csv'
    
    with open(dataFile) as csvfile:
        cpamreader = csv.reader(csvfile)
        inputData = list(cpamreader)
    
    with open(seqIndicesFile) as csvfile:
        cpamreader = csv.reader(csvfile)
        seqIndices = list(cpamreader)
    
    with open(timesPerSeqFile) as csvfile:
        cpamreader = csv.reader(csvfile)
        timesPerSeq = list(cpamreader)
    
    seqs = pd.read_csv(seqIDsFile)
    seqIDs = []
     
    for index, rows in seqs.iterrows(): 
        my_list =[rows.Row, rows.ID] 
        seqIDs.append(my_list) 

    to_be_deleted = []
    for i in range(len(inputData)):
        if inputData[i][2] == str(0) or  float(inputData[i][3]) <= 0:
            to_be_deleted.append(i)
    
    inputData = del_list_inplace(inputData, to_be_deleted)
    seqIndices =  del_list_inplace(seqIndices, to_be_deleted)
    timesPerSeq = del_list_inplace(timesPerSeq, to_be_deleted)

    kmersList = [inputData[i][0] for i in range(len(inputData))]
    
    # initialization
    clustData = [[0 for _ in range(len(kmersList))] for _ in range(len(seqIDs))]
    
    for i in range(len(seqIndices)):
        for j in range(len(seqIDs)):
            if str(j) in seqIndices[i]:
                clustData[j][i] = timesPerSeq[i][seqIndices[i].index(str(j))]
    
    clustData, seqIDs = shufflingBackwords(clustData, seqIDs)    
    this_seqIDs = [seqIDs[i][1] for i in range(len(seqIDs))]

    dfObj = pd.DataFrame(clustData, columns = kmersList, index= this_seqIDs)
    dfObj.to_csv(clusteringDataPath, index=this_seqIDs, columns=kmersList)
    
    # SOS: in order to read: data = pd.read_csv(path, index_col = 0)

    return clustData, dfObj

if __name__ == "__main__":
    
    # have to do sth
    clustData, dfObj = loadClusteringData()

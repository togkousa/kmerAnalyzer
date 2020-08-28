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
    
    with open(dataFile) as csvfile:
        cpamreader = csv.reader(csvfile)
        inputData = list(cpamreader)
    
    with open(seqIndicesFile) as csvfile:
        cpamreader = csv.reader(csvfile)
        seqIndices = list(cpamreader)
    
    with open(seqIDsFile) as csvfile:
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

    with open(clusteringDataPathNoHeaers, 'w') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        for i in range(len(clustData)):
            # Add contents of list as last row in the csv file
            csv_writer.writerow(clustData[i])
    
    return clustData, dfObj

if __name__ == "__main__":
    
    # have to do sth
    clustData, dfObj = loadClusteringData()

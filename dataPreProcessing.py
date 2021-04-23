# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 19:51:35 2019

@author: User
"""
from random import shuffle
import sys
import os
import csv
import pandas as pd


def read_fastq_file(fn):
    
    data = []
    seqIDs = []
    id_counter = 1

    with open(fn, 'r') as fh:
        
        lines = []
        
        for line in fh:
            if line[0] != '>':
                lines.append(line.rstrip())

            else:
                data.append(''.join(lines))
                seqIDs.append([line, 'ID-' + str(id_counter)])
                id_counter += 1
                lines = []

    
        data.append(''.join(lines))
    
    del data[0]
            
    return data, seqIDs


def save_results_to_file(data, name, outputFolder, seqIDs, seqName):
    
    currPath = os.getcwd()
    os.chdir(outputFolder)
    
    if os.path.exists(name):
        os.remove(name)
    
    if os.path.exists(seqName):
        os.remove(seqName)
    
    with open(name, 'w') as f:
        for item in data:
            f.write("%s\n" % item)
    
    rows = [seqIDs[i][0] for i in range(len(seqIDs))]
    IDs = [seqIDs[i][1] for i in range(len(seqIDs))]

    d = {'Row': rows, 'ID':IDs}

    df = pd.DataFrame.from_dict(d)
    df.to_csv(seqName, index=False)

    os.chdir(currPath)
    
    return


def file_len(fname):
    
    line_len = 0
    with open(fname) as f:
        for i, l in enumerate(f):
            if i == 0:
                line_len = len(l)
            pass
    
    return i+1, line_len


def createDictionaryUnshuffled(outputFolder, seqIDs, seqName):
    
    currPath = os.getcwd()
    os.chdir(outputFolder)
    
    if os.path.exists(seqName):
        os.remove(seqName)
    
    rows = [seqIDs[i][0] for i in range(len(seqIDs))]
    IDs = [seqIDs[i][1] for i in range(len(seqIDs))]

    d = {'Row': rows, 'ID':IDs}

    df = pd.DataFrame.from_dict(d)
    df.to_csv(seqName, index=False)

    os.chdir(currPath)
    
    return



if __name__ == "__main__":
    
    inputFolder = 'data'
    outputFolder = 'Input'

    files = os.listdir(inputFolder)
    

    for file in files:
        
        name = file[0:-6] + '.txt'
        seqName = file[0:-6] + '_sequencesIDs.csv'
        seqName_unshuffled = file[0:-6] + '_sequencesIDs_unshuffled.csv'

        data, seqIds = read_fastq_file(inputFolder + '/' + file)
        createDictionaryUnshuffled(outputFolder, seqIds, seqName_unshuffled)
        createDictionaryUnshuffled('Output', seqIds, 'heades_to_IDS.csv')

        #c = list(zip(seqIds, data))
        #shuffle(c)
        #seqIds, data = zip(*c)
        
        save_results_to_file(data, name, outputFolder, seqIds, seqName)

        del data
    
    
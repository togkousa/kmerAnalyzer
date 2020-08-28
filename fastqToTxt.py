# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 19:51:35 2019

@author: User
"""
from random import shuffle
import sys
import os
import csv

def read_fastq_file(fn):
    
    data = []
    seqIDs = []

    with open(fn, 'r') as fh:
        
        lines = []
        
        for line in fh:
            if line[0] != '>':
                lines.append(line.rstrip())
            else:
                data.append(''.join(lines))
                seqIDs.append(findSequenceId(line))
                lines = []
    
        data.append(''.join(lines))
    
    del data[0]
            
    return data, seqIDs

def findSequenceId(myline):
    
    index = myline.find('|')
    seqID = myline[1:index-1]

    return seqID

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

    with open(seqName, 'w') as write_obj:
        csv_writer = csv.writer(write_obj)
        csv_writer.writerow(seqIDs)
    
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

if __name__ == "__main__":
    
    inputFolder = 'data'
    outputFolder = 'Input'

    files = os.listdir(inputFolder)
    
    for file in files:
        
        name = file[0:-6] + '.txt'
        seqName = file[0:-6] + '_sequencesIDs.csv'
        data, seqIds = read_fastq_file(inputFolder + '/' + file)
        c = list(zip(seqIds, data))
        shuffle(c)
        seqIDs, data = zip(*c)
        save_results_to_file(data, name, outputFolder, seqIDs, seqName)
        del data
    
    
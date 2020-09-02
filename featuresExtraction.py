# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 11:15:44 2019

@author: User
"""

import TreeClass
import saveToFile
import os
import fastqToTxt

# Some useful global variables
mychars = ['A', 'C', 'G', 'T']
iterations = 10000000
training_perc = 0.4
kmin = 4
kmax = 7
#topN = 100000

# Adding all children in a node
def add_all_nodes(current, depth):
    
    chars = ['A', 'C', 'G', 'T']
    for char in chars:
        current.add_child(TreeClass.Node(char, current, depth+1))
    
    return current
    
# Initialization of the tree - depth = 4
def initialize_tree():
    
    tree = TreeClass.Tree()
    root = tree.root
    
    root = add_all_nodes(root , 0)
    for child1 in root.children:
        child1 = add_all_nodes(child1, 1)
        for child2 in child1.children:
            child2 = add_all_nodes(child2, 2)
            for child3 in child2.children:
                child3 = add_all_nodes(child3, 3)
    
    return root, tree

# Transforming a sequence into numbers in order to define a path
# For example the path of ATGC is 0321
def kmer2path(kmer):
    
    path = []
    
    for i in range(len(kmer)):
        
        if kmer[i] in mychars:
            path.append(mychars.index(kmer[i]))
        else:
            return -1
    
    return path

# routine_1 is only for the first scan of file
# This means that routine1 is used only for k = 4
def routine_1(fn, k, tree, numOfLines):
   
    # Open file
    with open(fn, 'r') as fh:
        
        # Initialization
        kmers_examined = 0
        count = 0
        
        for myline in fh:
            count += 1
            if count == iterations:
                break
            
            print(str(count) + ", " + str(k) )
            path = []
            
            for j in range(len(myline) - k):
                
                kmers_examined += 1
                this_kmer = myline[j:j+k]
                if j == 0 or path == -1:
                    path = kmer2path(this_kmer)
                    if path == -1:
                        continue
                else:
                    path.remove(path[0])
                    if this_kmer[-1] in mychars:
                        path.append(mychars.index(this_kmer[-1]))
                    else:
                        path = -1
                        continue
                
                if count < int(numOfLines * training_perc):
                    tree.find_in_tree(path, False , kmers_examined, k, False, sequenceIndex=count-1)
                else:
                    tree.find_in_tree(path, True, kmers_examined, k, False,  sequenceIndex=count-1)
            
    return tree
 
    

def routine_2(fn, k, tree, numOfLines):
    
    kmers_examined = 0
    
    with open(fn, 'r') as fh:
        
        count = 0
        
        for myline in fh:   
            count += 1
            if count == iterations:
                break
            
            print(str(count) + ", " + str(k) )
            path = []
            
            for j in range(len(myline) - k):
                
                kmers_examined += 1
                this_kmer = myline[j:j+k]
                
                if j == 0 or path == -1:
                    path = kmer2path(this_kmer)
                    if path == -1:
                        continue
                
                else:
                    path.remove(path[0])
                    if this_kmer[-1] in mychars:
                        path.append(mychars.index(this_kmer[-1]))
                    else:
                        path = -1
                        continue
                    
                if count < int(numOfLines * training_perc):
                    tree.find_in_tree(path, False , kmers_examined, k, True, sequenceIndex=count-1)
                else:
                    tree.find_in_tree(path, True, kmers_examined, k, True, sequenceIndex=count-1)
         
        TreeClass.check_tree(root, kmers_examined, k)
         
    return tree


# Present the tree using a list
def listTree(tree, curr_node, sequence, treelist, zipList, k):
    
    if  not curr_node.children:
        
        # Forming results
        this_seq = sequence + curr_node.char
        treelist.append([this_seq, curr_node.depth, curr_node.count ,curr_node.evaluation, curr_node.sequenceIndices, curr_node.timesPerSeq])
        
        # Calculating entropy
        if curr_node.depth == k:
            
            zipList.append([this_seq, curr_node.count])
        
        return treelist, zipList
    
    else:
        this_seq = sequence + curr_node.char
        
        for child in curr_node.children:
            
            newnode = tree.move_to_child(curr_node, mychars.index(child.char))
            treelist, zipList = listTree(tree, newnode, this_seq, treelist,zipList , k)
    
    return treelist, zipList


def listSortBasedOnEvaluation(sub_li): 
    
    return sorted(sub_li, key = lambda x: x[3], reverse = True)

if __name__ == "__main__":
    
    # preprocessing my input data
    exec(open("fastqToTxt.py").read())
    
    # Cd
    folder = 'Input'
    files = os.listdir(folder)
    
    # Start
    for filename in files:
        if filename.endswith('.txt'):

            file = folder + '/' + filename
            numOfLines, lenOfLine = fastqToTxt.file_len(file)
            filename = filename[0:-4]
            saveToFile.checkIfOutputFileExists(filename)
            
            # Range of k values
            kvals = list(range(kmin,kmax+1))
            
            # Initialization
            root, tree = initialize_tree()
            
            # Main loop 
            for k in kvals:
                
                # Only the first time run the first routine
                if k == 4:
                    tree = routine_1(file, k, tree, numOfLines)
                    treelist, zipList = listTree(tree, root, '', [], [], k)
                                    
                # Else run the second routine
                else:
                    tree = routine_2(file, k , tree, numOfLines)
                    
                    # Forming treelist
                    treelist, zipList = listTree(tree, root, '', [], [], k)
                    
                    # Check if it is time to exit
                    if len(zipList) <= 1:
                        break
        
            # Printing Tree and forming necessary lists
            treelist = listSortBasedOnEvaluation(treelist)
            #treelist = treelist[0:topN]
            
            seqIndices = [treelist[i][4] for i in range(len(treelist))]
            timesPerSeq = [treelist[i][5] for i in range(len(treelist))]
            
            # Saving the output
            saveToFile.createCsvOutput(filename, treelist)
            saveToFile.createCsvOutputForSeqIndices(filename,seqIndices, timesPerSeq)

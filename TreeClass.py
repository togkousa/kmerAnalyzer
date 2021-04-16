# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 20:51:10 2019

@author: User
"""
import math
import featuresExtraction

eval_factor = featuresExtraction.eval_factor

mychars = ['A', 'C', 'G', 'T']

class Node:
    """
    Class Node
    """
    # Initializing a Node
    def __init__(self, value = '' , par = [], dep = 0, this_appended = False, this_sequenceIndex = '', this_timesPerSeq = ''):
        
        self.parent = par
        
        self.children = []
        self.deleted = []
        
        self.childrenchars = ['A', 'C', 'G', 'T'] if dep < 4 else []
        
        self.char = value
        self.count = 0
        self.evaluation = 0
        self.stop = False
        self.depth = dep
        self.appended = this_appended
        self.sequenceIndices = [this_sequenceIndex]
        self.timesPerSeq  = [this_timesPerSeq]

    # Deleting a child when evaluation is poor and appending the corresponding numeric value into a list
    def deleteChild(self, num):
        
        this_char = self.children[num].char
        
        del self.children[num]
        del self.childrenchars[num]
        self.deleted.append(this_char)
       
        if not self.children and len(self.deleted) == 4:
            self.stop = True
        
        return
    
    # Evaluating a node and deleting it when evaluation is poor
    def evalueateNode(self, num_kmers_scanned, k):
        
        expected = float(num_kmers_scanned) / pow(4,self.depth)
        self.evaluation = math.log(self.count / expected)
        
        return
    
    # Adding a new child
    def add_child(self, newNode):
        
        self.children.append(newNode)
        
        return
    
class Tree:
    """
    Class tree will provide a tree as well as utility functions.
    """
     # Initialization
    def __init__(self):
        
        self.root = Node()
    
    # Moving to child
    def move_to_child(self, current_node, direction):
        
        return current_node.children[direction]
    
    # Moving to parent
    def move_to_parent(self, current_node):
        
        return current_node.parent
  
    # Walking accross the path to find and evaluate the latter node
    def find_in_tree(self, kmer, mycheck , num_kmers_scanned, k, append_if_needed, sequenceIndex = ''):
        
        current = self.root
        check = mycheck
        just_created = False
        test = False
        
        for this_char in kmer:
            
            try:
                direction = current.childrenchars.index(this_char)
                current = self.move_to_child(current, direction)
                
                if current.depth == k-1  and append_if_needed and not kmer[-1] in current.childrenchars and not kmer[-1] in current.deleted:
                    current.add_child(Node(kmer[-1], current, current.depth+1, True, this_sequenceIndex=sequenceIndex, this_timesPerSeq = 1))
                    current.childrenchars.append(kmer[-1])
                    just_created = True
                    check = False
                    continue

                if current.depth == k:
                    current.count += 1
                    

                    if not just_created:
                        if current.sequenceIndices[-1] == sequenceIndex:
                            current.timesPerSeq[-1] += 1
                        else:
                            current.sequenceIndices.append(sequenceIndex)
                            current.timesPerSeq.append(1)
                    else:
                        break
                
            except ValueError:
                return
        
        
        if check and not current.stop:
            current.evalueateNode(num_kmers_scanned,k)
            if current.evaluation <= eval_factor*current.parent.evaluation:
                this_char = current.char
                par = self.move_to_parent(current)
                par.deleteChild(par.childrenchars.index(this_char))
        
        return

def del_list_inplace(l, id_to_del):
    for i in sorted(id_to_del, reverse=True):
        del l[i]
    return

# Deleting useless nodes
def check_tree(current, num_kmers_scanned, k):
    
    if current.stop or not current.children:
        return
    else:
        if current.depth == k-1:
            to_be_deleted = []
            for child in current.children:
                child.evalueateNode(num_kmers_scanned,k)
                if child.count == 1 or not child.char in mychars:
                    to_be_deleted.append(current.children.index(child))
                else:
                    if child.evaluation <= eval_factor*current.evaluation:
                        to_be_deleted.append(current.children.index(child))
            if to_be_deleted:
                del_list_inplace(current.children, to_be_deleted)
                del_list_inplace(current.childrenchars, to_be_deleted)            
            return
        else:
            for child in current.children:
                check_tree(child, num_kmers_scanned, k)
            
    return



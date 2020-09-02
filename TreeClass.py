# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 20:51:10 2019

@author: User
"""
import math

mychars = ['A', 'C', 'G', 'T']

class Node:
    """
    Class Node
    """
    # Initializing a Node
    def __init__(self, value = '' , par = [], dep = 0, this_appended = False, this_sequenceIndex = '', this_timesPerSeq = ''):
        
        self.parent = par
        
        self.children = []
        self.childrenchars = []
        self.deleted = []
        
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
        
        goto = [current_node.children[i].char for i in range(len(current_node.children))].index(mychars[direction])
        
        return current_node.children[goto]
    
    # Moving to parent
    def move_to_parent(self, current_node):
        
        return current_node.parent
  
    # Walking accross the path to find and evaluate the latter node
    def find_in_tree(self, path, mycheck , num_kmers_scanned, k, append_if_needed, sequenceIndex = ''):
        
        current = self.root
        check = mycheck
        just_created = False
        
        for direction in path:
            
            if not current.stop and current.children and mychars[direction] in [current.children[i].char for i in range(len(current.children))]:
                
                current = self.move_to_child(current, direction)
                if current.depth == k:
                    current.count += 1

                    if not just_created:
                        if current.sequenceIndices[-1] == sequenceIndex:
                            current.timesPerSeq[-1] += 1                            
                            #current.timesPerSeq[current.sequenceIndices.index(sequenceIndex)] += 1
                        else:
                            current.sequenceIndices.append(sequenceIndex)
                            current.timesPerSeq.append(1)
                    else:
                        break
                
                if current.depth == k-1  and append_if_needed and not mychars[path[-1]] in [current.children[u].char for u in range(len(current.children))] and not mychars[path[-1]] in current.deleted:
                    current.add_child(Node(mychars[path[-1]], current, current.depth+1, True, this_sequenceIndex=sequenceIndex, this_timesPerSeq = 1))
                    just_created = True
                    check = False
            else:
                return
                
        if check:   
            current.evalueateNode(num_kmers_scanned,k)
            if current.evaluation <= current.parent.evaluation:
                this_char = current.char
                par = self.move_to_parent(current)
                par.deleteChild([par.children[i].char for i in range(len(par.children))].index(this_char))
        
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
                if child.count == 1:
                    to_be_deleted.append(current.children.index(child))
            if to_be_deleted:
                del_list_inplace(current.children, to_be_deleted)
            for child in current.children:
                child.evalueateNode(num_kmers_scanned,k)
                if child.evaluation <= current.evaluation:
                    this_char = child.char
                    current.deleteChild([current.children[i].char for i in range(len(current.children))].index(this_char))

            return
        else:
            for child in current.children:
                check_tree(child, num_kmers_scanned, k)
            
    return



#!/usr/bin/env python3
import os
import sys
import re
import subprocess
import collections
import olc
from main import start, stop, s, o_min, max_length, readList


#----------------------------------------------------
# Graph class
#----------------------------------------------------
#Initializes a Graph object: a dictionary will be used for storing the nodes and their corresponding neighbouring nodes
class Graph:
    #Constructor
    def __init__(self, graph_dict):          
        self._graph = graph_dict

    #Accessor
    def _get_graph(self):
        '''Method to be call when we want to access the attribute "graph"'''
        return self._graph

    #Property
    graph = property(_get_graph)

    #Method "__getattr__"
    def __getattr__(self, attr):
        '''If Python doesn't find the attribute "attr", it calls this method and print an alert'''
        print("WARNING: There is no attribute {} here !".format(attr))

    #Method "__delattr_"
    def __delattr_(self, attr):
        '''We can't delete an attribute, we raise the exception AttributeError'''
        raise AttributeError("You can't delete attributes from this class")

    #Method "add_node
    def add_node(self, node):                   
        '''Method to add the node 'node' to the graph if it's not already in the graph'''
        if node not in self.graph:
            self.graph[node] = []

    #Method "add_edge"
    def add_edge(self, edge): 
        '''Method to add an 'edge' between a node and its neighbours'''        
        (source_node, neighbour, overlap) = edge
        #add an edge between the source_node and its neighbour node, with their corresponding overlap length
        if source_node in self.graph:
            self.graph[source_node].append([neighbour, overlap])                                                 
        else:
            self.graph[source_node] = [[neighbour, overlap]]

    #Method "nodes"
    def nodes(self):
        '''Method to return a list of the graph' nodes'''
        nodes = []    
        for node in self.graph.keys():
            nodes.append(node)                     
        return nodes

    #Method "edges"
    def edges(self):
        '''Method to return a list of the graph' edges, represented as a set, with one node (a loop back to the node) or two nodes and their overlapping length'''
        edges = []
        for node in self.graph:
            for neighbour in self.graph[node]:
                if [node, neighbour[0], neighbour[1]] not in edges:
                    edges.append([node, neighbour[0], neighbour[1]])
        return edges

    #Method "create_graph_from_extensions"
    def create_graph_from_extensions(self, source_node, extGroup):
        '''Method to create or update the graph from the extension groups (e.g. reads overlapping with the source node that share the same extension)'''
        '''NB: add only the first read for each extension group, e.g. the read having the larger overlap'''
        self.add_node(source_node)
        for reads in extGroup.values():
            self.add_node(reads[0][0])
            self.add_edge((source_node, reads[0][0], len(source_node)-reads[0][1]))

    #return all the paths from the start_node to the end_node
    def find_all_paths(self, start_node, end_node, path, all_paths):          
        if start_node not in self.graph or end_node not in self.graph:
            return []
        #path found from end_node to start_node
        if end_node == start_node:
            path.reverse()
            all_paths.append(path)
            return
        #traverse the graph to find the path from end_node to start_node
        prev_nodes = []
        for node in self.graph:
            for neighbour in self.graph[node]:
                if end_node in neighbour:
                    prev_nodes.append(node)
        for node in prev_nodes:
            self.find_all_paths(start_node, node, path+[node], all_paths)
        return all_paths

    #Method "__repr__"
    def __repr__(self):
        return "Nodes graph: {}".format(self.graph)




#----------------------------------------------------
# reverse_complement function
#----------------------------------------------------
'''
To reverse complement a sequence:
    - it takes as input the sequence we want to reverse complement
    - it outputs the reverse complement's sequence of the input sequence
'''
def reverse_complement(S):  
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'} 
    return ''.join([complement[base] for base in S[::-1]])


#----------------------------------------------------
# index_read function
#----------------------------------------------------
'''
To index a read by its seed:
    - it takes as input the sequence of the read, its position in readList (list containing all reads' sequences), the sequence of its reverse complement
      and the seedDict dictionary it will output
    - it outputs a dictionary seedDict: key = seed's sequence ; value = list of positions of reads having this seed in readList
'''
def index_read(read, i, read_rc, seedDict):                
    #Index reads by their seed
    seed = read[:s]
    if seed in seedDict:
        seedDict[seed].append(str(i))           
    else:
        seedDict[seed] = [str(i)]

    #Index reverse complement of reads by their seed as well
    seed = read_rc[:s]
    if seed in seedDict:
        seedDict[seed].append("-"+str(i))
    else:
        seedDict[seed] = ["-"+str(i)]


#----------------------------------------------------
# find_overlapping_reads function
#----------------------------------------------------
'''
To find the reads overlapping with the current assembly S sequence
    - it takes as input the current assembly's sequence (S) and the length of the read from which we want to extend
    - it outputs a list 'overlapping_reads' containing all the overlapping reads' sequences, along with the index of the beginning of the overlap, 
      referenced as [read's sequence, index of beginning of overlap]
      (NB: list sorted automatically by smallest i, e.g. by larger overlap and so by smallest extension)
'''
def find_overlapping_reads(S, len_read):
    overlapping_reads = []

    #Get the putative reads (e.g. reads having a seed onto the S sequence)
    for i in range(len(S)-len_read+1, len(S)-o_min-s):
        seed = S[i:i+s]
        if seed in olc.seedDict:
            putative_reads = olc.seedDict[seed]

            #For each putative read, search for an overlap between the S sequence and the putative read
            for put_read in putative_reads:
                nb_substitutions = 0
                l = i + s
                j = s
                length_overlap = 0

                #get the sequence of the read
                if '-' in str(put_read):
                    read = reverse_complement(readList[int(put_read.split('-')[1])])
                else:
                    read = readList[int(put_read)]

                while l < len(S) and j < len(read):
                    #match
                    if S[l] == read[j]:
                        l += 1
                        j += 1
                        length_overlap += 1
                    #mismatch
                    elif nb_substitutions < 2:              #error in reads: we allow 2 substitutions maximum
                        l += 1
                        j += 1
                        length_overlap += 1
                        nb_substitutions += 1
                    else:
                        break

                #Overlap found
                if l == len(S):
                    overlapping_reads.append([read, i])

    return overlapping_reads


#----------------------------------------------------
# extend function
#----------------------------------------------------
'''
To extend a read's sequence with overlapping reads
    - it takes as input the current assembly's sequence (S), the length of the read from which we want to extend and the abundance threshold value
    - it outputs the gap-filled sequence (S) if found / or the reason why the gap-filling failed, and a Boolean variable representing the success of the gapfilling
extGroup = dictionary containing the extension's sequence as key, and the reads sharing this extension as value 
(value format: [read's sequence, index of beginning of overlap])
'''
'''
def extend(S, read, a, seedDict, graph):
'''
def extend(S, len_read, a):
    #Base cases
    if stop in S[-len_read:]:
        '''
        graph.add_node(stop)
        graph.add_edge((read, stop, 0))
        '''
        return S, True

    if len(S) > max_length:
        return "\nAbundance threshold value: {} \n|S| > max_length".format(a), False

    #Search for reads overlapping with the current assembly S sequence
    overlapping_reads = find_overlapping_reads(S, len_read)
    '''
    overlapping_reads = find_overlapping_reads(S, len(read), seedDict)
    '''
    if len(overlapping_reads) == 0:
        return "\nAbundance threshold value: {} \nNo overlapping reads".format(a), False

    #Group the overlapping reads by their extension
    extGroup = {}

    #add the smallest extension to extGroup
    i = overlapping_reads[0][1]
    min_ext = overlapping_reads[0][0][len(S)-i:]
    extGroup[min_ext] = [overlapping_reads[0]]

    #populate extGroup
    '''NB: overlapping_reads list sorted automatically by smallest extension'''
    added_to_extGroup = True
    for (read_seq, index) in overlapping_reads[1:]:
        if len(extGroup) == 1:
            if read_seq[len(S)-index:len(S)-index+len(min_ext)] == min_ext:
                extGroup[min_ext].append([read_seq, index])
            else:
                extGroup[read_seq[len(S)-index:]] = [[read_seq, index]]
        elif len(extGroup) > 1:
            for extension in extGroup:
                if read_seq[len(S)-index:len(S)-index+len(extension)] == extension:
                    extGroup[extension].append([read_seq, index])
                    added_to_extGroup = True
                    break
                else:
                    added_to_extGroup = False
            if added_to_extGroup == False:
                extGroup[read_seq[len(S)-index:]] = [[read_seq, index]]

    #Filter extGroup by the number of reads sharing an extension (argument '-a')
    for extension in list(extGroup.keys()):
        if len(extGroup[extension]) < a:
            del extGroup[extension]

    if len(extGroup) == 0:
        return "\nAbundance threshold value: {} \nNo extension".format(a), False

    #Sort extGroup by the maximum overlap (e.g. by the minimal extension)
    extGroup = collections.OrderedDict(sorted(extGroup.items(), key=lambda t: len(t[0])))

    #Create graph "à la volée"
    '''
    graph.create_graph_from_extensions(read, extGroup)
    '''

    #Iterative extension of the assembly's sequence S
    for extension in extGroup:
        res, success = extend(S+extension, len(extGroup[extension][0][0]), a)
        '''
        res, success = extend(S+extension, extGroup[extension][0][0], a, seedDict, graph)
        '''
        if success:
            return res, True 
    return res, False




#TODO: check that if extension already in graph, do not gapfill again: actually, do gapfill again because the seq before ext could be different so not the same reads overlapping
#TODO: do not gapfill again if search on same region as one previously done (same nodes with same window length)
#TODO: use path to find all possible sequences (with only one stop in graph) ??
#TODO: add time for each step
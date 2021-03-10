#!/usr/bin/env python3
"""Module 'helpers.py': classes and functions of the script OLC

The module 'helpers.py' contains the classes and functions used in the script OLC.
"""

import collections
from Bio.Seq import Seq
from main import STOP, input_seqName, seed_size, min_overlap, list_of_abundance_min, max_length, max_subs, readList
from ProgDynOptim_SuffixPrefix import DynamicMatrixOptim


#----------------------------------------------------
# Graph class
#----------------------------------------------------
class Graph:
    """The class 'Graph' contains all the attributes, properties and methods to create a Graph object.

    The class 'Graph' initializes a Graph object.
    A dictionary will be used for storing the nodes and their corresponding neighbouring nodes.
    """
    # Constructor.
    def __init__(self, graph_dict):
        self._graph = graph_dict

    # Accessor.
    def _get_graph(self):
        '''Method to be call when we want to access the attribute "graph"'''
        return self._graph

    # Property.
    graph = property(_get_graph)

    # Method "__getattr__".
    def __getattr__(self, attr):
        '''If Python doesn't find the attribute "attr", it calls this method and print an alert'''
        print("WARNING: There is no attribute {} here !".format(attr))

    # Method "__delattr_".
    def __delattr_(self, attr):
        '''We can't delete an attribute, we raise the exception AttributeError'''
        raise AttributeError("You can't delete attributes from this class")

    # Method "add_node".
    def add_node(self, node):
        '''Method to add the node 'node' to the graph if it's not already in the graph'''
        if node not in self.graph:
            self.graph[node] = []

    # Method "add_edge".
    def add_edge(self, edge):
        '''Method to add an 'edge' between a node and its neighbours'''
        (source_node, neighbour, overlap) = edge
        #add an edge between the source_node and its neighbour node, with their corresponding overlap length
        if source_node in self.graph:
            self.graph[source_node].append([neighbour, overlap])
        else:
            self.graph[source_node] = [[neighbour, overlap]]

    # Method "nodes".
    def nodes(self):
        '''Method to return a list of the graph' nodes'''
        nodes = []
        for node in self.graph.keys():
            nodes.append(node)
        return nodes

    # Method "edges".
    def edges(self):
        '''Method to return a list of the graph' edges, represented as a set, with one node (a loop back to the node) or two nodes and their overlapping length'''
        edges = []
        for node in self.graph:
            for neighbour in self.graph[node]:
                if [node, neighbour[0], neighbour[1]] not in edges:
                    edges.append([node, neighbour[0], neighbour[1]])
        return edges

    # Method "create_graph_from_extensions".
    def create_graph_from_extensions(self, source_node, extGroup):
        '''Method to create or update the graph from the extension groups (e.g. reads overlapping with the source node that share the same extension)'''
        '''NB: add only the first read for each extension group, e.g. the read having the larger overlap'''
        self.add_node(source_node)
        for reads in extGroup.values():
            self.add_node(reads[0][0])
            self.add_edge((source_node, reads[0][0], len(source_node)-reads[0][1]))

    # Method "find_all_paths".
    def find_all_paths(self, start_node, end_node, path, all_paths):
        '''Method to return all the paths from the start_node to the end_node'''
        if start_node not in self.graph or end_node not in self.graph:
            return []
        # Path found from end_node to start_node.
        if end_node == start_node:
            path.reverse()
            all_paths.append(path)
            return
        # Traverse the graph to find the path from end_node to start_node.
        prev_nodes = []
        for node in self.graph:
            for neighbour in self.graph[node]:
                if end_node in neighbour:
                    prev_nodes.append(node)
        for node in prev_nodes:
            self.find_all_paths(start_node, node, path+[node], all_paths)
        return all_paths

    # Method "__repr__".
    def __repr__(self):
        return "Nodes graph: {}".format(self.graph)




#----------------------------------------------------
# index_read function
#----------------------------------------------------
def index_read(read, i, read_rc, seedDict):
    """To index a read by its seed.
    It updates a dictionary 'seedDict': key = seed's sequence ; value = list of positions of reads having this seed in readList

    Args:
        - read: str
            sequence of the current read to index
        - i: str
            position of the current read in readList (list containing all reads' sequences)
        - read_rc: str
            sequence of the reverse complement of the current read
        - seedDict: dict
            this function will output a dictionary: key = seed's sequence ; value = list of positions of reads having this seed in readList

    Outputs:
        - seedDict: dict
            dictionary of reads indexed by their seed: key = seed's sequence ; value = list of positions of reads having this seed in readList
    """
    # Index read by its seed.
    seed = read[:seed_size]
    if seed in seedDict:
        seedDict[seed].append(str(i))
    else:
        seedDict[seed] = [str(i)]

    # Index reverse complement of read by its seed as well.
    seed = read_rc[:seed_size]
    if seed in seedDict:
        seedDict[seed].append("-"+str(i))
    else:
        seedDict[seed] = ["-"+str(i)]


#----------------------------------------------------
# find_overlapping_reads function
#----------------------------------------------------
def find_overlapping_reads(assembly, len_read, seedDict):
    """
    To find the reads overlapping with the current assembly's sequence S
    The list 'overlapping_reads' it returns is sorted automatically by smallest i, e.g. by largest overlap

    Args:
        - assembly: str
            current assembly's sequence
        - len_read: int
            length of the read from which we want to extend
        - seedDict: dict
            dictionary of reads indexed by their seed: key = seed's sequence ; value = list of positions of reads having this seed in readList

    Returns:
        - overlapping_reads: list
            list containing all the overlapping reads' sequences, along with the index of the beginning of the overlap, and the index of the beginning of extension
            referenced as [read's sequence, index of beginning of overlap, index of beginning of extension]
    """
    overlapping_reads = []

    # Get the putative reads (e.g. reads having a seed onto the current assembly's sequence).
    for i in range(len(assembly)-len_read+1, len(assembly)-min_overlap-seed_size):
        seed = assembly[i:i+seed_size]
        if seed in seedDict:
            putative_reads = seedDict[seed]

            # For each putative read, search for an overlap between the current assembly's sequence and the putative read, using the dynamic programmation.
            for put_read in putative_reads:

                # Get the sequence of the read.
                if '-' in str(put_read):
                    read = str(Seq(readList[int(put_read.split('-')[1])]).reverse_complement())
                else:
                    read = readList[int(put_read)]
                
                # Perform the alignment with the optimized dynamic programmation.
                ##NB: Allow 4 gaps/substitutions
                dm = DynamicMatrixOptim(assembly[i:], read, 4, -4)
                dist, posR = dm.getEditDistanceAndGenomePosition()
                
                # Overlap found.
                if dist != None:
                    overlapping_reads.append([read, i, posR])

                #TODO: filtrer par dist ?
                #TODO: verifier qu'on a bien un alignement suffix-prefix et que donc l'alignement se termine quand on atteint len(G)
                    #Normalement renvoie 'None, None, None' si la fin de G et R ne s'alignent pas

    return overlapping_reads


#----------------------------------------------------
# extend function
#----------------------------------------------------
def extend(assembly, len_read, seedDict, assemblyHash):
    """
    To extend a read's sequence with overlapping reads
    The Boolean value it returns represents the success of the gap-filling
    NB: extGroup is a dictionary containing the extension's sequence as key, and the reads sharing this extension as value
        (value format: [read's sequence, index of beginning of overlap])
    If we use the 'graph' module: def extend(S, read, a, seedDict, graph):

    Args:
        - assembly: str
            current assembly's sequence
        - len_read: int
            length of the read from which we want to extend
        - seedDict: dict
            dictionary of reads indexed by their seed: key = seed's sequence ; value = list of positions of reads having this seed in readList
        - assemblyHash = hashtable/dict
            hashtable/dictionary indicating if the search for overlapping reads has already been performed on the corresponding sequence (key):
            key = the last 70 bp of the current assembly's sequence ; value = Boolean value (0: overlapping reads search not performed / 1: overlapping reads search performed)

    Returns:
        str, Boolean
            - the gap-filled sequence (assembly) and a Boolean variable equal to True if a solution is found (e.g. we arrived to STOP kmer)
            OR
            - the current assembly's sequence updated and a Boolean variable equal to False if no solution is found but we extended a little bit the 
              current assembly's sequence (e.g. we didn't arrive to STOP kmer)
            OR
            - the reason why the gap-filling failed and a Boolean variable equal to False if no solution is found and we didn't extended the current assembly's sequence
    """
    tmp_solutions = "tmp_solutions.fasta"

    # Base cases.
    if STOP in assembly[-len_read:]:
        '''
        graph.add_node(stop)
        graph.add_edge((read, stop, 0))
        '''
        return assembly, True

    if len(assembly) > max_length:
        return "\n|S| > max_length", False

    if len(assembly) >= 70:
        # Check that we didn't already search for overlapping reads on this region (e.g. on the last 70 bp of the current assembly's sequence).
        if assemblyHash[assembly[-70:]] == 1:
            return "\nPath already explored: No solution", False
            
    # Search for reads overlapping with the current assembly's sequence.
    overlapping_reads = find_overlapping_reads(assembly, len_read, seedDict)
    if not overlapping_reads:
        with open(tmp_solutions, "a") as tmp_file:
            tmp_file.write(">" + input_seqName + " _ No_read_overlapping")
            tmp_file.write("\n"+str(assembly)+"\n")
        return "\nNo overlapping reads", False

    # Group the overlapping reads by their extension.
    extGroup = {}

    # Populate extGroup.
    '''NB: overlapping_reads list sorted automatically by smallest i, e.g. by largest overlap'''
    for (read_seq, posS, posR) in overlapping_reads:

        # If no extension, don't add it to extGroup.
        if read_seq[posR:] == "":
            continue

        # Add first extension to extGroup.
        if len(extGroup) == 0:
            extGroup[read_seq[posR:]] = [[read_seq, posS]]

        # Add all extensions to extGroup.
        elif len(extGroup) > 0:
            for extension in extGroup:

                # Check that current extension is smaller than the one(s) in extGroup.
                if len(read_seq[posR:]) < len(extension):
                    # Current extension already in extGroup.
                    if read_seq[posR:] == extension[:len(read_seq[posR:])]:
                        new_extension = read_seq[posR:]
                        extGroup[new_extension] = extGroup[extension]
                        extGroup[new_extension].append([read_seq, posS])
                        del extGroup[extension]
                        added_to_extGroup = True
                        break
                    # Current extension is partially in extGroup.
                    elif read_seq[posR] == extension[0]:
                        i = 1
                        while i < len(read_seq[posR:]):
                            if read_seq[posR+i] == extension[i]:
                                i += 1
                            else:
                                break
                        new_extension = extension[:i]
                        extGroup[new_extension] = extGroup[extension]
                        extGroup[new_extension].append([read_seq, posS])
                        del extGroup[extension]
                        added_to_extGroup = True
                        break
                    # Current extension not already in extGroup.
                    else:
                        added_to_extGroup = False

                # Current extension is not smaller than the one(s) in extGroup.
                else:
                    # Current extension already in extGroup.
                    if read_seq[posR:posR+len(extension)] == extension:
                        extGroup[extension].append([read_seq, posS])
                        added_to_extGroup = True
                        break
                    # Current extension is partially in extGroup.
                    elif read_seq[posR] == extension[0]:
                        i = 1
                        while i < len(extension):
                            if read_seq[posR+i] == extension[i]:
                                i += 1
                            else:
                                break
                        new_extension = extension[:i]
                        extGroup[new_extension] = extGroup[extension]
                        extGroup[new_extension].append([read_seq, posS])
                        del extGroup[extension]
                        added_to_extGroup = True
                        break
                    # Current extension not already in extGroup.
                    else:
                        added_to_extGroup = False
                        
            # Current extension not already in extGroup.
            if not added_to_extGroup:
                extGroup[read_seq[posR:]] = [[read_seq, posS]]

        # Sort extGroup by the smallest extension.
        extGroup = collections.OrderedDict(sorted(extGroup.items(), key=lambda t: len(t[0])))

    # Update 'assemblyHash' to indicate that we performed the search for overlapping reads on this region (e.g. on the last 70 bp of the current assembly's sequence).
    assemblyHash[assembly[-70:]] = 1

    # Filter extGroup by the number of reads sharing an extension (argument 'abundance_min').
    for abundance_min in list_of_abundance_min:
        extGroup_filtered = extGroup.copy()
        for extension in list(extGroup_filtered.keys()):
            if len(extGroup_filtered[extension]) < abundance_min:
                del extGroup_filtered[extension]
        
        # Iterate over the values of abundance_min only if number of reads sharing an extension < 'abundance_min'.
        if not extGroup_filtered:
            continue
        else:
            break

    # If number of reads sharing an extension < minimal 'abundance_min' provided, stop the extension.
    if not extGroup_filtered:
        with open(tmp_solutions, "a") as tmp_file:
            tmp_file.write(">" + input_seqName + " _ No_extGroup")
            tmp_file.write("\n"+str(assembly)+"\n")
        return "\nNo extension", False

    # Sort extGroup by the extension whose read has the largest overlap with the current assembly's sequence (smallest i). 
    '''NB: values of extGroup sorted by reads having the larger overlap'''
    extGroup_filtered = collections.OrderedDict(sorted(extGroup_filtered.items(), key=lambda  t: t[1][0][1]))

    # Create graph "a la volee".
    '''
    graph.create_graph_from_extensions(read, extGroup)
    '''

    # Iterative extension of the assembly's sequence S.
    for extension in extGroup_filtered:
        
        # Update 'assemblyHash' with the new region for which we will search for overlapping reads (with value '0' if search not already performed, or with value '1' if search already performed).
        if (assembly+extension)[-70:] in assemblyHash.keys():
            assemblyHash[(assembly+extension)[-70:]] = 1
        else:
            assemblyHash[(assembly+extension)[-70:]] = 0
        
        res, success = extend(assembly+extension, len(extGroup_filtered[extension][0][0]), seedDict, assemblyHash)
        '''
        res, success = extend(assembly+extension, extGroup_filtered[extension][0][0], seedDict, graph)
        '''
        if success:
            return res, True
    return res, False

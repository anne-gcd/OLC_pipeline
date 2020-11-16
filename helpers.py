"""Module 'helpers.py': classes and functions of the script OLC

The module 'helpers.py' contains the classes and functions used in the script OLC.
"""
#!/usr/bin/env python3
import collections
from Bio.Seq import Seq
from main import STOP, seed_size, min_overlap, max_length, readList


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
def find_overlapping_reads(S, len_read, seedDict):
    """
    To find the reads overlapping with the current assembly's sequence S
    The list 'overlapping_reads' it returns is sorted automatically by smallest i, e.g. by larger overlap and so by smallest extension

    Args:
        - S: str
            current assembly's sequence
        - len_read: int
            length of the read from which we want to extend
        - seedDict: dict
            dictionary of reads indexed by their seed: key = seed's sequence ; value = list of positions of reads having this seed in readList

    Returns:
        - overlapping_reads: list
            list containing all the overlapping reads' sequences, along with the index of the beginning of the overlap,
            referenced as [read's sequence, index of beginning of overlap]
    """
    overlapping_reads = []

    # Get the putative reads (e.g. reads having a seed onto the S sequence).
    for i in range(len(S)-len_read+1, len(S)-min_overlap-seed_size):
        seed = S[i:i+seed_size]
        if seed in seedDict:
            putative_reads = seedDict[seed]

            # For each putative read, search for an overlap between the S sequence and the putative read.
            for put_read in putative_reads:
                nb_substitutions = 0
                l = i + seed_size
                j = seed_size
                length_overlap = 0

                # Get the sequence of the read.
                if '-' in str(put_read):
                    read = str(Seq(readList[int(put_read.split('-')[1])]).reverse_complement())
                else:
                    read = readList[int(put_read)]

                while l < len(S) and j < len(read):
                    # Match.
                    if S[l] == read[j]:
                        l += 1
                        j += 1
                        length_overlap += 1
                    # Mismatch (error in reads: we allow 2 substitutions maximum).
                    elif nb_substitutions < 2:
                        l += 1
                        j += 1
                        length_overlap += 1
                        nb_substitutions += 1
                    else:
                        break

                # Overlap found.
                if l == len(S):
                    overlapping_reads.append([read, i])

    return overlapping_reads


#----------------------------------------------------
# extend function
#----------------------------------------------------
def extend(S, len_read, a, seedDict):
    """
    To extend a read's sequence with overlapping reads
    The Boolean value it returns represents the success of the gap-filling
    NB: extGroup is a dictionary containing the extension's sequence as key, and the reads sharing this extension as value
        (value format: [read's sequence, index of beginning of overlap])
    If we use the 'graph' module: def extend(S, read, a, seedDict, graph):

    Args:
        - S: str
            current assembly's sequence
        - len_read: int
            length of the read from which we want to extend
        - a: int
            abundance threshold value (e.g. number minimal of reads sharing an extension)
        - seedDict: dict
            dictionary of reads indexed by their seed: key = seed's sequence ; value = list of positions of reads having this seed in readList

    Returns:
        str, Boolean
            - the gap-filled sequence (S) and a Boolean variable equal to True if a solution is found (e.g. we arrived to STOP kmer)
            OR
            - the current assembly's sequence updated and a Boolean variable equal to False if no solution is found but we extended a little bit S
              (e.g. we didn't arrive to STOP kmer)
            OR
            - the reason why the gap-filling failed and a Boolean variable equal to False if no solution is found and we didn't extended S
    """
    tmp_solutions = "tmp_solutions.txt"

    # Base cases.
    if STOP in S[-len_read:]:
        '''
        graph.add_node(stop)
        graph.add_edge((read, stop, 0))
        '''
        return S, True

    if len(S) > max_length:
        return "\nAbundance threshold value: {} \n|S| > max_length".format(a), False

    # Search for reads overlapping with the current assembly S sequence.
    overlapping_reads = find_overlapping_reads(S, len_read, seedDict)
    if not overlapping_reads:
        with open(tmp_solutions, "a") as tmp_file:
            tmp_file.write(">No_read_overlapping")
            tmp_file.write("\n"+str(S)+"\n")
        return "\nAbundance threshold value: {} \nNo overlapping reads".format(a), False

    # Group the overlapping reads by their extension.
    extGroup = {}

    # Add the smallest extension to extGroup.
    i = overlapping_reads[0][1]
    min_ext = overlapping_reads[0][0][len(S)-i:]
    extGroup[min_ext] = [overlapping_reads[0]]

    # Populate extGroup.
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
            if not added_to_extGroup:
                extGroup[read_seq[len(S)-index:]] = [[read_seq, index]]

    # Filter extGroup by the number of reads sharing an extension (argument '-a').
    for extension in list(extGroup.keys()):
        if len(extGroup[extension]) < a:
            del extGroup[extension]

    if not extGroup:
        with open(tmp_solutions, "a") as tmp_file:
            tmp_file.write(">No_extGroup")
            tmp_file.write("\n"+str(S)+"\n")
        return "\nAbundance threshold value: {} \nNo extension".format(a), False

    # Sort extGroup by the maximum overlap (e.g. by the minimal extension).
    extGroup = collections.OrderedDict(sorted(extGroup.items(), key=lambda t: len(t[0])))

    # Create graph "a la volee".
    '''
    graph.create_graph_from_extensions(read, extGroup)
    '''

    # Iterative extension of the assembly's sequence S.
    for extension in extGroup:
        res, success = extend(S+extension, len(extGroup[extension][0][0]), a, seedDict)
        '''
        res, success = extend(S+extension, extGroup[extension][0][0], a, seedDict, graph)
        '''
        if success:
            return res, True
    return res, False

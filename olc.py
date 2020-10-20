#!/usr/bin/env python3
from __future__ import print_function
import os
import sys
import re
import subprocess
import time
from Bio import SeqIO
from operator import itemgetter
from helpers import Graph, reverse_complement, index_read, extend
from main import start, stop, input_seqName, s, o_min, list_of_a, max_length, readList, assembly_file

#Augmenter la taille maximale de la pile de recursion en Python:
sys.setrecursionlimit(50000)


#----------------------------------------------------
# Gapfilling with Seed-and-Extend approach
#----------------------------------------------------
'''
readList = list of all reads' sequences
seedDict = dictionary containing the seed's sequence as key, and the list of positions of reads having this seed in readList as value (-pos if revcomp of read)
readWithStart = list of all reads containing the full sequence of the kmer start, along with the index of the beginning of the kmer start's subsequence, 
                referenced as a sublist of the readWithStart list: [position of the read in readList, index of beginning of kmer start's subsequence]
'''
try:
    startTime = time.time()

    #Gap-filling with a specific abundance threshold
    for a in list_of_a:

        seedDict = {}                                   #value of dict is a list of position (i) of reads in readList (-pos if revcomp of read)
        readWithStart = []
        pos_read_in_readList = 0

        for read in readList:
            #Get the reverse complement of the read
            read_rc = reverse_complement(read)

            #Seed the read and update the dictionary of the seeds of the reads
            index_read(read, pos_read_in_readList, read_rc, seedDict)

            #Search if the read contains the whole kmer start's sequence and update the readWithStart list
            if start in read:
                readWithStart.append([str(pos_read_in_readList), read.index(start)])
            elif start in read_rc:
                readWithStart.append(["-"+str(pos_read_in_readList), read_rc.index(start)])

            pos_read_in_readList += 1

        #Sort the readWithStart list by the minimum extension size (e.g. the maximum index)
        readWithStart = sorted(readWithStart, key=itemgetter(1), reverse=True)
        if len(readWithStart) == 0:
            print("\nNo read in the dataset provided contains the kmer start... \nHence, tentative of gapfilling aborted...")
            sys.exit(1)

        #Create graph "à la volée"
        '''
        nodes = {}
        graph = Graph(nodes)
        graph.add_node(start)
        '''

        #Extend the reads containing the whole kmer start's sequence
        for (pos_read, index) in readWithStart:
            #get the sequence of the read
            if '-' in str(pos_read):
                read = reverse_complement(readList[int(pos_read.split('-')[1])])
            else:
                read = readList[int(pos_read)]

            #update the graph with this read
            '''
            graph.add_node(read)
            graph.add_edge((start, read, 0))
            '''

            #Extend the assembly sequence
            res, success = extend(read, len(read), a, seedDict)
            '''
            res, success = extend(read, read, a, seedDict, graph)
            '''
            if not success:
                print(res)
            if success:
                print("\nAbundance threshold value: {} \nSuccessful Gapfilling !".format(a))
                #Save the gapfilled sequence in the output_file
                with open(assembly_file, "a") as assemblyFile:           #TODO: append or write ?
                    assembly_startbeg = res.index(start)
                    assembly_stopbeg = res.index(stop)
                    seq = res[assembly_startbeg:assembly_stopbeg+len(stop)]
                    seq_name = "assembly." + input_seqName + "_a" + str(a) + " len " + str(len(seq))
                    assemblyFile.write(">" + seq_name)
                    assemblyFile.write("\n" + seq + "\n")
                break

    endTime = time.time()
    print("It took %f seconds" %(endTime-startTime))


except Exception as e:
    print("\nException-")
    print(e)
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, fname, exc_tb.tb_lineno)
    sys.exit(1)




#TODO: if readWithStart empty, search for reads overlapping towards the right, with min overlap (especially if long start)
    #--> modif lines 47-49 of olc.py
#TODO: in future: try with multiprocessing

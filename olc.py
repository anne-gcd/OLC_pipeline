#!/usr/bin/env python3
"""Module 'olc.py': development of the script OLC

The module 'olc.py' contains the pipeline of the gap-filling using an OLC method.
Three main variables are used in this pipeline:
- readList = list of all reads' sequences
- seedDict1 / seedDict2 = dictionaries containing the seed's sequence as key, and the list of positions (i) of reads having this seed in readList as value (-pos if revcomp of read)
- readWithStart = list of all reads containing the full sequence of the kmer start, along with the index of the beginning of the kmer start's subsequence,
                referenced as a sublist of the readWithStart list: [position of the read in readList, index of beginning of kmer start's subsequence]
- assemblyHash = hashtable/dictionary containing the last 70 bp of the current assembly's sequence as key, and a Boolean value (indicating if the search for overlapping reads was already performed) as value
"""

from __future__ import print_function
import os
import sys
from operator import itemgetter
from Bio.Seq import Seq
from main import START, STOP, input_seqName, readList, assembly_file
from helpers import index_read, extend

# Increase the maximum recursion depth in Python.
sys.setrecursionlimit(50000)


#----------------------------------------------------
# Gapfilling with Seed-and-Extend approach
#----------------------------------------------------
try:
    # Initiate the four main variables.
    seedDict1 = {}
    seedDict2 = {}
    readWithStart = []
    pos_read_in_readList = 0
    assemblyHash = {}

    # Iterate over the reads of 'readList' to obtain the 'seedDict1' and 'seedDict2' dictionaries, as well as the 'readWithStart' list.
    for read in readList:
        # Get the reverse complement of the current read.
        read_rc = str(Seq(read).reverse_complement())
        # Seed the read and update the 'seedDict1' and 'seedDict2' dictionaries.
        index_read(read, pos_read_in_readList, read_rc, seedDict1, seedDict2)
        # Search if the read contains the whole kmer START's sequence and update the 'readWithStart' list.
        if START in read:
            readWithStart.append([str(pos_read_in_readList), read.index(START)])
        elif START in read_rc:
            readWithStart.append(["-"+str(pos_read_in_readList), read_rc.index(START)])
        # Increment the position of the current read in 'readList'
        pos_read_in_readList += 1

    # Sort the 'readWithStart' list by the minimum extension size (e.g. by the maximum index).
    readWithStart = sorted(readWithStart, key=itemgetter(1), reverse=True)
    # If there is no read containing the kmer start, raise an exception.
    if not readWithStart:
        print("\nNo read in the dataset provided contains the kmer start... \nHence, tentative of gapfilling aborted...")
        sys.exit(1)

    # Extend the reads containing the whole kmer start's sequence.
    for (pos_read, index) in readWithStart:

        # Get the sequence of the read.
        if '-' in str(pos_read):
            read = str(Seq(readList[int(pos_read.split('-')[1])]).reverse_complement())
        else:
            read = readList[int(pos_read)]

        # Extend the assembly sequence (e.g. the current read containing the whole kmer start's sequence) using the function 'extend()'
        assemblyHash[read[-70:]] = 0
        res, success = extend(read, len(read), seedDict1, seedDict2, assemblyHash)

        # Case of unsuccessful gap-filling.
        if not success:
            print(res)
        # Case of successful gap-filling.
        if success:
            print("\nSuccessful Gapfilling !")
            # Save the gap-filled sequence in the output_file.
            with open(assembly_file, "a") as assemblyFile:
                assembly_startbeg = res.index(START)
                assembly_stopbeg = res.index(STOP)
                seq = res[assembly_startbeg:assembly_stopbeg+len(STOP)]
                seq_name = "assembly." + input_seqName + " len_" + str(len(seq))
                assemblyFile.write(">" + seq_name)
                assemblyFile.write("\n" + seq + "\n")
            break


except Exception as exc:
    print("\nException-")
    print(exc)
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, fname, exc_tb.tb_lineno)
    sys.exit(1)

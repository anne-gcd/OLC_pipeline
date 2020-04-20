#!/usr/bin/env python3
from __future__ import print_function
import os
import sys
import argparse
import re
import time
from Bio import SeqIO


#----------------------------------------------------
# Arg parser
#----------------------------------------------------
parser = argparse.ArgumentParser(prog="olc.py", usage="%(prog)s -in <input_sequence> -reads <reads_file> -s <seed_size> -o <minimum_overlap_size> -l <maximum_assembly_length> [options]", \
                                formatter_class=argparse.RawTextHelpFormatter, \
                                description=("Gapfilling, using an Overlap-Layout-Consensus (OLC) method"))

parser.add_argument('-in', action="store", dest="input", help="input sequence to gapfill (for example, kmers start and stop)", required=True)
parser.add_argument('-reads', action="store", dest="reads", help="file of reads", required=True)
parser.add_argument('-s', action="store", dest= "seed_size", type=int, help="seed size used for indexing the reads (bp)", required=True)
parser.add_argument('-o', action="store", dest="min_overlap", type=int, help="minimum overlapping size (bp)", required=True)
parser.add_argument('-a', action="store", dest="abundance_min", nargs='*', type=int, default=2, help="minimal abundance(s) of reads used for gapfilling ; extension's groups having less than this number of reads are discarded from the graph")
parser.add_argument('-l', action="store", dest="max_length", type=int, help="maximum assembly length (bp) (it could correspond to the length of the gap to fill (+length input sequences) OR it could be a very high length to prevent for searching indefinitely", required=True)
parser.add_argument('-out', action="store", dest="outdir", default="./olc_results", help="output directory for the results' files")

args = parser.parse_args()

if not re.match('^.*.fasta$', args.input) and not re.match('^.*.fa$', args.input):
    parser.error("The input file should be a FASTA file")
        
if not re.match('^.*.fasta$', args.reads) and not re.match('^.*.fa$', args.reads) and not re.match('^.*.fastq$', args.reads) and not re.match('^.*.fq$', args.reads):
    parser.error("The reads file should be a FASTA or FASTQ file")

#----------------------------------------------------
# Input files
#----------------------------------------------------
input_file = os.path.abspath(args.input)
if not os.path.exists(input_file):
    parser.error("\nThe path of the input file doesn't exist")
print("\nInput file: " + input_file)

reads_file = os.path.abspath(args.reads)
if not os.path.exists(reads_file):
    parser.error("The path of the reads' file doesn't exist")
print("Reads' file: " + reads_file)

#Get the inputs' sequences
with open(input_file, "r") as inputFile:
    for record in SeqIO.parse(inputFile, "fasta"):
        if record.id == "start" or "left" in record.description:
            start = str(record.seq)
            input_seqName = record.id
        if record.id == "stop" or "right" in record.description:
            stop = str(record.seq)
            if record.id != input_seqName:
                input_seqName += "-" + record.id

#----------------------------------------------------
# Directories for saving results
#----------------------------------------------------
cwd = os.getcwd() 

#outDir
if not os.path.exists(args.outdir):
    os.mkdir(args.outdir)
try:
    os.chdir(args.outdir)
except:
    print("Something wrong with specified directory. Exception-", sys.exc_info())
    print("Restoring the path")
    os.chdir(cwd)
outDir = os.getcwd()
print("\nThe results are saved in " + outDir)

#----------------------------------------------------
# Parameters
#----------------------------------------------------
s = args.seed_size
o_min = args.min_overlap
list_of_a = args.abundance_min
max_length = args.max_length

#----------------------------------------------------
# Output files for saving results
#----------------------------------------------------
#FASTA file containing all possible gapfilled sequences
output_file = "assembly.s{}.o{}.olc_gapfilling.fasta".format(s, o_min)
assembly_file = os.path.abspath(outDir +"/"+ output_file)

#----------------------------------------------------
# Save reads' sequences in a list
#----------------------------------------------------
'''
readList = list of all reads' sequences
'''
readList = []

with open(reads_file, "r") as readsFile:
    if re.match('^.*.fasta$', reads_file) or re.match('^.*.fa$', reads_file):
        readList = [str(read.seq) for read in SeqIO.parse(readsFile, "fasta")]
    elif re.match('^.*.fastq$', reads_file) or re.match('^.*.fq$', reads_file):
        readList = [str(read.seq) for read in SeqIO.parse(readsFile, "fastq")]




#TODO: take GFA as input
#TODO: input seq instead of kmers start and stop: get kmers start and stop from input seq to gapfill !
    #--> modif lines 59-68 of main.py
#TODO: tester sur format data ecoli
#TODO: tester sur format data supergene (same as for MTG-Link)
#TODO: try multiprocessing and see if runtime better (but also compare memory used)

"""Module 'main.py': initialization of the script OLC

The module 'main.py' enables to get the input parameters and creates the file and directory in which to save the results.
It creates as well the list 'readList' containing all reads' sequences.
"""
#!/usr/bin/env python3
from __future__ import print_function
import argparse
import os
import re
import sys
from Bio import SeqIO


#----------------------------------------------------
# Arg parser
#----------------------------------------------------
parser = argparse.ArgumentParser(prog="olc.py", usage="%(prog)s -in <input_sequences> -reads <reads_file> -s <seed_size> -o <minimum_overlap_size> -l <maximum_assembly_length> [options]", \
                                formatter_class=argparse.RawTextHelpFormatter, \
                                description=("Gapfilling, using an Overlap-Layout-Consensus (OLC) method"))

parser.add_argument('-in', action="store", dest="input", help="Input sequences to gapfill (for example, kmers start and stop)", required=True)
parser.add_argument('-reads', action="store", dest="reads", help="File of reads", required=True)
parser.add_argument('-s', action="store", dest="seed_size", type=int, help="Seed size used for indexing the reads (bp)", required=True)
parser.add_argument('-o', action="store", dest="min_overlap", type=int, help="Minimum overlapping size (bp)", required=True)
parser.add_argument('-a', action="store", dest="abundance_min", nargs='*', type=int, default=2, help="Minimal abundance(s) of reads used for gapfilling ; extension's groups having less than this number of reads are discarded from the graph")
parser.add_argument('-l', action="store", dest="max_length", type=int, help="Maximum assembly length (bp) (it could correspond to the length of the gap to fill (+length input sequences) OR it could be a very high length to prevent for searching indefinitely", required=True)
parser.add_argument('-out', action="store", dest="outdir", default="./olc_results", help="Output directory for the results' files")

args = parser.parse_args()

if not re.match('^.*.fasta$', args.input) and not re.match('^.*.fa$', args.input):
    parser.error("The input file should be a FASTA file.")
if not re.match('^.*.fasta$', args.reads) and not re.match('^.*.fa$', args.reads) and not re.match('^.*.fastq$', args.reads) and not re.match('^.*.fq$', args.reads):
    parser.error("The reads file should be a FASTA or FASTQ file.")

#----------------------------------------------------
# Input files
#----------------------------------------------------
# Get the file of input sequences to gap-fill.
input_file = os.path.abspath(args.input)
if not os.path.exists(input_file):
    parser.error("\nThe path of the input file doesn't exist.")
print("\nInput file: " + input_file)

# Get the reads file (FASTA or FASTQ).
reads_file = os.path.abspath(args.reads)
if not os.path.exists(reads_file):
    parser.error("The path of the reads' file doesn't exist.")
print("Reads' file: " + reads_file)

# Get the inputs' sequences.
with open(input_file, "r") as inputFile:
    for record in SeqIO.parse(inputFile, "fasta"):
        if record.id == "start" or "left" in record.description:
            START = str(record.seq)
            input_seqName = record.id
        if record.id == "stop" or "right" in record.description:
            STOP = str(record.seq)
            if record.id != input_seqName:
                input_seqName += "-" + record.id

#----------------------------------------------------
# Directories for saving results
#----------------------------------------------------
cwd = os.getcwd()

# Create the directory 'outDir', where the results will be saved.
if not os.path.exists(args.outdir):
    os.mkdir(args.outdir)
try:
    os.chdir(args.outdir)
except OSError:
    print("Something wrong with specified directory. Exception-", sys.exc_info())
    print("Restoring the path.")
    os.chdir(cwd)
outDir = os.getcwd()
print("\nThe results are saved in " + outDir)

#----------------------------------------------------
# Parameters
#----------------------------------------------------
seed_size = args.seed_size
min_overlap = args.min_overlap
list_of_abundance_min = args.abundance_min
max_length = args.max_length

#----------------------------------------------------
# Output file for saving results
#----------------------------------------------------
# Create the FASTA file containing all possible gapfilled sequences.
output_file = "assembly.s{}.o{}.olc_gapfilling.fasta".format(seed_size, min_overlap)
assembly_file = os.path.abspath(outDir +"/"+ output_file)

#----------------------------------------------------
# Save reads' sequences in a list
#----------------------------------------------------
# Create the list 'readList' containing all reads' sequences.
readList = []

with open(reads_file, "r") as readsFile:
    if re.match('^.*.fasta$', reads_file) or re.match('^.*.fa$', reads_file):
        readList = [str(read.seq) for read in SeqIO.parse(readsFile, "fasta")]
    elif re.match('^.*.fastq$', reads_file) or re.match('^.*.fq$', reads_file):
        readList = [str(read.seq) for read in SeqIO.parse(readsFile, "fastq")]

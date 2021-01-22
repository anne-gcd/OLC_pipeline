#!/usr/bin/env python3
from __future__ import print_function
import os
import sys
import argparse
import csv
import re
import subprocess
from pathos.multiprocessing import ProcessingPool as Pool
#from multiprocessing import Pool
import gfapy
from gfapy.sequence import rc
from Bio import SeqIO, Align
from helpers_pipeline import Gap, Scaffold, extract_barcodes, get_reads, stats_align, get_position_for_edges, get_output_for_gfa, update_gfa_with_solution


#----------------------------------------------------
# Arg parser
#----------------------------------------------------
parser = argparse.ArgumentParser(prog="olc_pipeline.py", usage="%(prog)s -gfa <input.gfa> -c <chunk_size> -bam <mapped.bam> -fastq <reads.fastq> -index <barcoded.shelve> [options]", \
                                description=("Gapfilling with linked read data, using MindTheGap in 'breakpoint' mode"))

parserMain = parser.add_argument_group("[Main options]")
parserOLC = parser.add_argument_group("[OLC option]")

parserMain.add_argument('-gfa', dest="input_gfa", action="store", help="Input GFA file (GFA 2.0) (format: xxx.gfa)", required=True)
parserMain.add_argument('-c', dest="chunk", action="store", type=int, help="Chunk size (bp)", required=True)
parserMain.add_argument('-bam', dest="bam", action="store", help="BAM file: linked reads mapped on current genome assembly (format: xxx.bam)", required=True)
parserMain.add_argument('-fastq', dest="reads", action="store", help="File of indexed reads (format: xxx.fastq | xxx.fq)", required=True)
parserMain.add_argument('-index', dest="index", action="store", help="Prefix of barcodes index file (format: xxx.shelve)", required=True)
parserMain.add_argument('-f', dest="freq", action="store", type=int, default=2, help="Minimal frequence of barcodes extracted in the chunk of size '-c' [default: 2]")
parserMain.add_argument('-out', dest="outDir", action="store", default="./olc_gapfilling", help="Output directory for the result's files [default './olc_gapfilling']")
parserMain.add_argument('-refDir', dest="refDir", action="store", help="Directory containing the reference sequences if any")
parserMain.add_argument('-line', dest="line", action="store", type=int, help="Line of GFA file input from which to start analysis (if not provided, start analysis from first line of GFA file input) [optional]")
parserMain.add_argument('-rbxu', dest="rbxu", action="store", help="Files containing the reads of the union of the corresponding gaps (if already extracted) [optional]")

parserOLC.add_argument('-s', dest="seed_size", action="store", type=int, help="Seed size used for indexing the reads (bp)", required=True)
parserOLC.add_argument('-o', dest="min_overlap", action="store", type=int, help="Minimum overlapping size (bp)", required=True)
parserOLC.add_argument('-a', dest="abundance_min", action="store", nargs='*', type=int, default=[3, 2], help="Minimal abundance(s) of reads used for gapfilling ; extension's groups having less than this number of reads are discarded from the graph")
parserOLC.add_argument('-ext', dest="extension", action="store", type=int, default=500, help="Extension size of the gap on both sides (bp); determine start/end of gapfilling [default: '500']")
parserOLC.add_argument('-l', dest="max_length", action="store", type=int, help="Maximum assembly length (bp) (it could correspond to the length of the gap to fill (+ length START/STOP + 2*ext) OR it could be a very high length to prevent for searching indefinitely", required=True)
parserOLC.add_argument('-subs', dest="max_subs", action="store", type=int, default=2, help="Maximum number of substitutions allowed in the inexact overlap between reads")

args = parser.parse_args()

if re.match('^.*.gfa$', args.input_gfa) is None:
    parser.error("Warning: The suffix of the GFA file should be: '.gfa'")

if re.match('^.*.bam$', args.bam) is None:
    parser.error("Warning: The suffix of the BAM file should be: '.bam'")

#----------------------------------------------------
# Input files and arguments
#----------------------------------------------------
#GFA 2.0 file
gfa_file = os.path.abspath(args.input_gfa)
if not os.path.exists(gfa_file):
    parser.error("Warning: The path of the GFA file doesn't exist")
gfa_name = gfa_file.split('/')[-1]
print("\nInput GFA file: " + gfa_file)

#BAM file: linked reads mapped on current genome assembly
bam_file = os.path.abspath(args.bam)
if not os.path.exists(bam_file): 
    parser.error("Warning: The path of the BAM file doesn't exist")
print("BAM file: " + bam_file)

#Reads file: file of indexed reads
reads_file = os.path.abspath(args.reads)
if not os.path.exists(reads_file):
    parser.error("Warning: The path of the file of indexed reads doesn't exist")
print("File of indexed reads: " + reads_file)

#Prefix of barcodes index file
index_file = os.path.abspath(args.index)
print("Barcodes index file (prefix): " + index_file)

#Directory containing the reference sequences if any
if args.refDir is not None:
    refDir = os.path.abspath(args.refDir)
    if not os.path.exists(refDir):
        parser.error("Warning: The path of the directory containing the reference sequences doesn't exist")

#Directory containing the union' reads files if any
if args.rbxu is not None:
    rbxuDir = os.path.abspath(args.rbxu)
    if not os.path.exists(rbxuDir):
        parser.error("Warning: The path of the directory containing the union' reads files doesn't exist")

#variable 'ext' is the size of the extension of the gap, on both sides [by default 500]
ext = args.extension


#----------------------------------------------------
# Directories for saving results
#----------------------------------------------------
cwd = os.getcwd() 

#outDir
if not os.path.exists(args.outDir):
    os.mkdir(args.outDir)
try:
    os.chdir(args.outDir)
except:
    print("Something wrong with specified directory. Exception-", sys.exc_info())
    print("Restoring the path")
    os.chdir(cwd)
outDir = os.getcwd()
print("\nThe results are saved in " + outDir + "\n")

#unionDir
unionDir = outDir + "/union"
os.mkdir(unionDir)

#olcDir
olcDir = outDir + "/olc_results"
os.mkdir(olcDir)

#contigDir
contigDir = outDir + "/contigs"
os.mkdir(contigDir)

#statsDir
statsDir = outDir + "/alignments_stats"


#----------------------------------------------------
# gapfilling function - Pipeline
#----------------------------------------------------
'''
To perform the gap-filling on a specific gap:
    - it takes as input the current gap on which we want to perform the gap-filling
    - it outputs the list 'union_summary' containing the gap ID, the names of the left and right flanking sequences, the gap size, the chunk size, and the number of barcodes and reads extracted on the chunks to perform the gap-filling
    - it outputs as well the list 'output_for_gfa' containing the gap-filled sequence's name, as well as its length, its sequence, the number of solution found, the beginning and ending positions of the overlap and the quality of the sequence
'''
def gapfilling(current_gap):

    os.chdir(outDir)

    #Open the input GFA file to get the corresponding Gap line ('G' line)
    gfa = gfapy.Gfa.from_file(gfa_file)
    for _gap_ in gfa.gaps:
        if str(_gap_) == current_gap:
            current_gap = _gap_
            #Create the object 'gap' from the class 'Gap'
            gap = Gap(current_gap)

    #Get some information on the current gap we are working on
    gap.info()
    gap_label = gap.label()

    #Create two objects ('left_scaffold' and 'right_scaffold') from the class 'Scaffold'
    left_scaffold = Scaffold(current_gap, gap.left, gfa_file)
    right_scaffold = Scaffold(current_gap, gap.right, gfa_file)

    #If chunk size larger than length of scaffold(s), set the chunk size to the minimal scaffold length
    #chunk_L
    if args.chunk > left_scaffold.slen:
        print("Warning for {}: The chunk size you provided is higher than the length of the left scaffold. Thus, for the left scaffold, the barcodes will be extracted on its whole length".format(gap_label))
        chunk_L = left_scaffold.slen
    else:
        chunk_L = args.chunk
    #chunk_R
    if args.chunk > right_scaffold.slen:
        print("Warning for {}: The chunk size you provided is higher than the length of the right scaffold. Thus, for the right scaffold, the barcodes will be extracted on its whole length".format(gap_label))
        chunk_R = right_scaffold.slen
    else:
        chunk_R = args.chunk

    #----------------------------------------------------
    # BamExtractor
    #----------------------------------------------------
    #Union output directory
    os.chdir(unionDir)
    
    #Initiate a dictionary to count the occurences of each barcode
    barcodes_occ = {}
    
    #Obtain the left barcodes that are extracted on the left region and store the barcodes and their occurences in the dict 'barcodes_occ'
    left_region = left_scaffold.chunk(chunk_L)
    extract_barcodes(bam_file, gap_label, left_region, barcodes_occ)

    #Obtain the right barcodes that are extracted on the right region and store the barcodes and their occurences in the dict 'barcodes_occ'
    right_region = right_scaffold.chunk(chunk_R)
    extract_barcodes(bam_file, gap_label, right_region, barcodes_occ)

    #Do the union of the barcodes on both left and right regions
    union_barcodes_file = "{}.{}.g{}.c{}.bxu".format(gfa_name, str(gap_label), gap.length, args.chunk)
    with open(union_barcodes_file, "w") as union_barcodes:
        #Filter barcodes by freq
        for (barcode, occurences) in barcodes_occ.items():
            if occurences >= args.freq:
                union_barcodes.write(barcode + "\n")

    #----------------------------------------------------
    # GetReads
    #----------------------------------------------------
    #If the reads of the union are already extracted, use the corresponding file
    if args.rbxu is not None:
        for file_ in os.listdir(rbxuDir):
            if (re.match('^.*.rbxu.fastq$', str(file_))) and (str(gap_label) in file_):
                union_reads_file = rbxuDir +"/"+ str(file_)
        if not os.path.isfile(union_reads_file):
            print("Warning: No union' reads file was found for this gap...")
    
    #Union: extract the reads associated with the barcodes
    else:
        union_reads_file = "{}.{}.g{}.c{}.rbxu.fastq".format(gfa_name, str(gap_label), gap.length, args.chunk)
        with open(union_reads_file, "w") as union_reads:
            get_reads(reads_file, index_file, gap_label, union_barcodes_file, union_reads)

    #----------------------------------------------------
    # Summary of union (barcodes and reads)
    #----------------------------------------------------
    bxu = sum(1 for line in open(union_barcodes_file, "r"))
    rbxu = sum(1 for line in open(union_reads_file, "r"))/4
    union_summary = [str(gap.identity), str(gap.left), str(gap.right), gap.length, args.chunk, bxu, rbxu]

    #Remove the barcodes files
    subprocess.run(["rm", union_barcodes_file])

    #----------------------------------------------------
    # OLC Gap-Filling
    #----------------------------------------------------        
    #Get flanking contigs sequences
    seq_L = str(left_scaffold.sequence())
    seq_R = str(right_scaffold.sequence())

    #Execute the OLC module on the union
    os.chdir(olcDir)

    #----------------------------------------------------
    # Input FASTA file, with offset of size ext removed
    #----------------------------------------------------
    #Get a FASTA file containing the input sequences (start and stop kmers)
    input_file = "{}.{}.g{}.c{}.start_stop.fa".format(gfa_name, str(gap_label), gap.length, args.chunk)
    with open(input_file, "w") as input_olc:

        #Start sequence
        line1 = ">ctg{}_start _ len_31_bp (left)\n".format(left_scaffold.scaffold)
        line2 = seq_L[(left_scaffold.slen - ext - 31):(left_scaffold.slen - ext)]
        
        #Stop sequence
        line3 = "\n>ctg{}_stop _ len_31_bp (right)\n".format(right_scaffold.scaffold)
        line4 = seq_R[ext:(ext + 31)]

        input_olc.writelines([line1, line2, line3, line4])

    #----------------------------------------------------
    # Gap-filling
    #----------------------------------------------------
    print("\nGap-filling of {} (union)".format(str(gap_label)))

    #Input arguments for OLC
    input_reads_file = os.path.join(unionDir, union_reads_file)
    seed_size = args.seed_size
    min_overlap = args.min_overlap
    list_of_abundance_min = args.abundance_min
    str_of_abundance_min = ' '.join(map(str, list_of_abundance_min))
    max_length = args.max_length
    max_subs = args.max_subs
    olc_outDir = "./s{}o{}".format(seed_size, min_overlap)
    output_file = "{}.g{}.c{}.s{}.o{}.olc_gapfilling.fasta".format(str(gap_label), gap.length, args.chunk, seed_size, min_overlap)

    #Perform the gap-filling with OLC
    olc_command = str(sys.path[0]) + "/olc.py -in " + input_file + " -reads " + input_reads_file + " -s " + str(seed_size) + " -o " + str(min_overlap) + " -a " + str_of_abundance_min + " -l " + str(max_length) + " -subs " + str(max_subs) + " -out " + olc_outDir + " -assembly " + output_file
    olcLog = str(gap_label) + "_olc.log"

    with open(olcLog, "a") as log:
        subprocess.run(olc_command, shell=True, stderr=log, stdout=log)

    #If one solution is found, perform qualitative evaluation of the gap-filled sequence(s)
    assembly_file = os.path.abspath(olcDir +"/"+ olc_outDir +"/"+ output_file)
    if os.path.exists(assembly_file):
        
        #----------------------------------------------------
        # Stats of the alignments query_seq vs reference_seq
        #----------------------------------------------------
        #Qualitative evaluation with the reference sequence
        if args.refDir is not None:
            for file_ in os.listdir(refDir):
                if str(gap_label) in file_:
                    ref_file = refDir +"/"+ str(file_)
            if not os.path.isfile(ref_file):
                print("Warning: No reference file was found for this gap. The qualitative evaluation will be performed with the flanking contigs information.")
    
        #Qualitative evalution with the flanking contigs information
        elif (args.refDir is None) or (ref_file is None):

            #Merge both left and right flanking contigs sequences into a unique file (ref_file)
            ref_file = contigDir +"/"+ str(gap_label) +".g"+ str(gap.length) + ".contigs.fasta"
            with open(ref_file, "w") as ref_fasta:

                #Left scaffold oriented '+'
                if left_scaffold.orient == "+":
                    ref_fasta.write(">" + left_scaffold.name + "_region:" + str(left_scaffold.slen-ext) + "-" + str(left_scaffold.slen) + "\n")
                    ref_fasta.write(seq_L[(left_scaffold.slen - ext):left_scaffold.slen])
                #Left scaffold oriented '-' ~ Right scaffold oriented '+'
                elif left_scaffold.orient == "-":
                    ref_fasta.write(">" + left_scaffold.name + "_region:0-" + str(ext) + "\n")
                    ref_fasta.write(str(rc(seq_L)[0:ext]))

                #Right scaffold oriented '+'
                if right_scaffold.orient == "+":
                    ref_fasta.write("\n>" + right_scaffold.name + "_region:0-" + str(ext) + "\n")
                    ref_fasta.write(seq_R[0:ext])
                #Right scaffold oriented '-' ~ Left scaffold oriented '+'
                elif right_scaffold.orient == "-":
                    ref_fasta.write("\n>" + right_scaffold.name + "_region:" + str(right_scaffold.slen-ext) + "-" + str(right_scaffold.slen) + "\n")
                    ref_fasta.write(str(rc(seq_R)[(right_scaffold.slen - ext):right_scaffold.slen]))

        if not os.path.isfile(ref_file):
            print("Warning: Something wrong with the specified reference file. Exception-", sys.exc_info())

        #Do statistics on the alignments of query_seq (found gapfill seq) vs reference
        else:
            prefix = "{}.s{}.o{}".format(str(gap_label), seed_size, min_overlap)
            stats_align(gap_label, assembly_file, ref_file, str(ext), prefix, statsDir)


        #----------------------------------------------------
        # Estimate quality of gapfilled sequence
        #----------------------------------------------------
        #Reader for alignment stats' files
        ref_qry_file = statsDir + "/" + prefix + ".ref_qry.alignment.stats"

        if not os.path.exists(ref_qry_file):
            print("Warning: The '{}' file doesn't exits".format(ref_qry_file))

        else:
            ref_qry_output = open(ref_qry_file)

            reader_ref_stats = csv.DictReader(ref_qry_output, \
                                            fieldnames=("Gap", "Len_gap", "Chunk", "Seed_size", "Min_overlap", "Len_Q", "Ref", "Len_R", \
                                                        "Start_ref", "End_ref", "Start_qry", "End_qry", "Len_alignR", "Len_alignQ", "%_Id", "%_CovR", "%_CovQ", "Frame_R", "Frame_Q", "Quality"), \
                                            delimiter='\t')

            #Obtain a quality score for each gapfilled seq
            output_for_gfa = []
            assembly_quality_file = os.path.abspath(olcDir +"/"+ olc_outDir +"/"+ output_file.split('.fasta')[0] + ".quality.fasta")
            bad_solutions_file = os.path.abspath(outDir + "/bad_solutions.fasta")

            with open(assembly_file, "r") as query, open(assembly_quality_file, "w") as qualified:
                for record in SeqIO.parse(query, "fasta"):
                    seq = record.seq
                    record_label = (record.id).split("assembly.ctg")[1].split('_start')[0] +"_"+ (record.id).split("-ctg")[1].split('_stop')[0]

                    #----------------------------------------------------
                    #Ref = reference sequence of simulated gap
                    #----------------------------------------------------
                    if args.refDir is not None:
                        #quality score for stats about the ref
                        quality_ref = []
                        for row in reader_ref_stats:
                            if (row["Gap"] == record_label):
                                quality_ref.append(row["Quality"])
                        
                        if quality_ref == []:
                            quality_ref.append('D')

                        ref_qry_output.seek(0)

                        #global quality score
                        quality_gapfilled_seq = min(quality_ref)
                        
                        record.description = "Quality " + str(quality_gapfilled_seq)
                        SeqIO.write(record, qualified, "fasta")

                        #Update GFA with only the good solutions (the ones having a good quality score)
                        if (len(seq) > 2*ext) and (re.match('^.*Quality [AB]$', record.description)):
                            gfa_output = get_output_for_gfa(record, ext, seed_size, min_overlap, gap.left, gap.right, left_scaffold, right_scaffold)
                            output_for_gfa.append(gfa_output)
                        #Add the bad solutions to a FASTA file containing all bad solutions
                        else:
                            output_for_gfa = []     #Pbm if several solutions output
                            with open(bad_solutions_file, "a") as bad_file:
                                bad_file.write("\n>" + str(gap.left) +":"+ str(gap.right) + "_gf.s" + str(seed_size) + ".o" + str(min_overlap) + "_fwd _ len_" + str(len(seq)) + "_qual_" + str(quality_gapfilled_seq) +"\n")
                                bad_file.write(str(seq))

                    #----------------------------------------------------
                    #Ref = flanking contigs' sequences
                    #----------------------------------------------------
                    else:
                        #quality score for stats about the extension
                        quality_ext_left = []
                        quality_ext_right = []
                        for row in reader_ref_stats:
                            if (row["Gap"] == record_label) and (row["Ref"] == left_scaffold.name):
                                quality_ext_left.append(row["Quality"])
                            elif (row["Gap"] == record_label) and (row["Ref"] == right_scaffold.name):
                                quality_ext_right.append(row["Quality"])
                        if quality_ext_left == []:
                            quality_ext_left.append('D')
                        if quality_ext_right == []:
                            quality_ext_right.append('D')

                        ref_qry_output.seek(0)

                        #global quality score
                        quality_gapfilled_seq = min(quality_ext_left) + min(quality_ext_right)

                        record.description = "Quality " + str(quality_gapfilled_seq)
                        SeqIO.write(record, qualified, "fasta")

                        #Update GFA with only the good solutions (the ones having a good quality score)
                        if (len(seq) > 2*ext) and (re.match('^.*Quality [AB]{2}$', record.description)):
                            gfa_output = get_output_for_gfa(record, ext, seed_size, min_overlap, gap.left, gap.right, left_scaffold, right_scaffold)
                            output_for_gfa.append(gfa_output)
                        #Add the bad solutions to a FASTA file containing all bad solutionss
                        else:
                            output_for_gfa = []         #Pbm if several solutions output
                            with open(bad_solutions_file, "a") as bad_file:
                                bad_file.write("\n>" + str(gap.left) +":"+ str(gap.right) + "_gf.s" + str(seed_size) + ".o" + str(min_overlap) + "_fwd _ len_" + str(len(seq)) + "_qual_" + str(quality_gapfilled_seq) +"\n")
                                bad_file.write(str(seq))

                qualified.seek(0)


    #If no solution found, set 'output_for_gfa' to empty list
    else:
        output_for_gfa = []

    #----------------------------------------------------
    # GFA output: case gap, no solution
    #----------------------------------------------------
    #Save the current G line into the variable 'output_for_gfa' only if this variable is empty 
    if len(output_for_gfa) == 0:
        output_for_gfa.append([str(current_gap)])


    #TODO: remove the flanking_contig.fasta files

    os.chdir(outDir)


    return union_summary, output_for_gfa


#----------------------------------------------------
# Gapfilling with MindTheGap
#----------------------------------------------------
try:
    #Open the input GFA file
    gfa = gfapy.Gfa.from_file(gfa_file)
    #Create the output GFA file
    out_gfa_file = str(gfa_name).split('.gfa')[0] + "_olc.gfa"

    #----------------------------------------------------
    # GFA output: case no gap
    #----------------------------------------------------
    #If no gap, rewrite all the lines into GFA output
    if len(gfa.gaps) == 0:
        with open(out_gfa_file, "w") as f:
            out_gfa = gfapy.Gfa()
            for line in gfa.lines:
                out_gfa.add_line(str(line))
            out_gfa.to_file(out_gfa_file)

    #----------------------------------------------------   
    # Fill the gaps
    #----------------------------------------------------
    #If gap, rewrite the H and S lines into GFA output
    if args.line is None:
        with open(out_gfa_file, "w") as f:
            out_gfa = gfapy.Gfa()
            out_gfa.add_line("H\tVN:Z:2.0")
            for line in gfa.segments:
                out_gfa.add_line(str(line))
            out_gfa.to_file(out_gfa_file)
        
    gaps = []
    gaps_label = []
    #If '-line' argument provided, start analysis from this line in GFA file input
    if args.line is not None:
        for _gap_ in gfa.gaps[(args.line - (len(gfa.segments)+2)):]:
            _gap_ = str(_gap_)
            gaps.append(_gap_)
    else:
        #Convert Gfapy gap line to a string to be able to use it with multiprocessing
        for _gap_ in gfa.gaps:
            _gap_ = str(_gap_)
            gaps.append(_gap_)

    p = Pool()

    with open("{}.union.sum".format(gfa_name), "w") as union_sum:
        legend = ["Gap_ID", "Left_scaffold", "Right_scaffold", "Gap_size", "Chunk_size", "Nb_barcodes", "Nb_reads"]
        union_sum.write('\t'.join(j for j in legend))

        for union_summary, output_for_gfa in p.map(gapfilling, gaps):
            #Write all union_summary (obtained for each gap) from 'gapfilling' into the 'union_sum' file
            union_sum.write("\n" + '\t'.join(str(i) for i in union_summary))

            #Output the 'output_for_gfa' results (obtained for each gap) from 'gapfilling' in the output GFA file
            print("\nCreating the output GFA file...")
            if len(output_for_gfa[0]) > 1:          #solution found for the current gap
                for output in output_for_gfa:
                    gapfill_file = update_gfa_with_solution(outDir, gfa_name, output, out_gfa_file)
                    success = True
            else:                                   #no solution found for the current gap
                out_gfa = gfapy.Gfa.from_file(out_gfa_file)
                out_gfa.add_line(output_for_gfa[0][0])
                out_gfa.to_file(out_gfa_file)
                success = False


        p.close()

    #Remove the raw files obtained from MindTheGap
    os.chdir(olcDir)


except Exception as e:
    print("\nException-")
    print(e)
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, fname, exc_tb.tb_lineno)
    sys.exit(1)


print("\nThe results from OLC pipeline are saved in " + outDir)
print("The results from OLC Gap-filling are saved in " + olcDir)
print("The statistics from OLC pipeline are saved in " + statsDir)
print("Summary of the union: " +gfa_name+".union.sum")
print("GFA output file: " + out_gfa_file)
if success == True:
    print("Corresponding file containing all gapfill sequences: " + gapfill_file + "\n")

#----------------------------------------------------
#Summary output
#----------------------------------------------------
gfa_output = gfapy.Gfa.from_file(outDir +"/"+ str(out_gfa_file))

#Total initials gaps
total_gaps = []
for g_line in gfa.gaps:
    gap_start = str(g_line.sid1) +"_"+ str(g_line.sid2) 
    total_gaps.append(gap_start)
nb_total_gaps = len(total_gaps)
print("------------------------------------------------------------------------------------------------------------------------\n")
print("Attempt to gap-fill {} gaps \n".format(nb_total_gaps))

#Gap(s) not gap-filled
no_gapfill = []
for g_line in gfa_output.gaps:
    gap_end = str(g_line.sid1) +"_"+ str(g_line.sid2) 
    no_gapfill.append(gap_end)
    print("The gap {} was not successfully gap-filled".format(gap_end))

nb_gapfill = len(total_gaps) - len(no_gapfill)
print("\nIn total, {} gaps were successfully gap-filled:\n".format(str(nb_gapfill)))


#Gaps gap-filled
out_fasta_file = outDir +"/"+ gapfill_file
gap_names = []
if (out_fasta_file) is not None:
    with open(out_fasta_file, "r") as gapfilled:
        for record in SeqIO.parse(gapfilled, "fasta"):
            gap_name = str(record.id).split('_')[0]

            #For a new gap
            if gap_name not in gap_names:
                gap_names.append(gap_name)
                k = str(record.id).split('.k')[-1].split('_')[0]
                print("\t* " + gap_name + "\tk" + k)

            #For all gaps
            orientation = str(record.id).split('_')[-1]
            length = str(record.description).split('_ len_')[1].split('_qual_')[0]
            quality = str(record.description).split('_qual_')[1]
            print("\t\t* " + orientation + "\t" + length + " bp\t" + quality)
           
print("\n")

#TODO: two modules, one when reference sequence provided (args.refDir), one when no reference sequence is provided (args.scaff)
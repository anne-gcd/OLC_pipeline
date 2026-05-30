#!/usr/bin/env python3
from __future__ import print_function
import os
import sys
import re
import csv
import argparse
import subprocess
import gfapy
from gfapy.sequence import rc
from Bio import SeqIO, Align
from Bio.Seq import Seq


#PairwiseAligner object
aligner = Align.PairwiseAligner()
aligner.match_score = 1
aligner.mismatch_score = -1
aligner.query_end_gap_score = 0
aligner.internal_open_gap_score = -1
aligner.internal_extend_gap_score = -0.5
aligner.target_end_open_gap_score = -1
aligner.target_end_extend_gap_score = -0.5


#----------------------------------------------------
# Arg parser
#----------------------------------------------------
parser = argparse.ArgumentParser(prog="stats_alignment.py", usage="%(prog)s -qry <query_sequences_file> -ref <reference_sequence> -ext <extension_size> -p <output_file_prefix> [options]", \
                                formatter_class=argparse.RawTextHelpFormatter, \
                                description=(''' \
                                Statistics about the inserted sequence obtained from MindTheGap (-qry)
                                Note: there are kmer flanking regions on the edges of the inserted sequence (which are included in '-ext' bp flanking regions)
                                '''))

parser.add_argument("-qry", "--query", action="store", help="File containing the inserted sequences obtained from MindTheGap (format: 'xxx.insertions.fasta')", required=True)
parser.add_argument("-ref", "--reference", action="store", help="File containing either the reference sequence or the flanking contigs sequences of the gap (format: 'xxx.fasta')", required=True)
parser.add_argument("-ext", "--ext", action="store", type=int, help="Extension size of the gap, on both sides; determine start/end of gapfilling", required=True)
parser.add_argument("-p", "--prefix", action="store", help="Prefix of output file to save the statistical results", required=True)
parser.add_argument("-out", "--outDir", action="store", default="./mtglink_results/alignments_stats", help="Output directory for saving results")

args = parser.parse_args()

if re.match('^.*.fasta$', args.query) is None:
    parser.error("Warning: Qualitative evaluation _ The suffix of the inserted sequences (query sequences) file should be: '.fasta'")

if re.match('^.*.fasta$', args.reference) is None:
    parser.error("Warning: Qualitative evaluation _ The suffix of the reference sequence file should be: '.fasta'")

#----------------------------------------------------
# Input files
#----------------------------------------------------
#Query file: inserted/gap-filled sequences file
qry_file = os.path.abspath(args.query)
if not os.path.exists(args.query):
    parser.error("Warning: Qualitative evaluation _ The path of the query file (inserted sequences file) doesn't exist")

#Reference file: containing either the reference sequence or the flanking contigs sequences
ref_file = os.path.abspath(args.reference)
if not os.path.exists(ref_file):
    parser.error("Warning: Qualitative evaluation _ The path of the reference file doesn't exist")

#----------------------------------------------------
# Directory for saving results
#----------------------------------------------------
cwd = os.getcwd()
if not os.path.exists(args.outDir):
    os.mkdir(args.outDir)
try:
    os.chdir(args.outDir)
except:
    print("Something wrong with specified directory. Exception-", sys.exc_info())
    print("Restoring the path")
    os.chdir(cwd)
outDir = os.getcwd()

try:
    #-----------------------------------------------------------------------------
    # Statistics about the Alignment Ref vs Qry
    #-----------------------------------------------------------------------------

    if not re.match('^.*.contigs.fasta$', args.reference):
        #----------------------------------------------------
        # Ref = reference sequence of simulated gap
        #----------------------------------------------------
        #Run NUCmer to obtain alignment of the reference sequence against the query's sequences
        prefix = args.prefix + ".ref_qry"

        log_file = str(prefix) + ".log"
        with open(log_file, "a") as log:
            log.write("Query file: " + str(qry_file) + "\n")
            log.write("Reference file" + str(ref_file) + "\n")
            log.write("The results are saved in " + outDir)

        nucmerLog = "{}_nucmer_ref_qry.log".format(args.prefix)
        delta_file = prefix + ".delta"
        coords_file = prefix + ".coords.unsorted" 

        nucmer_command = ["nucmer", "--maxmatch", "-p", prefix, ref_file, qry_file]
        coords_command = ["show-coords", "-rcdlT", delta_file]

        with open(coords_file, "w") as coords, open(nucmerLog, "a") as log:
            subprocess.run(nucmer_command, stderr=log)
            subprocess.run(coords_command, stdout=coords, stderr=log)

        #Sort the 'xxx.coords.unsorted' file for further analysis
        coords_sorted_file = prefix + ".coords"
        sort_command = ["sort", "-n", coords_file]
        with open(coords_sorted_file, "w") as coords_sorted:
            subprocess.run(sort_command, stdout=coords_sorted)

        #Output stats file of alignment query vs ref
        ref_qry_output = outDir + "/" + args.prefix + ".ref_qry.alignment.stats.unsorted"
        stats_legend = ["Gap", "Len_gap", "Chunk", "Seed_size", "Min_overlap", "Len_Q", "Ref", "Len_R", \
                        "Start_ref", "End_ref", "Start_qry", "End_qry", "Len_alignR", "Len_alignQ", "%_Id", "%_CovR", "%_CovQ", "Frame_R", "Frame_Q", "Quality"]
        
        #Get the gap size, chunk size, seed size and overlap min values
        gap_size = qry_file.split('.')[-6]
        g = int("".join(list(gap_size)[1:]))
        chunk_size = qry_file.split('.')[-5]
        c = int("".join(list(chunk_size)[1:]))
        seed_size = qry_file.split('.')[-4]
        s = int("".join(list(seed_size)[1:]))
        min_overlap = qry_file.split('.')[-3]
        o = int("".join(list(min_overlap)[1:]))
        qry_id = qry_file.split('/')[-1].split('.')[0]

        #Get output values from NUCmer:
        reader = csv.DictReader(open(coords_sorted_file), \
                                    fieldnames=("S1", 'E1', "S2", "E2", "LEN_1", "LEN_2", "%_IDY", "LEN_R", "LEN_Q", "COV_R", "COV_Q", "FRM_R", "FRM_Q", "TAG_1", "TAG_2"), \
                                    delimiter='\t')

        rows = list(reader)
        for row in rows[3:]:

            len_q = row["LEN_Q"]
            ref = row["TAG_1"]
            len_r = row["LEN_R"]
            start_r = row["S1"]
            end_r = row["E1"]
            start_q = row["S2"]
            end_q = row["E2"]
            len_align_r = row["LEN_1"]
            len_align_q = row["LEN_2"]
            identity = row["%_IDY"]
            cov_r = row["COV_R"]
            cov_q = row["COV_Q"]
            frame_r = row["FRM_R"]
            frame_q = row["FRM_Q"]

            # Estimate quality of gapfilled sequence
            ref_len = int(len_r)
            error_10_perc = int(0.1 * ref_len)
            qry_len = int(len_q) - 2*args.ext - 2*31

            #length of query sequence is equal +-10% of ref length
            if qry_len in range((ref_len - error_10_perc), (ref_len + error_10_perc)):
                #the gapfilled seq matches to the whole ref seq
                if int(len_align_q) == ref_len:
                    quality_rq = 'A'
                #the gapfilled seq matches to the ref seq +-10% of ref length
                elif int(len_align_q) in range((ref_len - error_10_perc), (ref_len + error_10_perc)):
                    quality_rq = 'B'
                #the gapfilled seq matches to the ref seq, but not along all their length (>= 50% of their length align)
                elif int(len_align_q) >= int(0.5*ref_len) and int(len_align_r) >= int(0.5*qry_len):
                    quality_rq = 'C'
                else:
                    quality_rq = 'D'

            else:
                quality_rq = 'D'

            #Write stats results in output file
            stats = [qry_id, g, c, s, o, len_q, ref, len_r, \
                    start_r, end_r, start_q, end_q, len_align_r, len_align_q, identity, cov_r, cov_q, frame_r, frame_q, quality_rq]

            if os.path.exists(ref_qry_output):
                with open(ref_qry_output, "a") as output:
                    output.write('\n' + '\t'.join(str(i) for i in stats))
            else:
                with open(ref_qry_output, "a") as output:
                    output.write('\t'.join(j for j in stats_legend))
                    output.write('\n'+'\n' + '\t'.join(str(i) for i in stats))

        #Sort the 'xxx.alignment.stats.unsorted' file for further analysis
        ref_qry_sorted = outDir + "/" + args.prefix + ".ref_qry.alignment.stats"
        order_command = ["sort", "-k6,7", "-k11,12n", "-r", ref_qry_output]
        with open(ref_qry_sorted, "w") as r_sorted:
            subprocess.run(order_command, stdout=r_sorted)


    elif re.match('^.*.contigs.fasta$', args.reference):
        #----------------------------------------------------
        # Ref = contigs' sequences
        #----------------------------------------------------
        #Run NUCmer to obtain alignment of extension portions (-ext) (of the flanking contigs) against the query's sequences
        prefix = args.prefix + ".ref_qry"

        log_file = str(prefix) + ".log"
        with open(log_file, "a") as log:
            log.write("Query file: " + str(qry_file) + "\n")
            log.write("Reference file" + str(ref_file) + "\n")
            log.write("The results are saved in " + outDir)

        nucmerLog = "{}_nucmer_ref_qry.log".format(args.prefix)
        delta_file = prefix + ".delta"
        coords_file = prefix + ".coords"

        nucmer_command = ["nucmer", "-p", prefix, ref_file, qry_file]
        coords_command = ["show-coords", "-rcdlT", delta_file]

        with open(coords_file, "w") as coords, open(nucmerLog, "a") as log:
            subprocess.run(nucmer_command, stderr=log)
            subprocess.run(coords_command, stdout=coords, stderr=log)

        #Output stats file of alignment query vs ref
        ref_qry_output = outDir + "/" + args.prefix + ".ref_qry.alignment.stats"
        stats_legend = ["Gap", "Len_gap", "Chunk", "Seed_size", "Min_overlap", "Len_Q", "Ref", "Len_R", \
                        "Start_ref", "End_ref", "Start_qry", "End_qry", "Len_alignR", "Len_alignQ", "%_Id", "%_CovR", "%_CovQ", "Frame_R", "Frame_Q", "Quality"]
        
        #Get the gap size, chunk size, seed size and overlap min values
        gap_size = qry_file.split('.')[-6]
        g = int("".join(list(gap_size)[1:]))
        chunk_size = qry_file.split('.')[-5]
        c = int("".join(list(chunk_size)[1:]))
        seed_size = qry_file.split('.')[-4]
        s = int("".join(list(seed_size)[1:]))
        min_overlap = qry_file.split('.')[-3]
        o = int("".join(list(min_overlap)[1:]))
        qry_id = qry_file.split('/')[-1].split('.')[0]

        #Get output values from NUCmer:
        reader = csv.DictReader(open(coords_file), \
                                    fieldnames=("S1", 'E1', "S2", "E2", "LEN_1", "LEN_2", "%_IDY", "LEN_R", "LEN_Q", "COV_R", "COV_Q", "FRM_R", "FRM_Q", "TAG_1", "TAG_2"), \
                                    delimiter='\t')

        rows = list(reader)
        for row in rows[3:]:
            if row["TAG_1"].split("_")[0] in str(qry_id):

                len_q = row["LEN_Q"]
                ref = row["TAG_1"].split("_")[0]
                len_r = row["LEN_R"]
                start_r = row["S1"]
                end_r = row["E1"]
                start_q = row["S2"]
                end_q = row["E2"]
                len_align_r = row["LEN_1"]
                len_align_q = row["LEN_2"]
                identity = row["%_IDY"]
                cov_r = row["COV_R"]
                cov_q = row["COV_Q"]
                frame_r = row["FRM_R"]
                frame_q = row["FRM_Q"]

                # Estimate quality of gapfilled sequence
                left = str(qry_id).split('_')[0]
                left_scaffold = left[:-1]
                right = str(qry_id).split('_')[1]
                right_scaffold = right[:-1]
                error_10_perc = int(0.1 * args.ext)
                error_50_perc = int(0.5 * args.ext)

                #ref = Left scaffold
                if ref == left_scaffold:
                    #extension of qry match perfectly as expected to ref
                    if ('+' in left and int(start_q) == 32 and int(end_q) == (args.ext + 31)) or ('-' in left and int(start_q) == (args.ext + 31) and int(end_q) == 32):
                        quality_rq = 'A'
                    #extension of qry almost match as expected to ref (+-10% of extension size) 
                    elif ('+' in left and int(start_q) in range(32, (32 + error_10_perc + 1)) and int(end_q) in range((args.ext + 31 - error_10_perc), (args.ext + 31 + error_10_perc + 1))) or ('-' in left and int(start_q) in range((args.ext + 31 - error_10_perc), (args.ext + 31 + error_10_perc + 1)) and int(end_q) in range(32, (32 + error_10_perc + 1))):
                        quality_rq = 'B'
                    #extension of qry almost match as expected to ref (+-50% of extension size) 
                    elif ('+' in left and int(start_q) in range(32, (32 + error_50_perc + 1)) and int(end_q) in range((args.ext + 31 - error_50_perc), (args.ext + 31 + error_50_perc + 1))) or ('-' in left and int(start_q) in range((args.ext + 31 - error_50_perc), (args.ext + 31 + error_50_perc + 1)) and int(end_q) in range(32, (32 + error_50_perc + 1))):
                        quality_rq = 'C'
                    else:
                        quality_rq = 'D'

                #ref = Right scaffold
                elif ref == right_scaffold:
                    #extension of qry match perfectly as expected to ref
                    if ('+' in right and int(start_q) == (int(len_q) - 31 - args.ext + 1) and int(end_q) == (int(len_q) - 31)) or ('-' in right and int(start_q) == (int(len_q) - 31) and int(end_q) == (int(len_q) - 31 - args.ext + 1)):
                        quality_rq = 'A'
                    #extension of qry almost match as expected to ref (+-10% of extension size) 
                    elif ('+' in right and int(start_q) in range((int(len_q) - 31 - args.ext + 1 - error_10_perc), (int(len_q) - 31 - args.ext + 1 + error_10_perc + 1)) and int(end_q) in range((int(len_q) - 31 - error_10_perc), (int(len_q) - 31 + 1))) or ('-' in right and int(start_q) in range((int(len_q) - 31 - error_10_perc), (int(len_q) - 31 + 1)) and int(end_q) in range((int(len_q) - 31 - args.ext + 1 - error_10_perc), (int(len_q) - 31 - args.ext + 1 + error_10_perc + 1))):
                        quality_rq = 'B'
                    #extension of qry almost match as expected to ref (+-50% of extension size) 
                    elif ('+' in right and int(start_q) in range((int(len_q) - 31 - args.ext + 1 - error_50_perc), (int(len_q) - 31 - args.ext + 1 + error_50_perc + 1)) and int(end_q) in range((int(len_q) - 31 - error_50_perc), (int(len_q) - 31 + 1))) or ('-' in right and int(start_q) in range((int(len_q) - 31 - error_50_perc), (int(len_q) - 31 + 1)) and int(end_q) in range((int(len_q) - 31 - args.ext + 1 - error_50_perc), (int(len_q) - 31 - args.ext + 1 + error_50_perc + 1))):
                        quality_rq = 'C'
                    else:
                        quality_rq = 'D'
                
                #Write stats results in output file
                stats = [qry_id, g, c, s, o, len_q, ref, len_r, \
                        start_r, end_r, start_q, end_q, len_align_r, len_align_q, identity, cov_r, cov_q, frame_r, frame_q, quality_rq]

                if os.path.exists(ref_qry_output):
                    with open(ref_qry_output, "a") as output:
                        output.write('\n' + '\t'.join(str(i) for i in stats))
                else:
                    with open(ref_qry_output, "a") as output:
                        output.write('\t'.join(j for j in stats_legend))
                        output.write('\n'+'\n' + '\t'.join(str(i) for i in stats))


    #Remove the raw file obtained from statistics ('.log', '.delta', '.coords', '.unsorted' files)
    subprocess.run(["rm", nucmerLog])
    subprocess.run(["rm", delta_file])
    subprocess.run(["rm", coords_file])
    if not re.match('^.*.contigs.fasta$', args.reference):
        subprocess.run(["rm", coords_sorted_file])  #only when refDir


except Exception as e:
    print("\nException-")
    print(e)
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, fname, exc_tb.tb_lineno)
    sys.exit(1)
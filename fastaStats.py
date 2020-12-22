#!/usr/bin/env python3

# Parse a fasta file, split into windows and compute stats

import sys
import argparse
import os
import shutil
import numpy as np

# local imports
from utils import fasta
from utils import windows
from utils import GC
from utils import kmers
from utils import dictionaries

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", type=str, help = "The input fasta", required=True)
parser.add_argument("-w", "--windowSize", type=int, help = "The size of windows to iterate over. Default 1000", default=1000)
parser.add_argument("-v", "--overlap", type=int, help = "The overlap of windows. Default is same as window size - i.e. no overlap", default=1000)
parser.add_argument("-k", "--kmerLength", type=int, help = "The length of kmer to calculate diversity of. Default is 4 - i.e. tetranucleotide diversity", default=4)
parser.add_argument("-o", "--outdir", type=str, help = "The output directory for GC stats", default="fastaStats_output")
args = parser.parse_args()

input_fasta = args.fasta
windowSize = args.windowSize
overlap = args.overlap
kmerLength = args.kmerLength
outdir = args.outdir

if os.path.exists(outdir):
		print("[WARNING] \tRemoving existing output directory ('" + outdir + "') and its contents.")
		shutil.rmtree(outdir)

os.mkdir(outdir)

# write a general info file
GeneralStats = open(outdir + '/' + 'GeneralStats.txt', 'w')
GeneralStats.write('Arguments used: ' + " ".join(sys.argv[1:]) + '\n')

# write to CSV
GCperWindowCSV = open(outdir + '/' + 'FastaStats.csv', 'w')
GCperWindowCSV.write('ID,bin,GCPercent,GCSkew,UniqueKmers\n')

# bin integer
bint = 0
# keep track of fastas being processed
fastaCount = 0
totalSeqLength = 0
seq_lengths = []
totalNucleotideCounts = {}

# go through each fasta
for header, sequence in fasta.parse_fastai(input_fasta):
    # general data
    fastaCount += 1
    seq_length = len(sequence)
    seq_lengths.append(seq_length)
    totalSeqLength += seq_length

    # begin sliding windows
    wins = windows.slidingWindow(sequence, size=windowSize, step=overlap, fillvalue="-")
    # let's count the windows
    i = 0
    for window in wins:
        if i % 100 == 0:
            print("[STATUS] \t" + str(i) + " windows processed for " + str(header) + ".", end = "\r")
        # get the sequence
        seq = ''.join(window)
        # calculate GC stats
        currentWindow = str(bint) + '-' + str(bint + windowSize)
        PerGC, GCSkew, nucleotideCounts = GC.GCStats(seq)
        # calculate kmer stats
        UniqueKmers = kmers.getUniqueKmers(seq, kmerLength)
        formattedHeader = header.replace(">.", "")
        GCperWindowCSV.write(formattedHeader + ',' + currentWindow + ',' + str(PerGC) + ',' + str(GCSkew) + ',' + str(UniqueKmers) + '\n')
        
        # keep running total of nucleotide counts
        totalNucleotideCounts = dictionaries.mergeDictionaries(totalNucleotideCounts, nucleotideCounts)

        if bint < seq_length - overlap:
            bint += overlap
        else:
            bint = 0
        i += 1
    # end sliding windows
    print("\n[STATUS] \t" + "contig " + str(fastaCount) + " processed.")

GCperWindowCSV.close()

# general stats about the genome
GeneralStats.write("Total number of contigs processed: " + str(fastaCount) + '\n')
GeneralStats.write("Total sequence length processed: " + str(totalSeqLength) + '\n')
# global GC
globalGC = (totalNucleotideCounts['G'] + totalNucleotideCounts['C']) / (totalNucleotideCounts['G'] + totalNucleotideCounts['C'] + totalNucleotideCounts['A'] + totalNucleotideCounts['T'])
GeneralStats.write("Global GC%: " + str(globalGC) + "\n")

# calculate L10-50 & N10-50
lengthsArray = np.array(seq_lengths)
sortedLengths = lengthsArray[np.argsort(-lengthsArray)]
cumulativeSums = np.cumsum(sortedLengths)
for level in [10, 20, 30, 40, 50]:
        nx = int(totalSeqLength * (level / 100))
        cSumN = min(cumulativeSums[cumulativeSums >= nx])
        l_level = int(np.where(cumulativeSums == cSumN)[0])
        n_level = int(sortedLengths[l_level])
        GeneralStats.write("L" + str(level) + ": " + str(l_level) + "\n")
        GeneralStats.write("N" + str(level) + ": " + str(n_level) + "\n")

GeneralStats.close()
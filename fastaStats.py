#!/usr/bin/env python3

# Parse a fasta file, split into windows and compute various metrics 

import sys
import argparse
import os
import shutil

#import GC as GC
#import fasta as fasta
#import windows as windows

from utils import fasta
from utils import windows
from utils import GC
from utils import kmers

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
		print("[WARNING] \tRemoving existing output directory ('" + outdir + "') and its contents...")
		shutil.rmtree(outdir)

os.mkdir(outdir)

fasta_seqs = fasta.parse_fasta(input_fasta, "")


GCperWindow = open(outdir + '/' + 'FastaStats.csv', 'w')
# write the CSV headers
GCperWindow.write('ID, bin, GC%, GCSkew, UniqueKmers\n')
# bin integer
bint = 0
# keep track of fastas being processed
fastaCount = 0
totalSeqLength = 0

for header, sequence in fasta_seqs.items():
    fastaCount += 1
    seq_length = len(sequence)
    totalSeqLength += seq_length
    wins = windows.slidingWindow(sequence, size=windowSize, step=overlap, fillvalue="-")

    # begin sliding windows
    for window in wins:
        seq = ''.join(window)
        currentWindow = str(bint) + '-' + str(bint + windowSize)
        PerGC, GCSkew = GC.GCStats(seq, currentWindow)
        UniqueKmers = kmers.getUniqueKmers(seq, kmerLength)
        formattedHeader = header.replace(">.", "")
        GCperWindow.write(header + ', ' + currentWindow + ', ' + str(PerGC) + ', ' + str(GCSkew) + ', ' + str(UniqueKmers) + '\n')
        
        if bint < seq_length - overlap:
            bint += overlap
        else:
            bint = 0
    # end sliding windows
    print("[STATUS] \t" + "contig " + str(fastaCount) + " processed.")

GCperWindow.close()

print("\nTotal number of contigs processed: " + str(fastaCount) + ".")
print("Total sequence length processed: " + str(totalSeqLength) + ".\n")
#!/usr/bin/env python3

# calculate tetranucleotide frequencies on a fasta file.
# adapted from https://github.com/CGATOxford/cgat/blob/cc3c5b332187b0f504d2906278a59fc0ccd2e504/CGAT/scripts/fasta2kmercontent.py

import itertools
import re
import sys

def getUniqueKmers(sequence, kmerLength):
    sequence = sequence.upper()
    if kmerLength > 7:
        sys.exit("kmer should be less than 8, this algorithm might be inefficient for large k.")

    nucleotides = []
    for nucleotide in ["A", "C", "T", "G"]:
        nucleotides = nucleotides + [x for x in itertools.repeat(nucleotide, kmerLength)]

    # get all kmer sequences to query
    kmers = set()
    for kmer in itertools.permutations(nucleotides, kmerLength):
        kmers.add(kmer)
    
    # store the results
    kmerCounts = {}
    for kmer in kmers:
            counts = [m.start()
                      for m in re.finditer("".join(kmer), sequence)]
            kmerCounts[kmer] = len(counts)
    # sum unique kmers. Information is chucked away here.
    uniqueKmers = 0
    for kmer_count in list(kmerCounts.values()):
        if kmer_count > 0:
            uniqueKmers += 1
    
    return uniqueKmers
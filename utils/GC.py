#!/usr/bin/env python3

# calculate the GC percentage on a strand of DNA.
# not adjusted for missing data..?
# Is GC skew important??

def GCStats(sequence):
    # make sure we use only upper case
    sequence = sequence.upper()
    # the four bases
    nucleotides = ['A', 'C', 'G', 'T'] 
    nucleotideCounts = {}

    # iterate over each base in the string
    for nucleotides in sequence:
        # count the bases
        count = sequence.count(nucleotides)
        nucleotideCounts[nucleotides] = count
        # catch missing key/value pairs
        if 'A' not in nucleotideCounts:
            nucleotideCounts['A'] = 0
        elif 'C' not in nucleotideCounts:
            nucleotideCounts['C'] = 0
        elif 'G' not in nucleotideCounts:
            nucleotideCounts['G'] = 0
        elif 'T' not in nucleotideCounts:
            nucleotideCounts['T'] = 0

    # GC percent in a sequence
    GC = nucleotideCounts['G'] + nucleotideCounts['C']
    sums = 0
    for key in nucleotideCounts:
        sums += nucleotideCounts[key]
    percentGC = (GC / sums) * 100

    # GC skew!
    numerator = nucleotideCounts['G'] - nucleotideCounts['C']
    denominator = nucleotideCounts['G'] + nucleotideCounts['C']
    if denominator == 0:
        GCskew = 0
    else:
        GCskew = numerator/denominator

    return percentGC, GCskew, nucleotideCounts
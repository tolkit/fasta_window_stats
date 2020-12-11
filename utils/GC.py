#!/usr/bin/env python3

# calculate the GC percentage on a strand of DNA.
# not adjusted for missing data..?
# Is GC skew important??

def GCStats(sequence, window):
    # make sure we use only upper case
    sequence = sequence.upper()
    # the four bases
    nucleotides = ['A', 'C', 'G', 'T'] 
    nucleotide_counts = {}

    # iterate over each base in the string
    for base in sequence:
        # catch missing key/value pairs
        if 'A' not in nucleotide_counts:
            print("[WARNING] \tbase A not present in window " + str(window))
            nucleotide_counts['A'] = 0
        elif 'C' not in nucleotide_counts:
            print("[WARNING] \tbase C not present in window " + str(window))
            nucleotide_counts['C'] = 0
        elif 'G' not in nucleotide_counts:
            print("[WARNING] \tbase G not present in window " + str(window))
            nucleotide_counts['G'] = 0
        elif 'T' not in nucleotide_counts:
            print("[WARNING] \tbase T not present in window " + str(window))
            nucleotide_counts['T'] = 0
        # and count the bases
        count = sequence.count(base)
        nucleotide_counts[base] = count

    # GC percent in a sequence
    GC = nucleotide_counts['G'] + nucleotide_counts['C']
    sums = 0
    for key in nucleotide_counts:
        sums += nucleotide_counts[key]
    percentGC = (GC / sums) * 100

    # GC skew!
    numerator = nucleotide_counts['G'] - nucleotide_counts['C']
    denominator = nucleotide_counts['G'] + nucleotide_counts['C']
    if denominator == 0:
        GCskew = 0
    else:
        GCskew = numerator/denominator

    return percentGC, GCskew
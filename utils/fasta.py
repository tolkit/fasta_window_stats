#!/usr/bin/env python3

from itertools import groupby

# parse a fasta file
# poached from https://github.com/MikeTrizna/assembly_stats/blob/master/assembly_stats/assembly_stats.py
def parse_fastai(fasta_file):
    fh = open(fasta_file)
    fa_iter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in fa_iter:
        # drop the ">"
        header = next(header)[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.upper().strip() for s in next(fa_iter))
        yield header, seq
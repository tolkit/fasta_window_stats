# Sequence statistics

## Main

A start on some basic statistics we might want to compute on fasta/q files. Implemented so far is GC%, GC skew, and unique kmer counts in sliding (and optionally overlapping) windows. It could be cool if people could look at the code and find ways to make it faster. Some bits have been poached from other sources, which have been indicated in the source files.

Total number of contigs processed and the total sequence length processed are also just printed out at the end.

```
usage: fastaStats.py [-h] -f FASTA [-w WINDOWSIZE] [-v OVERLAP] [-k KMERLENGTH] [-o OUTDIR]

optional arguments:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        The input fasta
  -w WINDOWSIZE, --windowSize WINDOWSIZE
                        The size of windows to iterate over. Default 1000
  -v OVERLAP, --overlap OVERLAP
                        The overlap of windows. Default is same as window size - i.e. no overlap
  -k KMERLENGTH, --kmerLength KMERLENGTH
                        The length of kmer to calculate diversity of. Default is 4 - i.e. unique tetranucleotide frequencies
  -o OUTDIR, --outdir OUTDIR
                        The output directory for fasta stats
```

Currently output is to a CSV in an output directory.

I've only tested this so far on the two test files in this repo.

## Dependencies

The main script uses itertools, sys, os, re, collections, argparse, shutil.

### Output 

If you clone this locally, the example.html should work... but I made a private example at https://observablehq.com/d/66426ebaca01559d. 
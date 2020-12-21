# Sequence statistics

## Main

A start on some basic statistics we might want to compute on fasta files. Implemented so far is GC%, GC skew, and unique kmer counts in sliding (and optionally overlapping) windows. It could be cool if people could look at the code and find ways to make it faster. Additional functionality hopefully soon, maybe some *actual* statistics in each window. Some bits have been poached from other sources, which have been indicated in the source files.

A general statistics file is also generated, containing number of contigs (chromosomes) processed, total sequence length, global GC%, and L/N10-50. e.g.

    Arguments used: -f Athaliana_genome/Athaliana_1_5_m_c.fasta
    Total number of contigs processed: 7
    Total sequence length processed: 119668634
    Global GC%: 0.3605598671007913
    L10: 0
    N10: 30427671
    L20: 0
    N20: 30427671
    L30: 1
    N30: 26975502
    L40: 1
    N40: 26975502
    L50: 2
    N50: 23459830

## Usage

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

Note that the script currently only works on unzipped fasta files.

## Speed

I re-ran the updated script on the genome of *Arabidopsis thaliana* again and all 5 chromosomes, plus the plastid/mitochondrion took:

```
real	15m47.850s
user	15m43.616s
sys	0m2.366s
```

A significant improvement to ~2 hours. 

### Output visualisation

If you clone this locally and fire up a server, the example.html should work... but I made a private example at https://observablehq.com/d/66426ebaca01559d. 
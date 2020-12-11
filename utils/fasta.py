#!/usr/bin/env python3

# parse a fasta file
# see https://github.com/lstevens17/busco2phylo-nf/blob/main/busco2fasta.py (thanks Lewis!)

def parse_fasta(fasta_file, prefix):
	with open(fasta_file) as fasta:
		fasta_dict = {}
		for line in fasta:
			if line.startswith(">"):
				header = ">" + prefix + "." + line.rstrip("\n").replace(">", "").replace(",", "")#.split(" ")[0]
				fasta_dict[header] = ''
			else:
				fasta_dict[header] += line.rstrip("\n")
	return fasta_dict
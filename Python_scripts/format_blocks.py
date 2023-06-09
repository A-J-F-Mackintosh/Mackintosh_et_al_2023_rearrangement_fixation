#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Usage: format_blocks.py -i <STR> -k <INT> -l <STR> -o <STR> [-h]

  [Options]
    -i, --input <STR>                          Tsv of blocks (seq, start, end, samples, [mutuple])
    -k, --kmax <INT>                           kmax
    -l, --lump <STR>                           Lump all blocks with any entry above kmax (choose from True or False)
    -o, --output <STR>                         Output prefix
    -h, --help                                 Show this message

"""

import sys
import gzip
from docopt import docopt
import collections
import re

def open_blocks(blocks):
	if blocks.endswith(".tsv"):
		with open(blocks, "r") as blocks_f:
			for line in blocks_f:
				yield line

	elif blocks.endswith(".tsv.gz"):
		with gzip.open(blocks, "rt") as blocks_f:
			for line in blocks_f:
				yield line

	else:
		sys.exit("[X] blocks file does not end in '.tsv' or '.tsv.gz'")

def update_blocks(blocks, kmax, lump):
	new_tally = []
	new_tally.append("seqs\tstart\tend\tsamples\tvariation")
	for line in open_blocks(blocks):
		line = line.rstrip()
		if line.startswith("seq"):
			pass
		else:
			seq, start, end, samples, bSFS = line.split("\t")

			new_bSFS_singletons = 0
			new_bSFS_doubletons = 0
			new_bSFS_tripletons = 0

			for i, mu in enumerate(map(int, re.findall(r'[0-9]+', bSFS))):
				if i == 0:
					new_bSFS_singletons += mu
				elif i == 1:
					new_bSFS_doubletons += mu
				elif i == 2:
					new_bSFS_tripletons += mu

			if lump == "False":
				if new_bSFS_singletons > kmax:
					new_bSFS_singletons = kmax + 1
				if new_bSFS_doubletons > kmax:
					new_bSFS_doubletons = kmax + 1
				if new_bSFS_tripletons > kmax:
					new_bSFS_tripletons = kmax + 1

			if lump == "True":
				if max(list([new_bSFS_singletons, new_bSFS_doubletons, new_bSFS_tripletons])) > kmax:
					new_bSFS_singletons = kmax + 1
					new_bSFS_doubletons = kmax + 1
					new_bSFS_tripletons = kmax + 1
					
			new_bSFS = "{},{},{}".format(new_bSFS_singletons, new_bSFS_doubletons, new_bSFS_tripletons)
			new_tally.append("{}\t{}\t{}\t{}\t{}".format(seq, start, end, samples, new_bSFS))

	return new_tally

def write_new_tally(new_tally, outprefix):
	with open(str(outprefix) + ".formatted_blocks.tsv", "w") as fout:
		fout.write("\n".join(new_tally))
		fout.write("\n")

if __name__ == '__main__':
	__version__ = '0.1'
	args = docopt(__doc__)

	print(args)

	if args["--lump"] == "True" or args["--lump"] == "False":

		new_tally = update_blocks(args["--input"], int(args["--kmax"]), args["--lump"])
		write_new_tally(new_tally, args["--output"])

	else:
		sys.exit("[X] --lump should be True or False")
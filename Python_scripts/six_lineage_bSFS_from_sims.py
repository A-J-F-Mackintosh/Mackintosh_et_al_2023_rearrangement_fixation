#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Usage: six_lineage_bSFS.py -v <STR> -c <STR> -q <STR> -0 <INT> -1 <INT> -o <STR> -b <INT> -m <INT> -s <STR> [-h]

  [Options]
    -v, --vcf <STR>                            VCF file (gzipped or not)
    -c, --bed <STR>                            Callable bed file
    -s, --samples <STR>                        Comma delimited list of RGs
    -q, --sequence <STR>                       Sequence to block
    -0, --start <INT>                          Start of interval
    -1, --end <INT>                            End of interval
    -o, --outprefix <STR>                      Outprefix prefix
    -b, --block_size <INT>                     Callable bases in block
    -m, --max_span <INT>                       Max span of block
    -h, --help                                 Show this message

"""

import sys
import gzip
from docopt import docopt
import collections
import itertools

def get_blocks(samples, sequence, bed, block_size, max_span, seq_start, seq_end):

	block_starts = {}
	block_ends = {}
	block_seqs = {}
	block_ID = 0

	current_start = max_span * -2
	current_bases = 0
	current_seq = None

	with open(bed, "r") as bed_f:
		for line in bed_f:
			i_seq, i_start, i_end, i_num_samples, i_sample_string = line.rstrip().split()
			i_start = int(i_start)
			i_end = int(i_end)
			i_samples = set(i_sample_string.split(","))  

			# if interval includes samples, is the right sequence and roughly the right positions
			if samples.issubset(i_samples) and i_seq == sequence and i_end >= seq_start and i_start <= seq_end: 

				# check if we can make a block by adding to what is leftover from before
				if i_seq == current_seq and current_bases > 0:
					block_starts, block_ends, block_seqs, block_ID, current_start, current_bases = finish_old_block(current_start, 
						i_start, i_end, i_seq, block_size, max_span, block_starts, block_ends, block_seqs, block_ID, current_bases)

				else: # make a fresh block
					current_start = i_start
					block_starts, block_ends, block_seqs, block_ID, current_start, current_bases = add_fresh_block(current_start, 
						i_end, i_seq, block_size, max_span, block_starts, block_ends, block_seqs, block_ID)

				current_seq = i_seq # update seq

	return block_starts, block_ends, block_seqs 

def finish_old_block(current_start, i_start, i_end, i_seq, block_size, max_span, block_starts, block_ends, 
	block_seqs, block_ID, current_bases):

	candidate_block_end = i_start + (block_size - current_bases)

	if candidate_block_end - current_start <= max_span: # gap between block wasn't too big
		if candidate_block_end <= i_end: # can use remaining bases in interval to make block
			block_starts[block_ID] = current_start
			block_ends[block_ID] = candidate_block_end
			block_seqs[block_ID] = i_seq
			block_ID += 1
			current_start = candidate_block_end
			current_bases = 0
			return add_fresh_block(current_start, i_end, i_seq, block_size, max_span, block_starts, block_ends, 
				block_seqs, block_ID)

		else: # cannot use remaining bases in interval to finish block, add to it instead
			current_bases += (i_end - i_start)
			return block_starts, block_ends, block_seqs, block_ID, current_start, current_bases

	else: # gap between block was too big, make a fresh one
		current_start = i_start
		return add_fresh_block(current_start, i_end, i_seq, block_size, max_span, block_starts, block_ends, 
			block_seqs, block_ID)

def add_fresh_block(current_start, i_end, i_seq, block_size, max_span, block_starts, block_ends, block_seqs, block_ID):
	current_bases = 0
	while current_start + block_size <= i_end: # keep making blocks until we hit the end of the interval
		block_starts[block_ID] = current_start
		block_ends[block_ID] = current_start + block_size
		block_seqs[block_ID] = i_seq
		block_ID += 1
		current_start += block_size

	current_bases = i_end - current_start # once we hit the end, take what we can towards the next block
	return block_starts, block_ends, block_seqs, block_ID, current_start, current_bases

# one current limitation is that once a block fails, we try to make a new one starting at the interval we are working with
# ideally, we should try and make a new one starting further back, e.g. 1 base forward from the previous start

# the other limitation is that the bed file and vcf need to be sorted in the same way

def get_variation(block_starts, block_ends, block_seqs, samples, sequence, vcf, seq_start, seq_end):

	block_mutuples = {} # make a mutuple dict and populate it with monomorphic blocks
	for block in block_starts:
		block_mutuples[block] = [0, 0, 0]

	sample_fields = [] # variables to help loop through variants and blocks
	block_position = 0
	total_blocks = len(block_starts)

	for line in open_vcf(vcf):
		line = line.rstrip()
		if line.startswith("#"):
			if line.startswith("#CHROM"): # connect samples to vcf columns
				sample_fields = get_sample_fields(line, samples)

		else:
			fields = line.split() # split vcf lines
			seq = fields[0]
			pos = int(fields[1])
			all_genotypes = []

			if seq == sequence and pos >= seq_start and pos <= seq_end: # only tally for target sequence


				for i in sample_fields: # for each sample of interest, get genotypes
					genotypes = get_genotypes(fields[i])
					all_genotypes.append(genotypes[0])
					all_genotypes.append(genotypes[1])

				if None not in all_genotypes: # if genotypes are all non-missing
					if len(set(all_genotypes)) == 2: # if variant is biallelic
						sorted_genotypes = sorted(all_genotypes)
						allele_0_count = sorted_genotypes.count(sorted_genotypes[0])
						allele_1_count = sorted_genotypes.count(sorted_genotypes[-1])
						allele_diff = abs(allele_0_count - allele_1_count) # allele_diff allows us to assign an iton

						for block in range(block_position, total_blocks): # find block that the variant contributes to
							if seq == block_seqs[block] and pos > block_starts[block] and \
							pos <= block_ends[block]:
								if allele_diff == 0:
									block_mutuples[block][2] += 1 # folded tripleton
								elif allele_diff == 2:
									block_mutuples[block][1] += 1 # folded doubleton
								elif allele_diff == 4:
									block_mutuples[block][0] += 1 # folded singleton
								block_position = block
								break # a variant can only contribute to one block, so, once we have found it, break

	return block_starts, block_ends, block_seqs, block_mutuples

def write_output(block_starts, block_ends, block_seqs, block_mutuples, samples, outprefix, seq_start, seq_end):
	output_list = []
	for block in block_starts:
		if seq_start <= block_starts[block] and seq_end >= block_ends[block]:
			sample_list = list(samples)
			sample_file_string = ".{}_{}_{}".format(sample_list[0], sample_list[1], sample_list[2])
			sample_string = "{},{},{}".format(sample_list[0], sample_list[1], sample_list[2])
			block_string = "{}\t{}\t{}\t{}\t{}".format(block_seqs[block], block_starts[block], 
				block_ends[block], sample_string, block_mutuples[block]) 
			output_list.append(block_string)
	with open(outprefix + sample_file_string + ".blocks.tsv", "w") as fout:
		fout.write("\n".join(output_list))
		fout.write("\n")

def open_vcf(vcf):
	if vcf.endswith(".vcf"):
		with open(vcf, "r") as vcf_f:
			for line in vcf_f:
				yield line
	elif vcf.endswith(".vcf.gz"):
		with gzip.open(vcf, "rt") as vcf_f:
			for line in vcf_f:
				yield line
	else:
		sys.exit("[X] VCF file does not end in '.vcf' or '.vcf.gz'")

def get_sample_fields(line, samples):
	sample_fields = []
	for i, entry in enumerate(line.split()):
		if entry in samples:
			sample_fields.append(i)
	return sample_fields

def get_genotypes(genotype_field):
	genotypes = genotype_field
	if genotypes == "./." or genotypes == ".":
		return [None, None]
	else:
		return [genotypes[0], genotypes[2]]

if __name__ == '__main__':
	__version__ = '0.1'
	args = docopt(__doc__)
	print(args)
	
	block_singletons = collections.defaultdict(int)
	block_doubletons = collections.defaultdict(int)
	block_stripletons = collections.defaultdict(int)

	vcf = args["--vcf"]
	bed = args["--bed"]
	samples = set(args["--samples"].split(","))
	sequence = args["--sequence"]
	seq_start = int(args["--start"])
	seq_end = int(args["--end"])
	outprefix = str(args["--outprefix"])
	block_size = int(args["--block_size"])
	max_span = int(args["--max_span"])

	for sample_combo in itertools.combinations(samples, 3):

		sample_combo_set = set(sample_combo)
			
		print("[+] Making blocks for {} ...".format(sample_combo))
		block_starts, block_ends, block_seqs  = get_blocks(sample_combo_set, sequence, bed, block_size, max_span, seq_start, seq_end)

		print("[+] Tallying variants for {} ...".format(sample_combo))
		block_starts, block_ends, block_seqs, block_mutuples = get_variation(block_starts, block_ends, 
			block_seqs, sample_combo_set, sequence, vcf, seq_start, seq_end)

		print("[+] Writing output for {} ...".format(sample_combo))
		write_output(block_starts, block_ends, block_seqs, block_mutuples, sample_combo_set, outprefix, 
			seq_start, seq_end)

	print("[=] Done")
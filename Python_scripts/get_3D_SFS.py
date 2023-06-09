#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Usage: apeVCF.py -v <STR> -s <STR> -o <STR> -c <INT> [-r <BOL> -h]

  [Options]
    -v, --vcf <STR>                            VCF file (gzipped or not)
    -s, --samples <STR>                        Tab delimited samples file containing ID\tPOP, where ID matches readgroups in the vcf and POP is P1/P2/P3/O
    -o, --outprefix <STR>                      Outprefix SFS
    -r, --ref_polarisation <BOL>               Use the ref allele as the ancestral allele. This overules any O samples. Choose from True or False [default: False]
    -c, --callable <INT>                       Number of callable sites
    -h, --help                                 Show this message

"""

import sys
import gzip
from docopt import docopt
import collections
from fractions import Fraction

def generate_pop_dict(samples):
	# generate a dictionary where the key is the population and RGs are the values
	pops = set()
	pop_dict = collections.defaultdict(set)
	with open(samples, "r") as sample_f:
		for line in sample_f:
			ID = line.rstrip().split()[0]
			pop = line.rstrip().split()[1]
			pop_dict[pop].add(ID)
			pops.add(pop)
	if pops != {"P1", "P2", "P3"}:
		sys.exit("[X] Samples file does not contain P1, P2, P3")
	return pop_dict

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

def update_pop_dict(chrom_line, pop_dict):
	new_pop_dict = {}
	for i, ID in enumerate(chrom_line.split()):
		for pop in pop_dict:
			if ID in pop_dict[pop]:
				new_pop_dict[str(i)] = pop
	return new_pop_dict

def get_genotypes(genotype_field):
	genotypes = genotype_field.split(":")[0]
	genotypes = {genotypes[0], genotypes[2]}
	return genotypes

def fraction_to_float(fraction):
	a_float = fraction.numerator / fraction.denominator
	return round(a_float, 10)

def parse_vcf(vcf, ref_polarisation, pop_dict):
	pop_dict_updated = False # use this to only set the CHROM info once
	# counters
	unpolarised = 0
	multiallelic = 0
	SFS_dict = collections.defaultdict(int) # collect SFS counts here
	for line in open_vcf(vcf):
		line = line.rstrip()
		if line.startswith("#"):
			if line.startswith("#CHROM") and not pop_dict_updated: # connect samples to vcf columns
				pop_dict = update_pop_dict(line, pop_dict)
				pop_dict_updated = True
		else:
			fields = line.split()
			O_genotypes = set() # outgroup genotypes
			I_genotypes = set() # ingroup genotypes
			for i in range(9, len(fields)):
				genotypes = get_genotypes(fields[i])
				if pop_dict[str(i)] == "O":
					O_genotypes = O_genotypes.union(genotypes)
				else:
					I_genotypes = I_genotypes.union(genotypes)
			# multiallelic
			if len(I_genotypes.union(O_genotypes)) > 2:
				multiallelic += 1
			# if outgroup polymorphic
			elif len(O_genotypes) > 1:
				if ref_polarisation == "True": # not a problem if using the ref allele
					P1_der, P2_der, P3_der = count_alleles({"0"}, fields, pop_dict)
					SFS_dict[(P1_der, P2_der, P3_der)] += 1
				else:
					unpolarised += 1
			# easy ones
			elif ref_polarisation == "True":
				P1_der, P2_der, P3_der = count_alleles({"0"}, fields, pop_dict)
				SFS_dict[(P1_der, P2_der, P3_der)] += 1
			else:
				P1_der, P2_der, P3_der = count_alleles(O_genotypes, fields, pop_dict)
				SFS_dict[(P1_der, P2_der, P3_der)] += 1
	return SFS_dict, multiallelic, unpolarised

def count_alleles(ancestral_allele, fields, pop_dict):
	# this function count derived alleles in each population				
	P1_der = 0
	P2_der = 0
	P3_der = 0
	for i in range(9, len(fields)):
		genotypes = get_genotypes(fields[i])
		if pop_dict[str(i)] == "P1":
			if len(genotypes.intersection(ancestral_allele)) == 1:
				if len(genotypes) == 2:
					P1_der += 1
			elif len(genotypes.intersection(ancestral_allele)) == 0:
				P1_der += 2
		elif pop_dict[str(i)] == "P2":
			if len(genotypes.intersection(ancestral_allele)) == 1:
				if len(genotypes) == 2:
					P2_der += 1
			elif len(genotypes.intersection(ancestral_allele)) == 0:
				P2_der += 2
		elif pop_dict[str(i)] == "P3":
			if len(genotypes.intersection(ancestral_allele)) == 1:
				if len(genotypes) == 2:
					P3_der += 1
			elif len(genotypes.intersection(ancestral_allele)) == 0:
				P3_der += 2
	return P1_der, P2_der, P3_der

def correct_count(Y, X):
	# this function corrects SFS counts for double hit mutations
	# currently the error rate is hardcoded
	e_rate = 0.0212
	y = (Y - (Y*e_rate) - (X*e_rate)) / (1 - (2*e_rate))
	x = (X - (y*e_rate)) / (1 - e_rate)
	if y <= 0:
		return 0
	if x <= 0:
		y += x
	return round(y)

def write_counts(SFS, pop_dict, multiallelic, unpolarised, callable_sites, outprefix):
	P1_allels = 2*len(pop_dict["P1"]) # number of total allele in P1
	P2_allels = 2*len(pop_dict["P2"]) # number of total allele in P2
	P3_allels = 2*len(pop_dict["P3"]) # number of total allele in P3
	polymorphic_sites = sum(SFS_dict.values()) - SFS_dict[(0, 0, 0)] - SFS_dict[(P1_allels, P2_allels, P3_allels)]
	monomorphic_sites = callable_sites - polymorphic_sites - multiallelic - unpolarised
	SFS[(0, 0, 0)] = monomorphic_sites
	SFS_dict[(P1_allels, P2_allels, P3_allels)] = 0
	formatted_SFS = []
	for i in range(0, P1_allels + 1):
		for j in range(0, P2_allels + 1):
			for k in range(0, P3_allels + 1):
				"""
				# here we do a correction for double hits, which is probably a bad idea
				corrected_count = correct_count(SFS[(i, j, k)], SFS[((P1_allels-i), (P2_allels-j), (P3_allels-k))])
				formatted_SFS.append("{}\t{}".format((i, j, k), corrected_count))
				"""
				formatted_SFS.append("{}\t{}".format((i, j, k), SFS[(i, j, k)]))
	with open(str(outprefix) + ".3D_SFS.txt", "w") as fout:
		fout.write("\n".join(formatted_SFS))
		fout.write("\n")
	print("Done")

if __name__ == '__main__':
	__version__ = '0.1'
	args = docopt(__doc__)
	#print(args)
	pop_dict = generate_pop_dict(args['--samples'])
	SFS_dict, multiallelic, unpolarised = parse_vcf(args['--vcf'], args["--ref_polarisation"], pop_dict)
	write_counts(SFS_dict, pop_dict, multiallelic, unpolarised, int(args["--callable"]), args["--outprefix"])
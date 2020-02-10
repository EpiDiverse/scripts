#!/usr/bin/env python

'''
Title: filter_bedgraph.py
Date: 20190803
Author: Adam Nunn
Description:
	This program reads a bedGraph file from std.in to filter based on context and
	coverage values specified by the user, printing output to std.out

List of functions:
	main()
	build_genome()
	iterator()

Procedure:
	1. Iterate through reference FASTA to build genome dictionary
	2. Iterate through input BEDGRAPH to retrieve data for OUTPUT formatting

Usage:
	<bedGraph> | ./filter_bedgraph.py [options] <fasta> 
eg. tail -q -n+2 *.bedGraph | ./filter_bedgraph.py --coverage 2 --CAA --CG in.fa
'''

import sys, string
import argparse

####### Function main
def main(options):

	# 1) Iterate through reference FASTA to build genome dictionary
	genome = build_genome(options.fasta)

	fcon = list()
	# Obtain list of contexts
	if options.CAA or options.CHH: fcon.append("CAA") 
	if options.CAC or options.CHH: fcon.append("CAC")
	if options.CAT or options.CHH: fcon.append("CAT")
	if options.CCA or options.CHH: fcon.append("CCA")
	if options.CCC or options.CHH: fcon.append("CCC")
	if options.CCT or options.CHH: fcon.append("CCT")
	if options.CTA or options.CHH: fcon.append("CTA")
	if options.CTC or options.CHH: fcon.append("CTC")
	if options.CTT or options.CHH: fcon.append("CTT")
	if options.CAG or options.CHG: fcon.append("CAG")
	if options.CCG or options.CHG: fcon.append("CCG")
	if options.CTG or options.CHG: fcon.append("CTG")
	if options.CpG: fcon.append("CG")	

	# 2) Iterate through sys.stdin to generate output bedGraph
	iterator(genome,options.coverage,fcon,options.tabix,options.header)


####### Function to build 'genome' dictionary from 'FASTA' file object
def build_genome(FASTA):

	############################################################################################
	## FASTA		= 	string of path to reference genome eg. "/path/to/reference.fa"        ##
	############################################################################################

	# declare an empty dictionary and stage the 'FASTA' file
	genome = {}
	with open(FASTA, 'r') as fasta:
		first = True

		# iterate through the lines of the 'FASTA' file
		for line in fasta:
			line = line.rstrip()

			# identify the first sequence
			if line.startswith(">") and (first == True):
				line = line.split(" ")
				chr = line[0][1:]
				sequence = ''
				first = False

			# identify all following sequences
			elif line.startswith(">") and (first == False):
				genome[chr] = sequence.upper()
				line = line.split(" ")
				chr = line[0][1:]
				sequence = ''

			# append the line to the concatenated sequence
			else: sequence += line

		# append the final sequence to the dictionary according to chr
		genome[chr] = sequence.upper()

	# return the constructed genome dictionary
	print("Genome indexing complete\n#######", file=sys.stderr)
	return genome


####### Function to iterate std.in and extract subcontexts from 'genome' dictionary
def iterator(genome,coverage=0,fcon=None,tabix=False,header=False):

	############################################################################################
	## genome		= 	dictionary object from build_genome()                                 ##
	## coverage 	= 	integer eg. 4                                                         ##
	## fcon		 	= 	list of contexts eg. ["CG", "CAA", ...]                               ##
	## tabix	 	= 	boolean decision eg. True                                             ##
	## header	 	= 	boolean decision eg. True                                             ##
	############################################################################################

	if header:
		description = " ".join(fcon)
		print("track type=\"bedGraph\" description=\"combined contexts: {}\"".format(description))

	# 3) Iterate through input stdin to retrieve data for OUTPUT formatting
	for line in sys.stdin:
		line = line.rstrip()
		column = line.split("\t")

		# Filter coverage
		cov = int(column[4]) + int(column[5])
		if cov <= coverage: continue

		# determine strand, retrieve sequence
		CorG = genome[column[0]][int(column[1]):int(column[2])]
		
		if CorG == "C":
			strand = "+"
			seq = genome[column[0]][int(column[1]):int(column[2])+2]
			tri = seq
		
		elif CorG == "G":
			strand = "-"
			seq = genome[column[0]][int(column[1])-2:int(column[2])]
			tri = seq.translate(str.maketrans("ACGT","TGCA"))[::-1]

		# skip erroneous
		else: continue

		# Determine methylation context
		if len(tri) == 3: 
			if tri[1] == "G": context = "CG"
			else:
				if tri[2] == "G": context = "CHG"
				else: context = "CHH"
		elif len(tri) == 2:
			if tri[1] == "G": context = "CG"
			else: continue
		else: continue

		# Do not filter subcontext
		if len(fcon) == 0: 
			# Print tabix format
			if tabix: print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(column[0],column[1],strand,column[4],column[5],context,tri))
			# Print default format
			else: print(line)

		# Filter subcontext
		else:
			if (context == "CG" and "CG" in fcon) or (context != "CG" and tri in fcon):
				# Print tabix format
				if tabix: print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(column[0],column[1],strand,column[4],column[5],context,tri))
				# Print default format
				else: print(line)


# run main
if __name__ == '__main__':

	# define argparse
	usage = ''' This program reads a bedGraph file from std.in to filter based on context and
		coverage values specified by the user, printing output to std.out '''
	
	parser = argparse.ArgumentParser(description=usage)
	parser.add_argument('fasta', metavar='<fasta>',
						help='[REQUIRED] Path to input FASTA file')
	parser.add_argument('--CAA', required=False, default=False, action='store_true',
						help='Filter CAA context')
	parser.add_argument('--CAC', required=False, default=False, action='store_true',
						help='Filter CAC context')
	parser.add_argument('--CAT', required=False, default=False, action='store_true',
						help='Filter CAT context')
	parser.add_argument('--CCA', required=False, default=False, action='store_true',
						help='Filter CCA context')
	parser.add_argument('--CCC', required=False, default=False, action='store_true',
						help='Filter CCC context')
	parser.add_argument('--CCT', required=False, default=False, action='store_true',
						help='Filter CCT context')
	parser.add_argument('--CTA', required=False, default=False, action='store_true',
						help='Filter CTA context')
	parser.add_argument('--CTC', required=False, default=False, action='store_true',
						help='Filter CTC context')
	parser.add_argument('--CTT', required=False, default=False, action='store_true',
						help='Filter CTT context')
	parser.add_argument('--CAG', required=False, default=False, action='store_true',
						help='Filter CAG context')
	parser.add_argument('--CCG', required=False, default=False, action='store_true',
						help='Filter CCG context')
	parser.add_argument('--CTG', required=False, default=False, action='store_true',
						help='Filter CTG context')
	parser.add_argument('--CHH', required=False, default=False, action='store_true',
						help='Filter CHH context')
	parser.add_argument('--CHG', required=False, default=False, action='store_true',
						help='Filter CHG context')
	parser.add_argument('--CpG', required=False, default=False, action='store_true',
						help='Filter CG context')
	parser.add_argument('--tabix', required=False, default=False, action='store_true',
						help='Print in tabix format [default: off]')
	parser.add_argument('--header', required=False, default=False, action='store_true',
						help='Print a header line for integration with DMR pipeline [default: off]')
	parser.add_argument('--coverage', type=int, required=False, default=0,
						help='Minimum coverage threshold [default: 0]')
	parsed, unparsed = parser.parse_known_args()
	main(parsed)
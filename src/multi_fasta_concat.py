#!/usr/bin/env python3.6

'''
I take 2 or more multi-fasta files and I concatenate them together according to the sequence order appearing in the files.
First come, first served.
I replace missing nucleotides with "-".
I print out concatenated multifasta.

Andrea Del Cortona
2022/10/01
'''



#------------------------------------------------------------------#
# LOAD LIBRARIES

import argparse
import glob
import re
import sys
from datetime import datetime
from itertools import accumulate



#------------------------------------------------------------------#
# INPUT PARSER

parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter,
	description = '''

	==================================================================
	I take 2 or more multi-fasta files and I concatenate them together according to the sequence order appearing in the files.
	First come, first served.
	I replace missing nucleotides with "-".
	I print out concatenated multifasta.

	__________________________________________________________________

	Usage:
	python3.8 multi_fasta_concat.py --infiles "fasta1.fa,fasta2.fa,fasta3.fa..." 

	==================================================================
	''',
	epilog = '''
	__________________________________________________________________
			 Andrea Del Cortona - andrea.delcortona@gmail.com
						   2022-10-01
	__________________________________________________________________
	''')

parser.add_argument("--infiles", metavar = "INFILES", action = "store",
	type = str, dest = "INFILES", required = True,
	help = 'List of multi-fasta to be concatenated, comma separated.')

args = parser.parse_args()



#------------------------------------------------------------------#
# FUNCTIONS

# main function
def main():
	'''
	I take 2 or more multi-fasta files and I concatenate them together according to the sequence order appearing in the files.
	First come, first served.
	I replace missing nucleotides with "-".
	I print out concatenated multifasta.
	'''

	# create file lists
	*file_list, = args.INFILES.split(",")

	# create sequences database
	seq_DB = {
		"seq_list"   : [],       # list of unqiue fasta sequence identified
		"seq_length" : [],       # length of alignments in residues [156788, 2222, 3435]
		"occupancy"  : {},       # [0, 1] presence of sequences for each multifasta alignment
		"sequences"  : {}        # placeholder for multifasta sequences
	}

	# iterate fasta files and import them
	for file in file_list:
		with open(file) as infile:
			for line in infile:

				# remove newline
				line = line.rstrip('\n')

				# check if it header
				if line[0] == ">":
					# add sequence length

					# get fasta header
					seq_name = line[1:]
					last_seq = seq_name
					# is sequence already in list?
					if seq_name not in seq_DB["seq_list"]:
						seq_DB["seq_list"].append(seq_name)
						# create a placeholder in the "sequences" database
						seq_DB["sequences"][seq_name] = ""
						seq_DB["occupancy"][seq_name] = []
						# populate with missing residues where necessary and occupancy
						if len(seq_DB["seq_length"]) >= 1:
							for k in range(len(seq_DB["seq_length"])):
								seq_DB["sequences"][seq_name] += "-" * seq_DB["seq_length"][k]
								seq_DB["occupancy"][seq_name].append("0")
						seq_DB["occupancy"][seq_name].append("1")
					
					else:
						# update occupancy
						seq_DB["occupancy"][seq_name].append("1")
				
				else:
					# alignment sequence
					seq_DB["sequences"][last_seq] += line

			# add sequence length last alignment
			if seq_DB["seq_length"] == []:
				seq_DB["seq_length"].append(len(seq_DB["sequences"][last_seq]))
			else:
				seq_DB["seq_length"].append(len(seq_DB["sequences"][last_seq]) - list(accumulate(seq_DB["seq_length"]))[-1])
			
			# iterate sequences in multifasta and add "-" at the end of the sequences if necessary
			longest_sequence = max(len(item) for item in seq_DB["occupancy"].values())
			for sequence in seq_DB["seq_list"]:
				if len(seq_DB["occupancy"][sequence]) < longest_sequence:
					seq_DB["occupancy"][sequence].append("0")
					seq_DB["sequences"][sequence] += "-" * seq_DB["seq_length"][-1]

	# print matrix of sequence occupancies to standard error
	for key, value in seq_DB["occupancy"].items():
		sys.stderr.write("\t".join([key, "\t".join(value), "\n"]))

	# print concatenated multifasta alignment to standard output
	for sequence in sorted(seq_DB["seq_list"]):
		print(">" + sequence)
		print(seq_DB["sequences"][sequence])



#------------------------------------------------------------------#
# RUN

# run script and give the running time
if __name__ == '__main__':
	t0 = datetime.now()
	main()
	dt = datetime.now() - t0
	sys.stderr.write("# Time elapsed: %s\n" % dt)

#!/usr/bin/env python3.6

"""
I take a list of 


Andrea Del Cortona
2022/12/23
"""



#------------------------------------------------------------------#
# LOAD LIBRARIES

import argparse
import itertools
import os
import subprocess
import sys
from datetime import datetime



#------------------------------------------------------------------#
# INPUT PARSER

parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter,
	description = """
	==================================================================
	I 
	__________________________________________________________________
	Usage:
	python3.8 combinatorial_ML_phylogeny.py --infolder "fasta1.fa,fasta2.fa,..." \
		--names --iqtree --n_genes --outdir
	==================================================================
	""",
	epilog = """
	__________________________________________________________________
			 Andrea Del Cortona - andrea.delcortona@gmail.com
						   2022-12-23
	__________________________________________________________________
	""")

parser.add_argument("--infolder", metavar = "INFOLDER", action = "store",
	type = str, dest = "INFOLDER", required = True,
	help = "List of folders with input genes, comma separated.")

parser.add_argument("--names", metavar = "NAMES", action = "store",
	type = str, dest = "NAMES", required = True,
	help = "List of names of genes groups, comma separated.")

parser.add_argument("--iqtree", metavar = "IQTREE", action = "store",
	type = str, dest = "IQTREE", required = True,
	help = "iqtree binary to run the ML trees.")

parser.add_argument("--n_genes", metavar = "N_GENES", action = "store",
	type = str, dest = "N_GENES", required = True,
	help = "Number of max genes to concatenate for the phylogeny [min = 2].")

parser.add_argument("--outdir", metavar = "OUTDIR", action = "store",
	type = str, dest = "OUTDIR", required = True,
	help = "Output directory.")

args = parser.parse_args()



#------------------------------------------------------------------#
# FUNCTIONS

# main function
def main():
	"""
	I .
	"""

	# create dir lists
	*dir_list, = args.INFOLDER.split(",")

	# create name lists
	*name_list, = args.NAMES.split(",")

	# create empty input database
	gene_lists = {}

	# structure input database
	for i in range(len(dir_list)):
		gene_lists[name_list[i]] =  {
			"DIR"        : dir_list[i],		# what I do
			"genes"      : list(),			# who I am
			"alignments" : dict()			# I think I can remove
		}

	# reads genes
	for DIR in gene_lists:
		gene_lists[DIR]["genes"] = os.listdir(gene_lists[DIR]["DIR"])
		# add path to file
		for k in range(len(gene_lists[DIR]["genes"])):
			gene_lists[DIR]["genes"][k] = "/".join([gene_lists[DIR]["DIR"], gene_lists[DIR]["genes"][k]])

	# cp only combinations
	run_single_organelles_combination(
		OUTDIR = args.OUTDIR,
		N_GENES = args.N_GENES,
		IQTREE = args.IQTREE,
		gene_lists = gene_lists,
		outfolder_name1 = "01_cp_combinatorial",
		outfolder_name2 = "cp_conc_",
		dir_name = "CP"
	)

	# mt only combinations
	run_single_organelles_combination(
		OUTDIR = args.OUTDIR,
		N_GENES = args.N_GENES,
		IQTREE = args.IQTREE,
		gene_lists = gene_lists,
		outfolder_name1 = "01_mt_combinatorial",
		outfolder_name2 = "mt_conc_",
		dir_name = "MT"
	)

	# cp + mt combinations
	run_both_organelles_combination(
		OUTDIR = args.OUTDIR,
		N_GENES = args.N_GENES,
		IQTREE = args.IQTREE,
		gene_lists = gene_lists,
		outfolder_name1 = "01_cp_mt_combinatorial",
		outfolder_name2 = "cp_mt_conc_"
	)
	


# concatenate and run IQtree on single organelles
def run_single_organelles_combination(
		OUTDIR,
		N_GENES,
		IQTREE,
		gene_lists,
		outfolder_name1,
		outfolder_name2,
		dir_name
	):

	"""
	I taB.

		---
	OUTDIR : str
		output directory
	N_GENES : int
		max number of genes to concatenate
	IQTREE : str
		path to iqtree binary
	gene_lists : dict
		database with list of genes to concatenate
	outfolder_name1 : str
		name of the output folder
	outfolder_name2 : str
		name of the output folder
	dir_name : str
		name of the directory
	"""	

	# generate output folder
	os.makedirs("/".join([OUTDIR, outfolder_name1]), exist_ok = True)

	# get max number of gene
	max_n = len(gene_lists[dir_name]["genes"]) + 1 if len(gene_lists[dir_name]["genes"]) + 1 < int(N_GENES) + 1 else int(N_GENES) + 1
	# generate genes combinations
	combinations_DB = {}
	for k in range(2, max_n):
		
		current_dir = "/".join([OUTDIR, outfolder_name1, outfolder_name2 + str(k) + "_genes/"])
		os.makedirs(current_dir, exist_ok = True)
		
		combinations_DB[k] = list(itertools.combinations(gene_lists[dir_name]["genes"], k))

		for combination in combinations_DB[k]:
			
			# create sequences database
			seq_DB = {
				"seq_list"   : [],       # list of unique fasta sequence identified
				"seq_length" : [],       # length of alignments in residues [156788, 2222, 3435]
				"occupancy"  : {},       # [0, 1] presence of sequences for each multifasta alignment
				"sequences"  : {}        # placeholder for multifasta sequences
			}

			# get files and gene names
			file_list = combination
			gene_list = [gene.replace(".aln.fa", "") for gene in file_list]
			gene_list = [gene.replace(gene_lists[dir_name]["DIR"] + "/", "") for gene in gene_list]
			
			# import files and process them
			seq_DB = concatenate_multifasta(seq_DB, file_list)

			# print matrix of sequence occupancies
			outfile = open("/".join([current_dir, "_".join(gene_list) + ".concat.log"]), "w+")
			outfile.write("\t".join(["# Sequences", "\t".join(file_list), "\n"]))
			for key, value in seq_DB["occupancy"].items():
				outfile.write("\t".join([key, "\t".join(value), "\n"]))
			outfile.close()

			# print oncatenated multifasta alignment
			outfile = open("/".join([current_dir, "_".join(gene_list) + ".aln.fa"]), "w+")
			for sequence in sorted(seq_DB["seq_list"]):
				outfile.write(">" + sequence + "\n")
				outfile.write(seq_DB["sequences"][sequence] + "\n")
			outfile.close()

			# run IQtree
			bash_command = IQTREE + " -s " + "/".join([current_dir, "_".join(gene_list) + ".aln.fa"]) \
				+ " -st DNA -pre " + "/".join([current_dir, "_".join(gene_list) + ".aln"]) \
				+  " -nt 6 -wbt -bb 1000 -alrt 1000 -m MFP+MERGE"

			process = subprocess.run(bash_command.split(), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
			print(process.stdout)
			print(process.stderr)



# concatenate and run IQtree on both organelles
def run_both_organelles_combination(
		OUTDIR,
		N_GENES,
		IQTREE,
		gene_lists,
		outfolder_name1,
		outfolder_name2
	):

	"""
	I taB.
		---
	seq_DB : dict
		seq_DB database with placeholders for list of fasta sequences, alignment length, occupancy matrix and final fasta sequence
	file_list : list
		list of fasta alignment to be concatenated
	"""	

	# generate output folder
	os.makedirs("/".join([OUTDIR, outfolder_name1]), exist_ok = True)

	# combine cp and mt genes
	allgenes = gene_lists["CP"]["genes"] + gene_lists["MT"]["genes"]

	# generate genes combinations
	combinations_DB = {}
	for k in range(2, int(N_GENES) + 1):
		
		current_dir = "/".join([OUTDIR, outfolder_name1, outfolder_name2 + str(k) + "_genes/"])
		os.makedirs(current_dir, exist_ok = True)
		
		combinations_DB[k] = list(itertools.combinations(allgenes, k))

		for combination in combinations_DB[k]:
			
			# create sequences database
			seq_DB = {
				"seq_list"   : [],       # list of unique fasta sequence identified
				"seq_length" : [],       # length of alignments in residues [156788, 2222, 3435]
				"occupancy"  : {},       # [0, 1] presence of sequences for each multifasta alignment
				"sequences"  : {}        # placeholder for multifasta sequences
			}

			# get files and gene names
			file_list = combination
			gene_list = [gene.replace(".aln.fa", "") for gene in file_list]
			gene_list = [gene.replace(gene_lists["CP"]["DIR"] + "/", "") for gene in gene_list]
			gene_list = [gene.replace(gene_lists["MT"]["DIR"] + "/", "") for gene in gene_list]
			
			# import files and process them
			seq_DB = concatenate_multifasta(seq_DB, file_list)

			# print matrix of sequence occupancies
			outfile = open("/".join([current_dir, "_".join(gene_list) + ".concat.log"]), "w+")
			outfile.write("\t".join(["# Sequences", "\t".join(file_list), "\n"]))
			for key, value in seq_DB["occupancy"].items():
				outfile.write("\t".join([key, "\t".join(value), "\n"]))
			outfile.close()

			# print oncatenated multifasta alignment
			outfile = open("/".join([current_dir, "_".join(gene_list) + ".aln.fa"]), "w+")
			for sequence in sorted(seq_DB["seq_list"]):
				outfile.write(">" + sequence + "\n")
				outfile.write(seq_DB["sequences"][sequence] + "\n")
			outfile.close()

			# run IQtree
			bash_command = IQTREE + " -s " + "/".join([current_dir, "_".join(gene_list) + ".aln.fa"]) \
				+ " -st DNA -pre " + "/".join([current_dir, "_".join(gene_list) + ".aln"]) \
				+  " -nt 6 -wbt -bb 1000 -alrt 1000 -m MFP+MERGE"

			process = subprocess.run(bash_command.split(), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
			print(process.stdout)
			print(process.stderr)



# import and concatenate the multifasta alignment
def concatenate_multifasta(
		seq_DB,
		file_list
	):

	"""
	I take the seq_DB and file_list.
	I iterate multifasta alignment.
	For each alignment, I check if the fasta sequence was already encountered or not.
	New fasta sequences are added to list of fasta sequences and occupancy matrix.
	I add the residues to the growing fasta sequences.
	I add trailing or starting missing residues where necessary ("-").
	I return an updated seq_DB.
		---
	seq_DB : dict
		seq_DB database with placeholders for list of fasta sequences, alignment length, occupancy matrix and final fasta sequence
	file_list : list
		list of fasta alignment to be concatenated
	"""		

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
				seq_DB["seq_length"].append(len(seq_DB["sequences"][last_seq]) - list(itertools.accumulate(seq_DB["seq_length"]))[-1])
			
			# iterate sequences in multifasta and add "-" at the end of the sequences if necessary
			longest_sequence = max(len(item) for item in seq_DB["occupancy"].values())
			for sequence in seq_DB["seq_list"]:
				if len(seq_DB["occupancy"][sequence]) < longest_sequence:
					seq_DB["occupancy"][sequence].append("0")
					seq_DB["sequences"][sequence] += "-" * seq_DB["seq_length"][-1]
					
	# return updated sequence database
	return(seq_DB)



#------------------------------------------------------------------#
# RUN

# run script and give the running time
if __name__ == '__main__':
	t0 = datetime.now()
	main()
	dt = datetime.now() - t0
	sys.stderr.write("# Time elapsed: %s\n" % dt)
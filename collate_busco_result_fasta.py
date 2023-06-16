#!/usr/bin/env python
#
# collate_busco_result_fasta.py created 2023-01-19

"""collate_busco_result_fasta.py  last modified 2023-01-19

./collate_busco_result_fasta.py -s species_to_dir_list.tab -o combined_buscos/

    species list is tab delimited of species and directory, like:

Aphrocallistes_vastus	Avas.v1.29_annotations.prot.fr.fasta.busco5
Amphimedon_queenslandica	Aqu2.1_Genes_proteins.fasta.busco5
Ephydatia_muelleri	Emu_augustus_sysnames_prots.fasta.busco5
Oopsacas_minuta	JAKMXF01.1.gbff.prot.fa.busco5

"""

import sys
import os
import time
import argparse
import glob

def get_species_list(species_file):
	"""read species list file and return list in order"""
	species_list = []
	main_dir_list = []
	sys.stderr.write("# Reading species list from {}  {}\n".format( species_file, time.asctime() ) )
	for line in open(species_file,'r'):
		line = line.strip()
		if line:
			lsplits = line.split("\t")
			if len(lsplits) != 2:
				sys.stderr.write("# ERROR: line does not have two columns\n{}\n".format( line ) )
			species_list.append(lsplits[0])
			main_dir_list.append(lsplits[1])
	sys.stderr.write("# Found {} species  {}\n".format( len(species_list), time.asctime() ) )
	return species_list, main_dir_list

def get_fasta_from_dir(species_name, main_dir, output_dir):
	"""glob the directory and take the names of all fasta files, and return this list"""
	full_busco_path = os.path.join(main_dir, "run_metazoa_odb10", "busco_sequences", "single_copy_busco_sequences", "*.faa")
	fasta_list = glob.glob(full_busco_path)
	write_count = 0

	sys.stderr.write("# Found {} .faa files for species {}  {}\n".format( len(fasta_list), species_name, time.asctime() ) )
	for faafile in fasta_list:
		output_fasta = os.path.join(output_dir, os.path.basename(faafile))
		with open(faafile, 'r') as rf , open(output_fasta,'a') as wf :
			for line in rf:
				if line[0]==">":
					outline = line.replace('>',">{}|".format(species_name))
				else:
					outline = line
				wf.write(outline)
		write_count += 1
	sys.stderr.write("# Wrote {} sequences for species {}  {}\n".format( write_count, species_name, time.asctime() ) )
	# no return

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-s','--species-list', help="tab-delimited file of species and BUSCO directories, one species per line, folders should be named XX.busco5 or similar")
	parser.add_argument('-o','--output-directory', default="./", help="directory to collect collated output files")
	args = parser.parse_args(argv)

	species_list, busco_dir_list = get_species_list(args.species_list)
	for species, busco_dir in zip(species_list,busco_dir_list):
		get_fasta_from_dir(species, busco_dir, args.output_directory)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)

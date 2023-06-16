#!/usr/bin/env python
#
# combine_busco_table_results.py created 2023-01-30

"""combine_busco_table_results.py  last modified 2023-01-30

./combine_busco_table_results.py -s species_to_dir_list.tab > group_combined_busco_table.tsv

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
from collections import defaultdict

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

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-s','--species-list', help="tab-delimited file of species and BUSCO directories, one species per line, folders should be named XX.busco5 or similar")
	args = parser.parse_args(argv)

	species_list, busco_dir_list = get_species_list(args.species_list)
	genes_by_species = defaultdict( lambda: defaultdict(list) ) # key is species, then key is gene ID, then value is completeness

	gene_order = [] # recorded once

	species_count = 0
	for species, busco_dir in zip(species_list,busco_dir_list):
		full_busco_path = os.path.join(busco_dir, "run_metazoa_odb10", "full_table.tsv")
		species_count += 1
		sys.stderr.write("# Reading result table {} files for species {}  {}\n".format( full_busco_path, species, time.asctime() ) )

		# BUSCO version is: 5.4.2 				
		# The lineage dataset is: metazoa_odb10 (Creation date: 2021-02-17, number of genomes: 65, number of BUSCOs: 954)				
		# Busco id	Status	Sequence	Score	Length
		#5951at33208	Complete	Aqu2.1.29871_001	4199.2	2006
		#9639at33208	Complete	Aqu2.1.41859_001	2346.4	1766
		#13816at33208	Complete	Aqu2.1.36111_001	1509.6	1450

		for line in open(full_busco_path, 'r'):
			line = line.strip()
			if line and line[0]!="#": # skip comments
				lsplits = line.split("\t") 
				gene_id = lsplits[0]
				if species_count < 2:
					gene_order.append(gene_id)
				found_status = lsplits[1]
				if found_status=="Complete":
					found_score = "1"
				elif found_status=="Duplicated":
					found_score = "100"
				elif found_status=="Fragmented":
					found_score = "0.01"
				elif found_status=="Missing":
					found_score = "0"
				else: # should never happen
					found_score = "-1"
			genes_by_species[species][gene_id] = found_score

	output_header = "gene_id\t{}\n".format( "\t".join( species_list ) )
	sys.stdout.write( output_header )

	write_count = 0
	for gene_id in gene_order:
		out_scores_list = []
		for species in species_list:
			out_score = genes_by_species.get(species).get(gene_id,"NA")
			out_scores_list.append(out_score)
		output_line = "{}\t{}\n".format( gene_id, "\t".join( out_scores_list ) )
		write_count += 1
	sys.stderr.write("# Wrote {} genes for {} species  {}\n".format( write_count, len(species_list), time.asctime() ) )

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)

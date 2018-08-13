#!/usr/bin/env python
#
# remove_constant_breaks.py created by WRF 2018-04-16

'''remove_constant_breaks.py v1.2 2018-08-13
tool to quickly check for abnormal sequences in fasta alignments

remove_constant_breaks.py -a matrix.aln -p partitions.txt -o matrix_no_breaks.aln > matrix_const_break_results.tab

    for partitioned alignments, formats (-f) include:
  clustal, fasta, nexus, phylip, phylip-relaxed, stockholm

    large matrices can be gzipped, as .gz

    by default, ? is used to fill removed sites
    this may cause problems for some downstream programs
    change to normal gaps with -C "-"
'''

import sys
import argparse
import time
import gzip
from collections import defaultdict,Counter
from Bio import AlignIO

def get_partitions(partitionfile):
	'''read comma-delimited partition information and return a list of tuples'''
	partitions = [] # list of tuples of intervals
	for line in open(partitionfile,'r'):
		line = line.strip()
		if line:
			blocks = line.split(",") # split "1:136,137:301,..." into ['1:136', '137:301',...]
			for block in blocks:
				alignindex = tuple( int(i) for i in block.split(":") ) # split '1:136' into ( 1,136 )
				partitions.append(alignindex)
	print >> sys.stderr, "# read {} partitions from {}".format(len(partitions), partitionfile), time.asctime()
	return partitions

def merge_intervals(rangelist, offset=1):
	'''take list of intervals and return a list of merged intervals'''
	# if intervals were: (1, 5), (10, 17), (19, 28), (29, 34), (41, 51)
	# should return: (1, 5), (10, 34), (41, 51)
	newintervals = []
	srtrangelist = sorted(rangelist) # sort list now
	interval = srtrangelist[0]
	for bounds in srtrangelist:
		# since it is sorted bounds[0] should always be >= interval[0]
		if bounds[0] > interval[1]+offset: # if the next interval starts past the end of the first + 1
			newintervals.append(interval) # add to the new list, and continue with the next
			interval = bounds
		else: # meaning bounds[0] is less than interval[1]+1, should be merged
			if bounds[1] > interval[1]: # bounds[1] <= interval[1] means do not extend
				interval = (interval[0], bounds[1]) # otherwise extend the interval
	else: # append last interval
		newintervals.append(interval)
	return newintervals

def remove_constant_breaks(sequence, intervals, replacechar):
	seqlist = list(sequence)
	for interval in intervals:
		for i in range(interval[0], interval[1]):
			seqlist[i] = replacechar
	return seqlist

def count_breaks(fullalignment, alignformat, partitions, outputname, MAXBREAKTHRES, SKIPTHRES, REPLACECHAR):
	'''read large alignment, return two dicts where key is species and values are number of unbroken sequences and sum of breaks'''
	species_breaks = defaultdict(int) # total constant breaks by species
	new_species_list = defaultdict(list) # list of amino acids by species

	if fullalignment.rsplit('.',1)[1]=="gz": # autodetect gzip format
		opentype = gzip.open
		print >> sys.stderr, "# reading alignment {} as gzipped".format(fullalignment), time.asctime()
	else: # otherwise assume normal open
		opentype = open
		print >> sys.stderr, "# reading alignment {}".format(fullalignment), time.asctime()
	alignedseqs = AlignIO.read(opentype(fullalignment), alignformat)

	numtaxa = len(alignedseqs)
	alignlength = alignedseqs.get_alignment_length()
	print >> sys.stderr, "# Alignment contains {} taxa for {} sites, including gaps".format( numtaxa, alignlength )

	
	print >> sys.stderr, "# finding breaks with score > {}, allowing {} intervening sites".format(MAXBREAKTHRES, SKIPTHRES)
	totalremoved = 0
	for part in partitions:
		aa_freq_by_site = {} # key is site number, value is Counter
		const_breaker_thres = {} # key is site number, value is int of number of species / number of AAs
		most_common_by_site = {} # key is site number, value is most frequent letter

		alignpart = alignedseqs[:, part[0]-1:part[1] ] # alignment of each partition only
		partlength = alignpart.get_alignment_length()
		# get base frequency of each site
		for j in range(partlength):
			aligncolumn = str(alignpart[:,j])
			aa_freq = Counter(aligncolumn.replace("X","-").replace("?","-"))
			aa_freq_by_site[j] = aa_freq
			num_aas = len(aa_freq)
			breakmax = numtaxa/num_aas
			const_breaker_thres[j] = breakmax

			most_common_by_site[j] = aa_freq.most_common(1)[0][0]

		# count number of times each taxa breaks constant site rules
		for i,seqrec in enumerate(alignpart):
			pause_break = 0
			runningscore = 0
			previous_letter = '-'
			previous_common = '-'
			break_start = 0
			break_end = 0
			bound_list = [] # list of intervals as tuples
			consensus_string = ""
			run_string = ""

			species = seqrec.id
			seqlen = len(seqrec.seq)
			# if species is all gaps, assign as negative 1, for absent
			gapcount = str(seqrec.seq).replace("X","-").replace("?","-").count("-")
			if gapcount==seqlen:
				fixedseq = ["-"]*seqlen
				new_species_list[species].extend(fixedseq)
				continue
			# check frequency of each letter, and sum of inverse frequencies of long runs
			for k,letter in enumerate(str(seqrec.seq)):
				if letter=="X" or letter=="?":
					letter = "-"
				# break at long gaps, either for this species, or if most common part is gap
				if (letter=="-" and letter==previous_letter) or (most_common_by_site[k]=='-' and previous_common=='-'):
					runningscore, consensus_string, run_string, pause_break = 0, "", "", 0
				# only keep running total if site is rare
				# meaning frequency of letter is less than N/n
				# where N is number of taxa, and n is number of AAs at that site
				elif aa_freq_by_site[k][letter] < const_breaker_thres[k]:
					if runningscore == 0: # indicate start of new run
						break_start = k
					inverse_freq = 1.0 / aa_freq_by_site[k][letter]
					runningscore += inverse_freq
					consensus_string += most_common_by_site[k]
					run_string += letter
					pause_break = 0
				# break at frequent or constant sites
				elif letter=="P" or letter=="C": # never break at proline and cysteine
					pause_break += 1
					inverse_freq = 1.0 / aa_freq_by_site[k][letter]
					runningscore += inverse_freq
					consensus_string += most_common_by_site[k]
					run_string += letter
				elif runningscore > (MAXBREAKTHRES - 1) and pause_break < SKIPTHRES:
					pause_break += 1
					runningscore += inverse_freq
					consensus_string += most_common_by_site[k]
					run_string += letter
				else:
					if runningscore >= (MAXBREAKTHRES + 1):
						break_end = k - pause_break
						trim_bounds = (break_start, break_end)
						bound_list.append(trim_bounds)
						consensus_string = consensus_string[:-pause_break]
						run_string = run_string[:-pause_break]
						print >> sys.stdout, "{}\t{}\t{}\t{}\t{}\t{}\t{}".format( species, part, len(run_string), trim_bounds, runningscore, consensus_string, run_string )
						totalremoved += len(run_string)
					runningscore, consensus_string, run_string, pause_break = 0, "", "", 0
				previous_letter = letter
				previous_common = most_common_by_site[k]
			if bound_list:
				bound_list = merge_intervals(bound_list)
				fixedseq = remove_constant_breaks( str(seqrec.seq), bound_list, REPLACECHAR)
			else:
				fixedseq = list( str(seqrec.seq) )
			new_species_list[species].extend(fixedseq)
	print >> sys.stderr, "# identified {} amino acids for removal".format(totalremoved)

	# begin writing new matrix
	if outputname:
		print >> sys.stderr, "# writing filtered matrix to {}".format(outputname)
		with open(outputname,'w') as fm:
			for seqrec in alignedseqs: # all species must be counted once
				species = seqrec.id
				print >> fm, ">{}".format(species)
				outputstring = "".join(new_species_list[species])
				print >> fm, outputstring
	# no return

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-a','--alignment', help="supermatrix alignment")
	parser.add_argument('-B','--break-limit', type=int, default=4, help="min consecutive score to count motif as constant breaker")
	parser.add_argument('-C','--character', default="?", help="character to replace motifs [?]")
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-H','--header', action="store_true", help="include header line")
	parser.add_argument('-p','--partition', help="partition file for splitting large alignments")
	parser.add_argument('-o','--output', help="name for filtered matrix")
	parser.add_argument('-S','--skip-limit', type=int, default=1, help="longest run of conserved AAs to skip [1]")
	args = parser.parse_args(argv)

	partitions = get_partitions(args.partition)
	count_breaks(args.alignment, args.format, partitions, args.output, args.break_limit, args.skip_limit, args.character)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)

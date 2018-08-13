#!/usr/bin/env python
#
# check_supermatrix_alignments.py created 2017-03-13

'''check_supermatrix_alignments.py v1.3 2018-04-17
tool to quickly check for abnormal sequences in fasta alignments

check_supermatrix_alignments.py -a matrix.phy -p partitions.txt

    for partitioned alignments, formats (-f) include:
  clustal, fasta, nexus, phylip, phylip-relaxed, stockholm

    large matrices can be gzipped, as .gz

    output consists of tab delimited fields:
Species  Partitions  Number-missing  Percent-missing  Number-partial  Percent-partial

    for optional occupancy matrix file:
check_supermatrix_alignments.py -a matrix.phy -p partitions.txt -m occupancy_matrix.tab
    where matrix consists of three values for:
    present (2), partial (1), and absent (0)

    percent coverage (instead of binary) can be instead printed with --percent

    to check for long runs of infrequent amino acids (likely from misassembly)
    add option -b, partitions must be specified with -p
    this cannot be run with --percent
check_supermatrix_alignments.py -b -a matrix.aln -p partitions.txt -m break_matrix.tab

    matrix order can be changed based on a phylogenetic tree in nexus format
check_supermatrix_alignments.py -a matrix.phy -p partitions.txt -m occupancy_matrix.tab -T tree.nex

    to optionally include gene names in the matrix, use --pair-stats
    pair stats file is the output of align_pair_stats.py, or any file as:
1-1000    GENE
    or can include species or SwissProt ID
Genus_species_1-1000    sp|123456|GENE
'''

import sys
import os
import argparse
import time
import gzip
from collections import defaultdict,Counter
from Bio import AlignIO
from Bio import Phylo

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

def parts_to_genes(pairstatsfile):
	'''read tabular pair-wise gene stats, return a dict where key is partition and value is gene'''
	part_to_gene = {}
	print >> sys.stderr, "# reading partitions and gene names from {}".format(pairstatsfile), time.asctime()
	for line in open(pairstatsfile,'r'):
		line = line.strip()
		if line:
			lsplits = line.split("\t")
			if lsplits[0]=="partition": # skip first line
				continue
			gene = lsplits[1].split("|")[-1].replace("_HUMAN","")
			partition = tuple( int(i) for i in lsplits[0].split("_")[-1].split("-") )
			part_to_gene[partition] = gene
	print >> sys.stderr, "# found {} gene names".format(len(part_to_gene)), time.asctime()
	return part_to_gene

def check_alignments(fullalignment, alignformat, partitions, makematrix=False, writepercent=False):
	'''read large alignment, return the dict where key is species and value is number of gap-only sequences'''
	gapdict = {} # must set all values to zero in order to not skip full taxa
	halfgapdict = defaultdict(int)
	totaloccs = defaultdict(int)

	occmatrix = [] if makematrix else None

	if fullalignment.rsplit('.',1)[1]=="gz": # autodetect gzip format
		opentype = gzip.open
		print >> sys.stderr, "# reading alignment {} as gzipped".format(fullalignment), time.asctime()
	else: # otherwise assume normal open
		opentype = open
		print >> sys.stderr, "# reading alignment {}".format(fullalignment), time.asctime()
	alignedseqs = AlignIO.read(opentype(fullalignment), alignformat)
	print >> sys.stderr, "# Alignment contains {} taxa for {} sites, including gaps".format( len(alignedseqs), alignedseqs.get_alignment_length() )
	for seqrec in alignedseqs: # all species must be counted once
		gapdict[seqrec.id] = 0
		if occmatrix is not None:
			occmatrix.append([seqrec.id]) # initate each species as a new list with ID
	if partitions is None: # for cases where partitions are not given
		partitions = [ tuple([1,alignedseqs.get_alignment_length()]) ]
	for part in partitions:
		alignpart = alignedseqs[:, part[0]-1:part[1] ] # alignment of each partition only
		for i,seqrec in enumerate(alignpart):
			species = seqrec.id
			seqlen = len(seqrec.seq)
			lettercounts = Counter(str(seqrec.seq).replace("X","-"))
			gapcount = lettercounts.get("-",0) + lettercounts.get("?",0)
			if writepercent: # write integer of percent covered by each gene
				occupancyscore = 100 - 100*gapcount/seqlen
				gapdict[species] += gapcount
				halfgapdict[species] += seqlen
			else: # meaning stick with absent-partial-present scheme
				occupancyscore = 2 # by default is present, reassign if absent or partial
				if gapcount == seqlen: # seq is all gaps, so no seq
					gapdict[species] += 1
					occupancyscore = 0 # set to 0 if all gaps
				elif gapcount >= seqlen * 0.5: # partial means half or more of sequence is gaps
					halfgapdict[species] += 1
					occupancyscore = 1 # set to 1 if partial
			totaloccs[occupancyscore] += 1

			if occmatrix: # if building matrix, add that value to the matrix at sequence i
				occmatrix[i].append(str(occupancyscore))

	totalspots = sum(totaloccs.values())
	if writepercent:
		print >> sys.stderr, "# Matrix has {} ({:.2f}%) complete and {} ({:.2f}%) empty out of {} total genes".format( totaloccs[100], 100.0*totaloccs[100]/totalspots, totaloccs[0], 100.0*totaloccs[0]/totalspots, totalspots ), time.asctime()
	else:
		print >> sys.stderr, "# Matrix has {} ({:.2f}%) complete and {} ({:.2f}%) partial out of {} total genes".format( totaloccs[2], 100.0*totaloccs[2]/totalspots, totaloccs[1], 100.0*totaloccs[1]/totalspots, totalspots ), time.asctime()
	return gapdict, halfgapdict, occmatrix

def count_breaks(fullalignment, alignformat, partitions, makematrix=False, MAXBREAKTHRES=2):
	'''read large alignment, return two dicts where key is species and values are number of unbroken sequences and sum of breaks'''
	species_breaks = defaultdict(int) # total constant breaks by species
	species_corrects = defaultdict(int) # number of genes with no breaks by species
	occmatrix = [] if makematrix else None

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

	for seqrec in alignedseqs: # all species must be counted once
		if occmatrix is not None:
			occmatrix.append([seqrec.id]) # initate each species as a new list with ID

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
			longest_break = 0 # float of sum of 1/f for the longest segment of uncommon AAs
			runningscore = 0
			total_breaks = 0
			previous_letter = '-'
			previous_common = '-'

			species = seqrec.id
			seqlen = len(seqrec.seq)
			# if species is all gaps, assign as negative 1, for absent
			gapcount = str(seqrec.seq).replace("X","-").replace("?","-").count("-")
			if gapcount==seqlen:
				occmatrix[i].append("-1")
				continue
			# check frequency of each letter, and sum of inverse frequencies of long runs
			for k,letter in enumerate(str(seqrec.seq)):
				if letter=="X" or letter=="?":
					letter = "-"
				# break at long gaps, either for this species, or if most common part is gap
				if (letter=="-" and letter==previous_letter) or (most_common_by_site[k]=='-' and previous_common=='-'):
					if runningscore > longest_break:
						longest_break = runningscore
					total_breaks += runningscore
					runningscore = 0
				# only keep running total if site is rare
				# meaning frequency of letter is less than N/n
				# where N is number of taxa, and n is number of AAs at that site
				elif aa_freq_by_site[k][letter] < const_breaker_thres[k]:
					inverse_freq = 1.0 / aa_freq_by_site[k][letter]
					runningscore += inverse_freq
				# break at frequent or constant sites
				elif letter=="P" or letter=="C": # never break at proline and cysteine
					inverse_freq = 1.0 / aa_freq_by_site[k][letter]
					runningscore += inverse_freq
				else:
					if runningscore > longest_break:
						longest_break = runningscore
					total_breaks += runningscore
					runningscore = 0
				previous_letter = letter
				previous_common = most_common_by_site[k]
			species_breaks[species] += total_breaks
			occmatrix[i].append("{:.2f}".format(longest_break))
			if total_breaks <= MAXBREAKTHRES:
				species_corrects[species] += 1
	return species_corrects, species_breaks, occmatrix

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-a','--alignment', help="supermatrix alignment")
	parser.add_argument('-b','--breaks', action="store_true", help="score consecutive breaks instead of counting coverage")
	parser.add_argument('-B','--break-limit', type=int, default=2, help="max consecutive score to count gene as unbroken")
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-H','--header', action="store_true", help="include header line")
	parser.add_argument('-m','--matrix-out', help="name for optional matrix-occupancy output file")
	parser.add_argument('-D','--matrix-delimiter', default="\t", help="delimiter for matrix file [default is tab]")
	parser.add_argument('-T','--matrix-tree', help="optional Nexus-format tree to reorder matrix")
	parser.add_argument('-p','--partition', help="partition file for splitting large alignments")
	parser.add_argument('--percent', action="store_true", help="report values in matrix by percent complete")
	parser.add_argument('--pair-stats', help="pair stats file for gene names, to use alternate matrix header format")
	args = parser.parse_args(argv)

	if args.partition is None: # turn on percent mode if there is only
		print >> sys.stderr, "# no partitions given, calculating percent of total dataset"
		args.percent = True
		partitions = None
	else:
		partitions = get_partitions(args.partition)
	genenames = parts_to_genes(args.pair_stats) if args.pair_stats else None

	# check for breaks
	if args.breaks:
		if args.percent:
			sys.exit("WARNING: PERCENT MODE --percent CANNOT BE USED WITH BREAKS -b")
		maindict, secondarydict, occmatrix = count_breaks(args.alignment, args.format, partitions, args.matrix_out, args.break_limit)
	# otherwise check for coverage
	else:
		maindict, secondarydict, occmatrix = check_alignments(args.alignment, args.format, partitions, args.matrix_out, args.percent)

	if args.matrix_out and occmatrix:
		print >> sys.stderr, "# writing matrix to {}".format(args.matrix_out), time.asctime()
		with open(args.matrix_out,'w') as mo:
			# generate header line
			if args.pair_stats:
				headerline = ["Species"] + ["{}--{}".format(part[0], genenames.get(part,"None")) for part in partitions]
			elif args.partition is None:
				headerline = ["Species","PercentOccupancy"]
			else:
				headerline = ["Species"] + ["{}-{}".format(*part) for part in partitions]
			print >> mo, args.matrix_delimiter.join(headerline)
			# print occupancy by each species
			if args.matrix_tree: # use order from a rooted nexus-format tree
				occdict = {} # convert occmatrix from list of lists to dict of lists
				for occrow in occmatrix:
					occdict[occrow[0]] = occrow[1:]
				tree = Phylo.read(args.matrix_tree,"nexus")
				for clade in tree.get_terminals():
					cleanname = str(clade.name).replace("'","").replace('"','')
					try:
						print >> mo, args.matrix_delimiter.join( [cleanname] + occdict[cleanname] )
					except KeyError:
						print >> sys.stderr, "WARNING: CANNOT FIND TAXA {} IN ALIGNMENT, SKIPPING".format(cleanname)
			else: # just use default order from the alignment
				for occbysplist in occmatrix:
					print >> mo, args.matrix_delimiter.join(occbysplist)

	if args.percent: # just print total percentage of sites
		if args.header:
			#                  0           1           2    3
			print >> wayout, "Species\tSites\tGaps\tG%"
		for k,v in sorted(maindict.iteritems(), reverse=True, key=lambda x: x[1]):
			print >> wayout, "{}\t{}\t{}\t{:.2f}".format(k, secondarydict[k], v, v*100.0/secondarydict[k])
	elif args.breaks:
		if args.header:
			#                  0           1           2    3     4
			print >> wayout, "Species\tPartitions\tComplete\tM%\tBreaksum"
		numparts = len(partitions)
		for k,v in sorted(maindict.iteritems(), reverse=True, key=lambda x: x[1]):
			print >> wayout, "{}\t{}\t{}\t{:.2f}\t{}".format(k, numparts, v, v*100.0/numparts, secondarydict[k])
	else: # only print normal output if not in percent mode
		if args.header:
			#                  0           1           2    3     4      5
			print >> wayout, "Species\tPartitions\tMissing\tM%\tPartial\tP%"
		numparts = len(partitions)
		for k,v in sorted(maindict.iteritems(), reverse=True, key=lambda x: x[1]):
			print >> wayout, "{}\t{}\t{}\t{:.2f}\t{}\t{:.2f}".format(k, numparts, v, v*100.0/numparts, secondarydict[k], secondarydict[k]*100.0/numparts)
if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)

#! /usr/bin/env python
# add_taxa_to_align.py v1.0 created 2016-12-08

'''
add_taxa_to_align.py v1.5 2018-04-17
    add new taxa to an existing untrimmed alignment
    requires Bio Python library
    get hmmbuild and hmmscan (from hmmer package at http://hmmer.org/)

  IT IS ADVISED TO COLLECT ALL STDERR AS A LOG FILE WITH 2>

    to add proteins from species1 and species2 to alignments prot1 and prot2:
add_taxa_to_align.py -a prot1.aln prot2.aln -t species1.fasta species2.fasta
    more generally as:
add_taxa_to_align.py -a *.aln -t new_transcriptomes/ ...

    for very long lists of new taxa, instead use -X (--tabular-taxa)
    where options of fasta files (-t) and taxon names (-T) are both omitted
    and instead the information is put into a tab-separated file, as:
species1.fasta    Homo_sapiens
species2.fasta    Mus_musculus
    comment lines within the file are allowed with #

    to use a supermatrix with partitions as the input instead of alignments:
add_taxa_to_align.py -i partitions.txt -a big_alignment.aln ...

    for partitioned alignments, formats include:
  clustal, fasta, nexus, phylip, phylip-relaxed, stockholm

    partition file is comma-delimited text, as single or multiple lines
1:136,137:301,...

    add taxon names with -T, in the same order as -t
add_taxa_to_align.py -t species1.fasta species2.fasta -T Homo_sapiens Mus_musculus

    specify the new super matrix name with -U
add_taxa_to_align.py -U big_alignment_w_Hs_Mm.aln ...

    no e-value threshold is set for hmmsearch, as results are filtered later
    e-value cutoff will be determined automatically (default)
    this can be statically assigned with --ev-threshold,
    though that is not recommended

    for highly trimmed alignments
    (meaning each individual gene alignment is many discontinuous pieces)
    use the option --fragmented to set an alternate gap penalty for hmmbuild
    it may also be advisable to relax the --ev-correction to 1e20 or 1e25
'''

import sys
import os
import argparse
import time
import subprocess
from collections import defaultdict
from glob import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO

def get_partitions(partitionfile, errorlog):
	'''read comma-delimited partition information and return a list of tuples'''
	partitions = [] # list of tuples of intervals
	for line in open(partitionfile,'r'):
		line = line.strip()
		if line:
			blocks = line.split(",") # split "1:136,137:301,..." into ['1:136', '137:301',...]
			for block in blocks:
				alignindex = tuple( int(i) for i in block.split(":") ) # split '1:136' into ( 1,136 )
				partitions.append(alignindex)
	print >> errorlog, "# read {} partitions from {}".format(len(partitions), partitionfile), time.asctime()
	return partitions

def tabular_taxa_to_lists(tabulartaxafile):
	'''read tab-separated file and return two lists, one of fasta files, the other of taxon names'''
	fastafiles = []
	taxonnames = []
	for line in open(tabulartaxafile,'r'):
		line = line.strip()
		if line and line[0]!="#": # ignore blank and comment lines
			lsplits = line.split("\t")
			if len(lsplits)!=2:
				sys.exit("ERROR: INCORRECT FORMAT IN LINE:\n{}".format(line))
			fastafiles.append(os.path.expanduser(lsplits[0]))
			taxonnames.append(lsplits[1])
	return fastafiles, taxonnames

def make_alignments(fullalignment, alignformat, partitions, partitiondir, errorlog):
	'''split large alignment into individual alignments, return the list of files'''
	splitalignments = [] # list of filenames
	alignedseqs = AlignIO.read(fullalignment, alignformat)
	for part in partitions:
		alignpartname = os.path.join(partitiondir, "{}_{}_{}_part.aln".format( os.path.splitext(os.path.basename(fullalignment))[0], part[0], part[1] ) )
		alignpart = alignedseqs[:, part[0]-1:part[1] ]
		with open(alignpartname, 'w') as ao:
			AlignIO.write(alignpart, ao, "fasta")
		splitalignments.append(alignpartname)
	print >> errorlog, "# split alignment by partitions", time.asctime()
	return splitalignments

def run_hmmbuild(HMMBUILD, alignmentfile, errorlog, fragmented=False):
	'''generate HMM profile from multiple sequence alignment and return HMM filename'''
	# filename contains relative path, so should include partitions folder
	hmm_output = "{}.hmm".format(os.path.splitext(alignmentfile)[0] )
	if fragmented:
		hmmbuild_args = [HMMBUILD, "--pextend", "0.8", hmm_output, alignmentfile]
	else:
		hmmbuild_args = [HMMBUILD, hmm_output, alignmentfile]
	print >> errorlog, "#TIME {}\n{}".format(time.asctime(), " ".join(hmmbuild_args) )
	subprocess.call(hmmbuild_args, stdout=errorlog)
	print >> errorlog, "# hmmbuild completed", time.asctime()
	if os.path.isfile(hmm_output):
		return hmm_output
	else:
		raise OSError("Cannot find expected output file {}".format(hmm_output) )

def run_hmmsearch(HMMSEARCH, hmmprofile, fastafile, threadcount, hmmlog, hmmdir, errorlog):
	'''search fasta format proteins with HMM profile and return formatted-table filename'''
	hmmtbl_output = os.path.join(hmmdir, os.path.basename("{}_{}.tab".format(os.path.splitext(fastafile)[0], os.path.splitext(os.path.basename(hmmprofile))[0] ) ) )
	hmmsearch_args = [HMMSEARCH,"--cpu", str(threadcount), "--domE", "0.0001", "--domtblout", hmmtbl_output, hmmprofile, fastafile]
	print >> errorlog, "#TIME {}\n{}".format(time.asctime(), " ".join(hmmsearch_args) )
	if hmmlog:
		with open(hmmlog) as hmmstdout:
			subprocess.call(hmmsearch_args, stdout=hmmstdout)
	else:
		DEVNULL = open(os.devnull, 'w')
		subprocess.call(hmmsearch_args, stdout=DEVNULL)
	print >> errorlog, "# hmmsearch completed", time.asctime()
	if os.path.isfile(hmmtbl_output):
		return hmmtbl_output
	else:
		raise OSError("Cannot find expected output file {}".format(hmmtbl_output) )

def hmmtable_to_seqids(hmmtable, evaluecutoff, bitlencutoff, seqdict, verbose ):
	'''parse hits from hmm tblout and return a list of kept protein IDs'''
	# here domain-wise table output is processed, so that length of the hit can be directly checked
	# such information is not kept in the regular table output
	seqids_to_keep = {} # key is seqid, value is bits per length
	maxscore = 0
	maxbpl = 0
	lastevalue = 1e-300
	lastspan = 1 # span of HMM
	lastbits = 0
	hitcounter = 0
	for line in open(hmmtable, 'r').readlines():
		line = line.strip()
		if not line or line[0]=="#": # skip comment lines
			continue # also catch for empty line, which would cause IndexError
		lsplits = line.split(None,22)
		hitcounter += 1
		targetname = lsplits[0]
	#	evalue = float(lsplits[4])
	#	bitscore = float(lsplits[5])
		querylen = int(lsplits[5])
		evalue = float(lsplits[11])
		bitscore = float(lsplits[13])
		alignlength = int(lsplits[18]) - int(lsplits[17]) + 1
	#	alignlength = len(seqdict[targetname].seq)
		hmmspan = int(lsplits[16]) - int(lsplits[15]) + 1
		bitsperlen = bitscore/alignlength
		if bitscore > maxscore:
			maxscore = bitscore
		if bitsperlen > maxbpl:
			maxbpl = bitsperlen
		if bitsperlen + 0.1 < maxbpl:
			if verbose:
				print >> sys.stderr, "# REMOVING {}: BpL {} + 0.1 < MAX {}".format(targetname, bitsperlen, maxbpl )
			continue
		if evalue > evaluecutoff and bitsperlen < bitlencutoff:
			if verbose:
				print >> sys.stderr, "# REMOVING {}: EVALUE {} > {} AND BpL {} < {}".format( targetname, evalue, evaluecutoff, bitsperlen, bitlencutoff )
			continue
		if bitscore < maxscore*0.5: # discard any seqs that score less than half max
			if verbose:
				print >> sys.stderr, "# REMOVING {}: BITS {} < {} * 0.5".format( targetname, bitscore, maxscore )
			continue
		if evalue==0: # change values of 0 to 1e-300 for comparisons
			evalue = 1e-300
		# should better resolve paralogs
		if len(seqids_to_keep) >= 1 and evalue > lastevalue * 1e50: # if second best seq is 1e50 worse, skip
			if verbose:
				print >> sys.stderr, "# REMOVING {}: EVALUE {} > LAST {} * 1e50".format( targetname, evalue, lastevalue )
			continue
		# should correctly choose a longer splice variant if seqs are this close in evalue and bitscore
		if len(seqids_to_keep) >= 1 and hmmspan < lastspan and bitscore < lastbits - 100:
			if verbose:
				print >> sys.stderr, "# REMOVING {}: SPAN {} < LAST {} AND BITS {} < LAST {} - 100".format( targetname, hmmspan, lastspan, bitscore, lastbits )
			continue
		lastevalue = evalue
		lastspan = hmmspan
		lastbits = bitscore
		if verbose:
			print >> sys.stderr, "# CANDIDATE {}: EVALUE {}, BITS {}, BpL {}".format( targetname, evalue, bitscore, bitsperlen )
		if targetname in seqids_to_keep: # since multiple domains are allowed, keep highest scoring one
			if seqids_to_keep.get(targetname,0) < bitsperlen:
				seqids_to_keep[targetname] = bitsperlen
		else:
			seqids_to_keep[targetname] = bitsperlen
	# sort IDs by highest bits-per-length
	sorted_ids = [ k for k,v in sorted(seqids_to_keep.items(), key=lambda x: x[1], reverse=True ) ]
	if seqids_to_keep: # meaning if any hits
		print >> sys.stderr, "# retaining {}/{} seqs from {}, highest bits/len is {:.4f} for {}".format( len(seqids_to_keep) , hitcounter, os.path.basename(hmmtable), maxbpl, sorted_ids[0] )
		if maxbpl != seqids_to_keep[sorted_ids[0]]: # in case best reasonable seq is not the max
			print >> sys.stderr, "# WARNING: MAX {} DOES NOT MATCH TOP SEQ {} WITH {}".format( maxbpl, sorted_ids[0] , seqids_to_keep[sorted_ids[0]] )
	else: # meaning no hits
		print >> sys.stderr, "# keeping NONE of {} hits, highest score from {} was {:.4f} W/ cutoff {:.4f}".format( hitcounter, os.path.basename(hmmtable), maxbpl, bitlencutoff )
	return sorted_ids

def make_seq_length_dict(sequencefile):
	print >> sys.stderr, "# Parsing target sequences from {}".format(sequencefile), time.asctime()
	lengthdict = {}
	for seqrec in SeqIO.parse(sequencefile,'fasta'):
		lengthdict[seqrec.id] = len(seqrec.seq)
	print >> sys.stderr, "# Found {} sequences".format(len(lengthdict)), time.asctime()
	return lengthdict

def get_evalue_from_hmm(HMMSEARCH, hmmprofile, alignment, threadcount, hmmevaluedir, errorlog, evcorrection, bplcorrection):
	'''get dynamic evalue and bitscore/length cutoff from alignment and hmm'''
#                                                               --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# 0                  1          2                    3            4        5      6      7        8      9
# target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#10 11  12   13 14  15  16  17  18
#------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------

	# DEGAP FASTA FILE
	fasta_unaligned = os.path.join( hmmevaluedir, "{}_no_gaps.fasta".format(os.path.splitext(os.path.basename(alignment))[0] ) )
	unalign_sequences(fasta_unaligned, alignment, notrim=True, calculatemedian=False, removeempty=True)
	seqlendict = make_seq_length_dict(fasta_unaligned)

	# RUN HMMSEARCH
	# regular table output is used here, as all seqs are assumed to be correct
	hmmtbl_output = os.path.join( hmmevaluedir, os.path.basename( "{}_self_hmm.tab".format( os.path.splitext(alignment)[0] ) ) )
	hmmsearch_args = [HMMSEARCH,"--cpu", str(threadcount), "-E", "0.0001", "--tblout", hmmtbl_output, hmmprofile, fasta_unaligned]
	print >> errorlog, "{}\n{}".format(time.asctime(), " ".join(hmmsearch_args) )
	DEVNULL = open(os.devnull, 'w')
	subprocess.call(hmmsearch_args, stdout=DEVNULL)
	print >> errorlog, "# self-hmmsearch completed", time.asctime()

	# PROCESS RESULTS
	bitslenlist = [] # list of bits per length for long sequences using hmm results
	if os.path.isfile(hmmtbl_output):
		maxevalue = 0.0
		for line in open(hmmtbl_output, 'r').readlines():
			line = line.strip()
			if not line or line[0]=="#": # skip comment lines
				continue # also catch for empty line, which would cause IndexError
			lsplits = line.split(None,18)
			targetname = lsplits[0]
			evalue = float(lsplits[4])
			bitscore = float(lsplits[5])
			bitsperlength = bitscore / seqlendict[targetname]
			bitslenlist.append(bitsperlength)
			if evalue > maxevalue: # larger number, so smaller -log(x)
				maxevalue = evalue
		if maxevalue==0.0: # if all evalues are 0.0, then set to 1e-300
			maxevalue = 1e-300
		else:
			maxevalue = maxevalue * evcorrection # correct by constant margin of error, should be 1e10
		print >> errorlog, "# calculated e-value for {} as {:.3e}".format(alignment, maxevalue)
		bplcutoff = min(bitslenlist) - bplcorrection # should be 0.1 by default
		if bplcutoff < 0.9:
			bplcutoff = 0.9
			print >> errorlog, "# using minimum bits-per-length cutoff for {} of {:.2f}".format(alignment, bplcutoff), time.asctime()
		else:
			print >> errorlog, "# calculated bits-per-length cutoff for {} as {:.3f}".format(alignment, bplcutoff), time.asctime()
		return maxevalue, bplcutoff
	else:
		raise OSError("Cannot find expected output file {}".format(hmmtbl_output) )

def unalign_sequences(unalignedtaxa, alignment, notrim, calculatemedian, removeempty):
	'''remove gaps from alignments, and possibly return median ungapped length'''
	sizelist = []
	with open(unalignedtaxa,'w') as notaln:
		for seqrec in SeqIO.parse(alignment,"fasta"):
			gappedseq = str(seqrec.seq)
			degappedseq = Seq(gappedseq.replace("-","").replace("X",""))
			seqrec.seq = degappedseq
			if calculatemedian:
				sizelist.append(len(degappedseq))
			if removeempty and len(degappedseq) < 0.75*len(gappedseq):
				print >> sys.stderr, "# PARTIAL SEQUENCE {} REMOVED FROM HMM EVALUE TEST {}".format(seqrec.id, os.path.basename(alignment) )
				continue
			if notrim:
				notaln.write( seqrec.format("fasta") )
	if calculatemedian:
		median = sorted(sizelist)[len(sizelist)/2]
		return median
	else: # essentially void, just generates unaligned fasta file
		return None

def collect_sequences(unalignednewtaxa, alignment, hitlistolists, lengthcutoff, speciesnames, maxhits, dosupermatrix, notrim ):
	'''write sequences from old alignment and new hits to file'''
	speciescounts = defaultdict(int) # key is species, value is number of written seqs per species
	median = unalign_sequences(unalignednewtaxa, alignment, notrim, calculatemedian=True, removeempty=False)
	with open(unalignednewtaxa,'a') as notaln:
		# hitlistolists is a list of lists, so that order of species is preserved
		for i,hitlist in enumerate(hitlistolists):
			writeout = 0
			if not hitlist:
				print >> sys.stderr, "# NO HITS FOR {} IN {}".format(speciesnames[i], os.path.basename(alignment) )
				print >> notaln, ">{}".format(speciesnames[i])
				continue
			for seqrec in hitlist: # sublist, each item is a SeqRecord object
				if writeout==maxhits: # if already have enough candidates
					break
				if len(seqrec.seq) >= median*lengthcutoff: # remove short sequences
					if dosupermatrix:
						if seqrec.id==speciesnames[i]: # check if seq was already used, so dict entry was renamed
							print >> sys.stderr, "WARNING: REDUNDANT SEQ {} FOR {}".format(seqrec.name, os.path.basename(alignment) )
						print >> sys.stderr, "# using seq {} for {}".format( seqrec.id, speciesnames[i] )
						seqrec.id = str(speciesnames[i])
						seqrec.description = ""
					notaln.write( seqrec.format("fasta") )
					writeout += 1
				else: # meaning too short
					print >> sys.stderr, "# SEQ {} TOO SHORT FOR {} IN {}".format(seqrec.name, speciesnames[i], os.path.basename(alignment) )
			if writeout==0: # all hits missed the cut or had no hits, give a dummy entry
				print >> notaln, ">{}".format(speciesnames[i])
				print >> sys.stderr, "# ALL HITS TOO SHORT FOR {} IN {}".format(speciesnames[i], os.path.basename(alignment) )
	# no return

def run_mafft(MAFFT, rawseqsfile, errorlog):
	'''generate multiple sequence alignment from fasta and return MSA filename'''
	aln_output = "{}.aln".format(os.path.splitext(rawseqsfile)[0] )
	aligner_args = [MAFFT, "--maxiterate", "1000", "--localpair", "--quiet", rawseqsfile]
	print >> errorlog, "#TIME {}\n{} > {}".format(time.asctime(), " ".join(aligner_args), aln_output )
	with open(aln_output, 'w') as msa:
		subprocess.call(aligner_args, stdout=msa)
	print >> errorlog, "# alignment of {} completed".format(aln_output), time.asctime()
	if os.path.isfile(aln_output):
		return aln_output
	else:
		raise OSError("Cannot find expected output file {}".format(aln_output) )

def run_mafft_addlong(MAFFT, oldalignment, rawseqsfile, errorlog):
	'''generate new MSA from fasta and old MSA and return MSA filename'''
	aln_output = "{}.aln".format(os.path.splitext(rawseqsfile)[0] )
#	aligner_args = [MAFFT, "--quiet", "--keeplength", "--auto", "--addlong", rawseqsfile, oldalignment]
	aligner_args = [MAFFT, "--maxiterate", "1000", "--localpair", "--quiet", "--addlong", rawseqsfile, oldalignment]
	print >> errorlog, "#TIME {}\n{} > {}".format(time.asctime(), " ".join(aligner_args), aln_output )
	with open(aln_output, 'w') as msa:
		subprocess.call(aligner_args, stdout=msa)
	print >> errorlog, "# alignment of {} completed".format(aln_output), time.asctime()
	if os.path.isfile(aln_output):
		return aln_output
	else:
		raise OSError("Cannot find expected output file {}".format(aln_output) )

def run_tree(FASTTREEMP, alignfile, errorlog):
	'''generate tree from alignment'''
	tree_output = "{}.tree".format(os.path.splitext(alignfile)[0] )
	fasttree_args = [FASTTREEMP, "-quiet", alignfile]
	print >> errorlog, "#TIME {}\n{}".format(time.asctime(), " ".join(fasttree_args) )
	with open(tree_output, 'w') as tree:
		subprocess.call(fasttree_args, stdout=tree)
	if os.path.isfile(tree_output):
		return tree_output
	else:
		raise OSError("Cannot find expected output file {}".format(tree_output) )

def main(argv, wayout, errorlog):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-a','--alignments', nargs="*", help="alignment files or directory")
	parser.add_argument('-d','--directory', default="new_taxa", help="directory for new alignments [autonamed]")
	parser.add_argument('--ev-threshold', help="static E-value cutoff for hmmsearch, otherwise determined automatically")
	parser.add_argument('--ev-correction', type=float, default=1e15, help="E-value offset for automatic detection [1e15]")
	parser.add_argument('--bpl-threshold', type=float, default=0.9, help="bits-per-length minimum threshold, otherwise determined automatically")
	parser.add_argument('--bpl-correction', type=float, default=0.1, help="bits-per-length offset for automatic detection [0.1]")
	parser.add_argument('-E','--hmm-evalue-dir', default="hmm_vs_self", help="temporary directory for hmm self-search to dynamically determine E-value [./hmm_vs_self]")
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-i','--partition', help="optional partition file for splitting large alignments")
	parser.add_argument('-I','--partition-dir', default="partitions", help="temporary directory for partitioned alignments [./partitions]")
	parser.add_argument('-l','--length', type=float, default=0.4, help="minimum length cutoff compared to median [0.4]")
	parser.add_argument('-m','--max-hits', type=int, default=1, help="max number of allowed protein hits [1]")
	parser.add_argument('-p','--processors', type=int, default=1, help="number of processors [1]")
	parser.add_argument('-r','--no-trim', action="store_true", help="do not trim to original alignment in mafft")
	parser.add_argument('-s','--hmm-results', default=None, help="optional filename for hmmsearch output")
	parser.add_argument('-S','--hmm-dir', default="hmm_hits", help="temporary directory for hmm hits [./hmm_hits]")
	parser.add_argument('-t','--taxa', nargs="*", help="new taxa as fasta files of proteins (such as multiple translated transcriptomes) or directory")
	parser.add_argument('-T','--taxa-names', nargs="*", help="optional species names for supermatrix (can contain underscores, no spaces)")
	parser.add_argument('-X','--tabular-taxa', help="optional tab-separated file for fasta files and taxon names (must not use -t or -T)")
	parser.add_argument('-U','--supermatrix', help="name for optional supermatrix output")
	parser.add_argument('-v','--verbose', action="store_true", help="verbose output of HMM hit removal")
	parser.add_argument('--mafft', default="mafft", help="path to mafft binary [default is in PATH]")
	parser.add_argument('--hmmbin', default="", help="path to hmm binaries, should be a directory containing hmmbuild and hmmsearch [default is ./]")
	parser.add_argument('--fasttree', default="FastTreeMP", help="path to fasttree binary [default is in PATH]")
	parser.add_argument('--fragmented', action="store_true", help="alignments are highly-trimmed fragments, use relaxed gap properties")
	parser.add_argument('--no-gene-trees', action="store_true", help="skip the FastTree step")
	args = parser.parse_args(argv)

	ALIGNER = os.path.expanduser(args.mafft)
	HMMBUILD = os.path.expanduser(os.path.join(args.hmmbin, "hmmbuild"))
	HMMSEARCH = os.path.expanduser(os.path.join(args.hmmbin, "hmmsearch"))
	FASTTREEMP = os.path.expanduser(args.fasttree)

	starttime = time.time()
	print >> errorlog, "# Script called as:\n{}".format( ' '.join(sys.argv) )

	### PROTEIN FILES FOR NEW TAXA
	if args.taxa is not None:
		if os.path.isdir(args.taxa[0]):
			print >> errorlog, "# Finding protein files from directory {}".format(args.taxa[0]), time.asctime()
			globstring = "{}*".format(args.taxa[0])
			newtaxafiles = glob(globstring)
		elif os.path.isfile(args.taxa[0]):
			newtaxafiles = args.taxa
	### IN CASE TABULAR TAXA FILE IS USED
	elif args.tabular_taxa is not None:
		newtaxafiles, args.taxa_names = tabular_taxa_to_lists(args.tabular_taxa)
	else: # otherwise no taxa are given
		raise OSError("ERROR: Unknown new protein files, exiting")

	### CHECK IF ALL FILES EXIST
	missingfiles = []
	for protfile in newtaxafiles:
		if not os.path.isfile(protfile):
			missingfiles.append(protfile)
	if missingfiles: # if any are missing
			raise OSError("ERROR: Cannot find new protein files:\n{}\nexiting".format("\n".join(missingfiles)))

	### APPLY NAMES OF TAXA
	if args.taxa_names:
		if len(args.taxa_names)!=len(newtaxafiles):
			for tfile, tname in zip(newtaxafiles, args.taxa_names):
				print >> errorlog, "{}\t{}".format(tfile, tname)
			raise ValueError("ERROR: number of taxa names ({}) does not match number of files ({}), exiting".format(len(args.taxa_names),len(newtaxafiles) ) )
		else:
			print >> errorlog, "# Using {} pairs of files and taxon names".format(len(newtaxafiles))
			for tfile, tname in zip(newtaxafiles, args.taxa_names):
				print >> errorlog, "{}\t{}".format(tfile, tname)
	if len(args.taxa_names) != len(set(args.taxa_names)):
		raise ValueError("ERROR: duplicate taxa names, check -T, exiting")

	### TIME STRING FOR ALL NEW DIRECTORIES
	timestring = time.strftime("%Y%m%d-%H%M%S")

	### SINGLE PARTITIONED ALIGNMENT TO BE EXTENDED
	if args.partition: # if partitioning, do this first
		if len(args.alignments) > 1:
			raise OSError("ERROR: Expecting one alignment to partition, found {}, exiting".format( len(args.alignments) ) )
		elif not os.path.isfile(args.alignments[0]):
			raise OSError("ERROR: Cannot find {} alignment for partitions, exiting".format(args.alignments[0]) )
		else:
			partitions = get_partitions(args.partition, errorlog)
			part_dir = os.path.abspath( "{}_{}".format( timestring, args.partition_dir ) )
			if not os.path.isdir(part_dir):
				if os.path.isfile(part_dir):
					raise OSError("ERROR: Directory {} exists as a file, cannot create, exiting".format(part_dir) )
				else:
					os.mkdir(part_dir)
					print >> errorlog, "# Creating directory {}".format(part_dir), time.asctime()
			alignfiles = make_alignments(args.alignments[0], args.format, partitions, part_dir, errorlog)
	### ALIGNMENTS TO BE EXTENDED AS MULTIPLE FILES
	else: # otherwise treat alignments as normal, either directory or single file
		if os.path.isdir(args.alignments[0]):
			print >> errorlog, "# Reading alignments from directory {}".format(args.alignments[0]), time.asctime()
			globstring = "{}*".format(args.alignments[0])
			alignfiles = glob(globstring)
		elif os.path.isfile(args.alignments[0]):
			alignfiles = args.alignments
		else:
			raise OSError("ERROR: Unknown alignment files, exiting")

	### DIRECTORY FOR NEW OUTPUT
	new_aln_dir = os.path.abspath( "{}_{}".format( timestring, args.directory ) )
	if not os.path.exists(new_aln_dir):
		os.mkdir(new_aln_dir)
		print >> errorlog, "# Creating directory {}".format(new_aln_dir), time.asctime()
	elif os.path.isdir(new_aln_dir):
		print >> errorlog, "# Using directory {}".format(new_aln_dir), time.asctime()

	### DIRECTORY FOR HMM RESULTS ###
	new_hmm_dir = os.path.abspath( "{}_{}".format( timestring, args.hmm_dir ) )
	if not os.path.isdir(new_hmm_dir):
		if os.path.isfile(new_hmm_dir):
			raise OSError("ERROR: Directory {} exists as a file, cannot create, exiting".format(new_hmm_dir) )
		else:
			print >> errorlog, "# Creating directory {}".format(new_hmm_dir), time.asctime()
			os.mkdir(new_hmm_dir)

	### DIRECTORY FOR HMM SELF SEARCH ###
	new_self_dir = os.path.abspath( "{}_{}".format( timestring, args.hmm_evalue_dir ) )
	if not os.path.isdir(new_self_dir):
		if os.path.isfile(new_self_dir):
			raise OSError("ERROR: Directory {} exists as a file, cannot create, exiting".format(new_self_dir) )
		else:
			print >> errorlog, "# Creating directory {}".format(new_self_dir), time.asctime()
			os.mkdir(new_self_dir)

	### MAIN LOOP
	supermatrix = None
	partitionlist = [] # empty list for new partition file
	runningsum = 0

	hmmprofilelist = [] # store hmms as list of files to run for each taxa
	aligntohmm = {} # key is alignment file, value is name of hmm profile
	evalueforhmmdict = {} # evalue is calculated once, stored as value where hmm name is key
	bitlenforhmmdict = {} # bitscore per length is calculated once, hmm name is key

	print >> errorlog, "# Building HMM profiles", time.asctime()
	if args.ev_threshold is None:
		print >> errorlog, "# Determining evalue for each HMM, using offset of {:.1e}".format(args.ev_correction)
	else:
		try:
			args.ev_threshold = float(args.ev_threshold)
		except ValueError:
			print >> errorlog, "# WARNING: {} CANNOT BE CONVERTED TO FLOAT, LEAVE --ev-threshold BLANK TO ENABLE AUTO-DETERMINATION OF EVALUE".format(args.ev_threshold)
			args.ev_threshold = None

	# GET EVALUE AND BPL FOR EACH ALIGNMENT
	for alignment in alignfiles:
		# make hmm profile from alignment, test against original set to determine evalue cutoff for new seqs
		hmmprofile = run_hmmbuild(HMMBUILD, alignment, errorlog, args.fragmented)
		aligntohmm[alignment] = hmmprofile
		hmmprofilelist.append(hmmprofile)
		# check if threshold was assigned, otherwise determine for each hmm profile
		# by running hmmsearch against the original proteins
		filtered_evalue, filtered_bitlen = [args.ev_threshold, args.bpl_threshold] if args.ev_threshold is not None else get_evalue_from_hmm(HMMSEARCH, hmmprofile, alignment, args.processors, new_self_dir, errorlog, args.ev_correction, args.bpl_correction)
		# either use the static values, or the dynamic values determined above
		evalueforhmmdict[hmmprofile] = filtered_evalue
		bitlenforhmmdict[hmmprofile] = filtered_bitlen

	# RUN HMMSEARCH FOR EACH SPECIES
	# search each species sequentially, keep matches as SeqRecords in list of lists
	# to keep memory profile low, sec dict is remade for each species, then all HMMs are run
	seqrecs_to_add = defaultdict( list ) # key is hmm, value is list of lists of SeqRecords by species
	print >> errorlog, "# Beginning search for orthologs in new taxa", time.asctime()
	speciesnames = args.taxa_names if args.taxa_names else [os.path.basename(newspeciesfile) for f in newtaxafiles]
	for newspeciesfile in newtaxafiles:
		seqids_to_add = [] # build list of lists by species
		print >> errorlog, "# Reading proteins from {} into memory".format(newspeciesfile), time.asctime()
		seqdict = SeqIO.to_dict( SeqIO.parse(newspeciesfile, "fasta") )
		for hmmprof in hmmprofilelist:
			hmmtableout = run_hmmsearch(HMMSEARCH, hmmprof, newspeciesfile, args.processors, args.hmm_results, new_hmm_dir, errorlog)
			# append seqrecord to list for each hmm
			seqids_to_add = [ seqdict[hitseqid] for hitseqid in hmmtable_to_seqids(hmmtableout, evalueforhmmdict[hmmprof], bitlenforhmmdict[hmmprof] , seqdict, args.verbose) ] # is list of seqrecords
			# add this list for each hmm in the dict 'seqrecs_to_add'
			seqrecs_to_add[hmmprof].append(seqids_to_add)

	# GET BEST HITS AND MAKE NEW ALIGNMENT
	# then reiterate over each alignment, make the new fasta file, alignment, tree, and add to supermatrix
	for alignment in alignfiles:
		nt_unaligned = os.path.join(new_aln_dir, "{}.fasta".format(os.path.splitext(os.path.basename(alignment))[0] ) )
		collect_sequences(nt_unaligned, alignment, seqrecs_to_add[aligntohmm[alignment]], args.length, speciesnames, args.max_hits, args.supermatrix, args.no_trim)
		if args.no_trim: # use original method, which allows gaps from new sequences
			nt_aligned = run_mafft(ALIGNER, nt_unaligned, errorlog)
		else: # use --keeplength in mafft
			nt_aligned = run_mafft_addlong(ALIGNER, alignment, nt_unaligned, errorlog)

	# BUILD SUPERMATRIX
		# generate supermatrix from alignments
		newaligned = AlignIO.read(nt_aligned, "fasta")
		alignment_length = newaligned.get_alignment_length()
		partitionlist.append("{}:{}".format(runningsum+1, runningsum+alignment_length) )
		runningsum += alignment_length
		if supermatrix is None:
			supermatrix = newaligned
		else:
			supermatrix += newaligned
		if not args.no_gene_trees:
			nt_tree = run_tree(FASTTREEMP, nt_aligned, errorlog)

	### WRITE SUPERMATRIX AND PARTITIONS
	if args.supermatrix:
		AlignIO.write(supermatrix, args.supermatrix, "fasta")
		print >> errorlog, "# Supermatrix written to {}".format(args.supermatrix), time.asctime()
		with open("{}.partition.txt".format(args.supermatrix),'w') as pf:
			print >> pf, ",".join(partitionlist)
	print >> errorlog, "# Process completed in {:.1f} minutes".format( (time.time()-starttime)/60 )

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout, sys.stderr)

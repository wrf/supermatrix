#! /usr/bin/env python
# add_taxa_to_align.py v1.0 created 2016-12-08

'''
add_taxa_to_align.py v1.1 2017-05-17
    add new taxa to an existing untrimmed alignment

    to add proteins from species1 and species2 to alignments prot1 and prot2:
add_taxa_to_align.py -a prot1.aln prot2.aln -t species1.fasta species2.fasta

    more generally as:
add_taxa_to_align.py -a *.aln -t new_transcriptomes/

    requires Bio Python library

    get hmmbuild and hmmscan (from hmmer package at http://hmmer.org/)

    for partitioned alignments, formats include:
  clustal, fasta, nexus, phylip, phylip-relaxed, stockholm

    partition file is comma-delimited text, as single or multiple lines
1:136,137:301,...
'''

import sys
import os
import argparse
import time
import subprocess
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

def run_hmmbuild(HMMBUILD, alignmentfile, errorlog):
	'''generate HMM profile from multiple sequence alignment and return HMM filename'''
	hmm_output = "{}.hmm".format(os.path.splitext(alignmentfile)[0] )
	hmmbuild_args = [HMMBUILD, hmm_output, alignmentfile]
	print >> errorlog, "{}\n{}".format(time.asctime(), " ".join(hmmbuild_args) )
	subprocess.call(hmmbuild_args, stdout=errorlog)
	print >> errorlog, "# hmmbuild completed", time.asctime()
	if os.path.isfile(hmm_output):
		return hmm_output
	else:
		raise OSError("Cannot find expected output file {}".format(hmm_output) )

def run_hmmsearch(HMMSEARCH, hmmprofile, fastafile, threadcount, hmmlog, hmmdir, errorlog):
	'''search fasta format proteins with HMM profile and return formatted-table filename'''
	hmmtbl_output = os.path.join(hmmdir, os.path.basename("{}_{}.tab".format(os.path.splitext(fastafile)[0], os.path.splitext(os.path.basename(hmmprofile))[0] ) ) )
	hmmsearch_args = [HMMSEARCH,"--cpu", str(threadcount), "--tblout", hmmtbl_output, hmmprofile, fastafile]
	print >> errorlog, "{}\n{}".format(time.asctime(), " ".join(hmmsearch_args) )
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

def hmmtable_to_seqids(hmmtable, evaluecutoff, scorecutoff=0.5):
	'''parse hits from hmm tblout and return dict where value is list of kept protein IDs'''
#                                                               --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# 0                  1          2                    3            4        5      6      7        8      9
# target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#10 11  12   13 14  15  16  17  18
#------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
	seqids_to_keep = []
	maxscore = 0
	for line in open(hmmtable, 'r').readlines():
		line = line.strip()
		if not line or line[0]=="#": # skip comment lines
			continue # also catch for empty line, which would cause IndexError
		lsplits = line.split(None,18)
		targetname = lsplits[0]
		evalue = float(lsplits[4])
		bitscore = float(lsplits[5])
		if bitscore > maxscore:
			maxscore = bitscore
		if evalue <= evaluecutoff and bitscore > maxscore*scorecutoff:
			seqids_to_keep.append(targetname)
	return seqids_to_keep

def collect_sequences(unalignednewtaxa, alignment, hitlist, sequencedict, lengthcutoff, speciesnames, maxhits, dosupermatrix, notrim, verbose=False):
	'''write sequences from old alignment and new hits to file'''
	sizelist = []
	with open(unalignednewtaxa,'w') as notaln:
		for seqrec in SeqIO.parse(alignment,"fasta"):
			gappedseq = str(seqrec.seq)
			degappedseq = Seq(gappedseq.replace("-","").replace("X",""))
			seqrec.seq = degappedseq
			sizelist.append(len(degappedseq))
			if notrim:
				notaln.write( seqrec.format("fasta") )
		median = sorted(sizelist)[len(sizelist)/2]
		for i,hits in enumerate(hitlist):
			writeout = 0
			for seqid in hits:
				if writeout==maxhits: # already have enough candidates
					break
				hitrecord = sequencedict[seqid]
				if len(hitrecord.seq) >= median*lengthcutoff: # remove short sequences
					if dosupermatrix:
						hitrecord.id = str(speciesnames[i])
						hitrecord.description = ""
					notaln.write( hitrecord.format("fasta") )
					writeout += 1
			if writeout==0: # all hits missed the cut or had no hits, give a dummy entry
				print >> notaln, ">{}".format(speciesnames[i])
				if verbose:
					print >> sys.stderr, "NO HITS FOR {} IN {}".format(speciesnames[i], alignment)
	# no return

def run_mafft(MAFFT, rawseqsfile, errorlog):
	'''generate multiple sequence alignment from fasta and return MSA filename'''
	aln_output = "{}.aln".format(os.path.splitext(rawseqsfile)[0] )
	aligner_args = [MAFFT, "--auto", "--quiet", rawseqsfile]
	print >> errorlog, "{}\n{}".format(time.asctime(), " ".join(aligner_args) )
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
	aligner_args = [MAFFT, "--quiet", "--keeplength", "--auto", "--addlong", rawseqsfile, oldalignment]
	print >> errorlog, "{}\n{}".format(time.asctime(), " ".join(aligner_args) )
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
	print >> errorlog, "{}\n{}".format(time.asctime(), " ".join(fasttree_args) )
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
	parser.add_argument('-d','--directory', default="./new_taxa", help="directory for new alignments [autonamed]")
	parser.add_argument('-e','--evalue', type=float, default=1e-20, help="Evalue cutoff for hmmsearch [0.01]")
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-i','--partition', help="optional partition file for splitting large alignments")
	parser.add_argument('-I','--partition-dir', default="partitions", help="temporary directory for partitioned alignments [./partitions]")
	parser.add_argument('-l','--length', type=float, default=0.5, help="minimum length cutoff compared to median [0.5]")
	parser.add_argument('-m','--max-hits', type=int, default=1, help="max number of allowed protein hits [1]")
	parser.add_argument('-p','--processors', type=int, default=1, help="number of processors [1]")
	parser.add_argument('-r','--no-trim', action="store_true", help="do not trim to original alignment in mafft")
	parser.add_argument('-s','--hmm-results', default=None, help="optional filename for hmmsearch output")
	parser.add_argument('-S','--hmm-dir', default="hmm_hits", help="temporary directory for hmm hits [./hmm_hits]")
	parser.add_argument('-t','--taxa', nargs="*", help="new taxa as fasta files of proteins (such as multiple translated transcriptomes) or directory")
	parser.add_argument('-T','--taxa-names', nargs="*", help="optional species names for supermatrix (can contain underscores, no spaces)")
	parser.add_argument('-U','--supermatrix', help="name for optional supermatrix output")
	parser.add_argument('--mafft', default="mafft", help="path to mafft binary [default is in PATH]")
	parser.add_argument('--hmmbin', default="", help="path to hmm binaries, should be a directory containing hmmbuild and hmmsearch [default is ./]")
	parser.add_argument('--fasttree', default="FastTreeMP", help="path to fasttree binary [default is in PATH]")
	args = parser.parse_args(argv)

	ALIGNER = os.path.expanduser(args.mafft)
	HMMBUILD = os.path.expanduser(os.path.join(args.hmmbin, "hmmbuild"))
	HMMSEARCH = os.path.expanduser(os.path.join(args.hmmbin, "hmmsearch"))
	FASTTREEMP = os.path.expanduser(args.fasttree)

	### PROTEIN FILES FOR NEW TAXA
	if os.path.isdir(args.taxa[0]):
		print >> errorlog, "# Reading protein files from directory {}".format(args.taxa[0]), time.asctime()
		globstring = "{}*".format(args.taxa[0])
		newtaxafiles = glob(globstring)
	elif os.path.isfile(args.taxa[0]):
		newtaxafiles = args.taxa
	else:
		raise OSError("ERROR: Unknown new protein files, exiting")

	if args.taxa_names and len(args.taxa_names)!=len(args.taxa):
		raise ValueError("ERROR: number of taxa does not match number of files, exiting")

	### SINGLE PARTITIONED ALIGNMENT TO BE EXTENDED
	if args.partition: # if partitioning, do this first
		if len(args.alignments) > 1:
			raise OSError("ERROR: Expecting 1 alignment, found {}, exiting".format( len(args.alignments) ) )
		elif not os.path.isfile(args.alignments[0]):
			raise OSError("ERROR: Cannot find {} alignment for partitions, exiting".format(args.alignments[0]) )
		else:
			partitions = get_partitions(args.partition, errorlog)
			if not os.path.isdir(args.partition_dir):
				if os.path.isfile(args.partition_dir):
					raise OSError("ERROR: Cannot create directory {}, exiting".format(args.partition_dir) )
				else:
					os.mkdir(args.partition_dir)
			alignfiles = make_alignments(args.alignments[0], args.format, partitions, args.partition_dir, errorlog)
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

	### PROTEIN DICT
	seqdict = {}
	print >> errorlog, "# Reading proteins into memory", time.asctime()
	for protsfile in newtaxafiles:
		seqdict.update(SeqIO.to_dict(SeqIO.parse(protsfile, "fasta")))

	### DIRECTORY FOR NEW OUTPUT
	new_aln_dir = os.path.abspath("{}_{}".format(args.directory, time.strftime("%Y%m%d-%H%M%S") ) )
	if not os.path.exists(new_aln_dir):
		os.mkdir(new_aln_dir)
		print >> errorlog, "# Creating directory {}".format(new_aln_dir), time.asctime()
	elif os.path.isdir(new_aln_dir):
		print >> errorlog, "# Using directory {}".format(new_aln_dir), time.asctime()

	### DIRECTORY FOR HMM RESULTS ###
	if not os.path.isdir(args.hmm_dir):
		if os.path.isfile(args.hmm_dir):
			raise OSError("ERROR: Cannot create directory {}, exiting".format(args.hmm_dir) )
		else:
			os.mkdir(args.hmm_dir)

	### MAIN LOOP
	supermatrix = None
	partitionlist = [] # empty list for new partition file
	runningsum = 0
	for alignment in alignfiles:
		seqids_to_add = [] # build list of lists by species
		hmmprofile = run_hmmbuild(HMMBUILD, alignment, errorlog)
		for newspeciesfile in newtaxafiles:
			hmmtableout = run_hmmsearch(HMMSEARCH, hmmprofile, newspeciesfile, args.processors, args.hmm_results, args.hmm_dir, errorlog)
			seqids_to_add.append(hmmtable_to_seqids(hmmtableout, args.evalue))
		nt_unaligned = os.path.join(new_aln_dir, "{}.fasta".format(os.path.splitext(os.path.basename(alignment))[0] ) )
		speciesnames = args.taxa_names if args.taxa_names else [os.path.basename(newspeciesfile) for f in newtaxafiles]
		collect_sequences(nt_unaligned, alignment, seqids_to_add, seqdict, args.length, speciesnames, args.max_hits, args.supermatrix, args.no_trim)
		if args.no_trim: # use original method, which allows gaps from new sequences
			nt_aligned = run_mafft(ALIGNER, nt_unaligned, errorlog)
		else: # use --keeplength in mafft
			nt_aligned = run_mafft_addlong(ALIGNER, alignment, nt_unaligned, errorlog)

		# generate supermatrix from alignments
		newaligned = AlignIO.read(nt_aligned, "fasta")
		alignment_length = newaligned.get_alignment_length()
		partitionlist.append("{}:{}".format(runningsum+1, runningsum+alignment_length) )
		runningsum += alignment_length
		if supermatrix is None:
			supermatrix = newaligned
		else:
			supermatrix += newaligned
		nt_tree = run_tree(FASTTREEMP, nt_aligned, errorlog)

	### BUILD SUPERMATRIX
	if args.supermatrix:
		AlignIO.write(supermatrix, args.supermatrix, "fasta")
		print >> errorlog, "# Supermatrix written to {}".format(args.supermatrix), time.asctime()
		with open("{}.partition.txt".format(args.supermatrix),'w') as pf:
			print >> pf, ",".join(partitionlist)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout, sys.stderr)

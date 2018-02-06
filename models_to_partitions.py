#!/usr/bin/env python
#
# models_to_partitions.py  created 2017-10-16

'''models_to_partitions.py  last modified 2017-11-28
  convert supermatrix model blocks to simple gene partition pairs
  assumes format of model blocks is:
  LG, p1 = 1-406, 213446-213686

  and converted to comma/colon in numerical order:
  1:406,407:623,...

models_to_partitions.py partition2.txt > partitions_short.txt
'''

import sys

if len(sys.argv)<2:
	print >> sys.stderr, __doc__
else:
	partsum = 0
	partitionlist = [] # empty list, tuples will be added then sorted
	for line in open(sys.argv[1],'r'):
		line = line.strip()
		if line:
			partsplits = line.split('=')[1].replace(";","").split(',')
			for part in partsplits:
				part = part.strip()
				partsublist = [int(p) for p in part.split('-')]
				partitionlist.append(partsublist)
				partsum += partsublist[1] - partsublist[0] + 1
	partitionlist.sort()
	# last value should be equal to sum of parts
	lastvalue = partitionlist[-1][1]
	if lastvalue!=partsum:
		print >> sys.stderr, "WARNING: SUM {} DOES NOT EQUAL LAST PARTITION {}".format(partsum, lastvalue)
	print >> sys.stdout, ",".join( [ "{}:{}".format(*p) for p in partitionlist ] )

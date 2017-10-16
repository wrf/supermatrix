#!/usr/bin/env python
#
# models_to_partitions.py  created 2017-10-16

'''models_to_partitions.py  last modified 2017-10-16
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
	partitionlist = [] # empty list, tuples will be added then sorted
	for line in open(sys.argv[1],'r'):
		line = line.strip()
		if line:
			partsplits = line.split('=')[1].split(',')
			for part in partsplits:
				part = part.strip()
				partsublist = [int(p) for p in part.split('-')]
				partitionlist.append(partsublist)
	partitionlist.sort()
	print >> sys.stdout, ",".join( [ "{}:{}".format(*p) for p in partitionlist ] )

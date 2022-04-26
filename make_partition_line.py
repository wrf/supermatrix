#!/usr/bin/env python
#
# make_partition_line.py
# python3 update  2022-04-26

'''make_partition_line.py  last modified 2022-04-26

    make fasta format line to add to an alignment, 
    indicating where partitions start and stop for easy viewing
    appearing like:

>Partitions
-------------XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX------------------XXXXXXXXXXXXXX

make_partition_line.py partitions.txt > partition_line.fasta
'''

import sys

if len(sys.argv) < 2:
	sys.exit( __doc__ )
else:
	partitions = [] # list of tuples of intervals
	partition_file = sys.argv[1]

	sys.stderr.write( "# reading partitions from {}\n".format( partition_file ) )
	for line in open(partition_file,'r'):
		line = line.strip()
		if line:
			blocks = line.split(",") # split "1:136,137:301,..." into ['1:136', '137:301',...]
			for block in blocks:
				alignindex = tuple( int(i) for i in block.split(":") ) # split '1:136' into ( 1,136 )
				partitions.append(alignindex)
	partition_count = len(partitions)
	sys.stderr.write( "# found {} partitions\n".format( partition_count ) )
	fastastring = ""
	for i,part in enumerate(partitions):
		if partition_count < 11:
			symbol = str(i)
		else:
			if i % 2: # meaning odd
				symbol = "X"
			else:
				symbol = "-"
		partlength = part[1] - part[0] + 1
		fastastring += partlength * symbol

	sys.stdout.write( ">Partitions_{}\n{}\n".format(partition_count, fastastring) )
#

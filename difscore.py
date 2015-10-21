#!/usr/bin/python
#Emerson Glassey
#7-28-2012

#------------------------------------------------------------------------------#
#Difscore algorithm adapted from pluskal isotope pattern matching.
#------------------------------------------------------------------------------#

import numpy
import sys

print """An implementation of consensus clustering. kmeans clustering is used
for values of k between 1 and n (number of values to be clustered in
dataset). If two things cluster together, their value in the nxn matrix
is increased by one. The larger the value, the closer the things.
To use:
declare object:			obj = kmeans(cdtfile)
preform clustering:		obj.consensus()
NOTE: The data matrix can be reached at nx
"""

def readcdt(cdtfile):
	lines = []
	labels =[]
	with open(cdtfile) as f:
		for line in f:
			if "GENE" in line and "nan" not in line:
				lin = line.strip().split("\t")
				labels.append(lin[2])
				lines.append(map(float, lin[4:]))
	return labels, numpy.array(lines, dtype=numpy.float32)

def writecsv(labels, nx):
	with open(sys.argv[1][:-4] + "_nxn_difscore.csv", 'w') as f:
		f.write("NAME\t" + "\t".join(labels) + "\n")
		for ind1 in range(len(labels)):
			f.write(labels[ind1] + "\t" + "\t".join(map(str, nx[ind1])) + "\n")

def difscore(x, y):
	if not len(x) == len(y):
		print "lengths not same"
		return None
	dif = 1.0
	for xi, yi in zip(x, y):
		dif *= 1.0 - float(xi - yi) ** 2
	return dif
	
def main(labels, data):
	nxn = numpy.zeros((len(labels), len(labels)))
	for ind1 in range(len(labels)):
		for ind2 in range(len(labels)):
			nxn[ind1][ind2] = difscore(data[ind1], data[ind2])
		print ind1, " of ", len(labels)
	return nxn

if __name__ == "__main__":

	labels, data = readcdt(sys.argv[1])
	nx = main(labels, data)
	writecsv(labels, nx)
#!/usr/bin/python
#Emerson Glassey
#7-28-2012

#------------------------------------------------------------------------------#
#Though the name of this module is kmeans, this is an implementation of
# consensus clustering using the kmeans clustering algorithm with euclidean
# distance. The idea is:
# As you increase the k value (number of clusters to divide dataset), if things
# continue to cluster together then they must be very similar. As such, the
# supplied dataset is subject to k values between 1 and n (the number of things
# to cluster in the dataset). Every time to things (IDruns) cluster together,
# a 1 is added to their shared value in the nxn matrix. The larger the value,
# the closer the IDruns. 
#------------------------------------------------------------------------------#

import numpy
import Pycluster, pp, sys, logging, time


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


def clusters(labels, data, k):
	kclus = Pycluster.kcluster(data, k, npass=1)[0]
	nx = numpy.zeros((len(labels), len(labels)), dtype=numpy.float32)
	for ind1 in range(len(labels)):
		for ind2 in range(len(labels)):
			if kclus[ind1] == kclus[ind2]:
				nx[ind1][ind2] = 1
	print k, " of ", len(labels)
	return nx

def writecsv(labels, nx):
	with open(sys.argv[1][:-4] + "_nxn_kmean.csv", 'w') as f:
		f.write("NAME\t" + "\t".join(labels) + "\n")
		for ind1 in range(len(labels)):
			f.write(labels[ind1] + "\t" + "\t".join(map(str, nx[ind1])) + "\n")

def consensus(labels, data):
	print labels
	print data
	jobs = []
	nx = numpy.zeros((len(labels),len(labels)), dtype=numpy.float32)
	for run in range(1, len(labels)):
		jobs.append(job_server.submit(clusters, (labels, data, run,), (), ("Pycluster", "numpy",)))
	num = 0
	while len(jobs) > 0:
		nx = nx + jobs.pop(0)()
		if num%100 == 10:
			print nx
		if num%100 == 20:
			print nx[len(labels)/2]
		num += 1
	numpy.round(nx, 3)
	return nx / (len(labels) - 1)

def hosts(hostfile):
	ppservers = []
	with open(hostfile) as f:
		line = f.readline().strip()
		while "master" not in line:
			line = f.readline().strip()
		try:
			ppservers.append(line.split()[1] + ":35000")
		except IndexError:
			pass
		while True:
			try:
				line = f.readline().strip()
				ppservers.append(line.split()[1] + ":35000")
			except IndexError:
				break
	if len(ppservers) == 0:
		ppservers = ()
	return tuple(ppservers)

if __name__ == "__main__":
	ppservers = ()
#	ppservers = hosts("/etc/hosts")
	try:
		ppservers = tuple(sys.argv[2].split(","))
		ppservers = tuple([serve + ":35000" for serve in ppservers])
	except IndexError:
		ppservers = ()
	job_server = pp.Server(ppservers=ppservers, secret="acetone")
	print ppservers, job_server.get_ncpus()
	labels, data = readcdt(sys.argv[1])
	nx = consensus(labels, data)
	writecsv(labels, nx)
# Written by Emerson Glassey
# on ~ 7-26-2012
#
# Edited by Emerson Glassey
# 7-30-2012
# to include a negative/anticluster function
#------------------------------------------------------------------------------#
# Python function for displaying information about clusters displayed in an nxn
# matrix that is or isn't clustered (but is more visually appealing if it is). 
# The nxn matrix can be generated by the pearson.py or kmeans.py programs in the
# form of a .csv file.
#------------------------------------------------------------------------------#

import numpy as np


class nbyn():
	"""Python class that defines an object with a cdtfile from cluster.
	To define: 			obj = cluster(cdtfile)
	To cluster:			obj.cluster(idrun, min, max)
	To access cluster list:		obj.cruns
	NOTE: obj.cruns is reset every time the cluster function is called.
	"""
	def __init__(self, cdtfile):
		#Save the connection to the file
		self.f = open(cdtfile)
		self.labels = self.f.readline().strip().split("\t")
		self.cdt = cdtfile
	
	def anticluster(self, idrun, kmax = 0, fraction = 5, limit = None):
		self.aruns = []
		self.f.seek(0)
		line = self.f.readline().strip().split("\t")
		while idrun not in line[0]:
			line = self.f.readline().strip().split("\t")
		if kmax == 0:
			kmax = self.mean - 0.1*self.std
		for pind in range(1, len(line)):
			if float(line[pind]) <= kmax and self.labels[pind][:-4] not in set([id[:-4] for id in self.cruns]):
				self.aruns.append(self.labels[pind])
		if limit:
			return self.aruns[::int(len(self.aruns)/limit)]
		return self.aruns[::fraction]
				

	def cluster(self, idrun, min_tolerance=0.8, max_tolerance=0.5, min_factor=None, max_factor=None, view=False):
		#Always start at the beginning of the file (but skip first line)
		self.f.seek(0)
		self.f.readline()
		if max_tolerance != "":
			self.limitruns = []
			self.cruns = []
		#read until the idrun row is found
		line = self.f.readline().strip()
		while idrun not in line:
			line = self.f.readline().strip()
		line = line.split("\t")
		if line == [""]:
			print "IDrun not present"
			return
		#Plot data if wanted, comment out if only interested in clusters
		x = []
		y = []
		ind = 2
		for pear in map(float, line[2:]):
			x.append(ind)
			y.append(pear)
			ind += 2
		if min_factor and max_factor:
			mean, std = (np.mean(np.array(y)), np.std(np.array(y)))
			min_tolerance = mean + std * min_factor
			if max_tolerance != "":
				max_tolerance = mean + std * max_factor
				self.mean, self.std = mean, std
		
		if max_tolerance != "" and view:
			import matplotlib
			matplotlib.use("pdf")
			import matplotlib.mlab as mlab
			import matplotlib.pyplot as plt
			plt.scatter(np.array(x), np.array(y))
			plt.axis([-10, max(x)+30, -0.025, 1.1])
#			plt.text(len(self.labels) / 3, 0.7, "mean: {0}, std: {1}".format(self.mean, self.std))
			plt.savefig(idrun + ".pdf")
			plt.close()
		# Iterate through all the idruns (columns) in the row
		for pind in range(2, len(line)):
			#If it's the first iteration, find the limit runs and the runs to iterate over
			if max_tolerance != "":
				if float(line[pind]) > max_tolerance:
					self.limitruns.append(self.labels[pind])
				if float(line[pind]) > min_tolerance:
					self.cruns.append(self.labels[pind])
					self.cluster(self.labels[pind], max_tolerance="", min_factor=min_factor, max_factor=max_factor)
			#If it's not the first iteration, find more runs to iterate over
			#Only iterate over runs that haven't been found and are within the limit
			else:
				if float(line[pind]) > min_tolerance and self.labels[pind] not in self.cruns and self.labels[pind] in self.limitruns:
					self.cruns.append(self.labels[pind])
					self.cluster(self.labels[pind], max_tolerance="", min_factor=min_factor, max_factor=max_factor)
		return 
		
	def score(self, sample1, sample2):
		self.f.seek(0)
		self.f.readline()
		line = self.f.readline().strip()
		while sample1 not in line:
			line = self.f.readline().strip()
		line = line.split("\t")
		if line == [""]:
			print "IDrun not present"
			return
		for ind, label in enumerate(self.labels):
			if sample2 in label:
				return float(line[ind])

	def classes(self):
		self.f.seek(0)
		classes = dict()
		for name in self.labels[1:]:
			clas = name.split(" - ")[-1].split("_")[0].strip()
			classes.setdefault(clas, [])
			classes[clas].append(name)
		return classes
	
	def reorder(self):
		classes = classes()
		new = np.empty((len(self.labels) - 1, len(self.labels) - 1))
		labels = [name for clas in classes.values() for name in clas if len(clas) > 1]
		indices = []
		for label in self.labels[1:]:
			try:
				indices.append(labels.index(label))
			except ValueError:
				continue
		self.f.readline()
		for line in self.f:
			line = line.strip().split("\t")
			try:
				ind = labels.index(line[0])
			except ValueError:
				continue
			for val in range(len(line[1:])):
				try:
					new[ind][indices[val]] = line[val + 1]
				except IndexError:
					continue
		with open(self.cdt[:-4] + "_reorder.csv", 'w') as g:
			g.write("NAME\t" + "\t".join(labels) + "\n")
			for ind1 in range(len(labels)):
				g.write(labels[ind1] + "\t" + "\t".join(map(str, new[ind1])) + "\n")


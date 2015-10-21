

"""Module contains classes defining objects for parsing and indexing cytological profiling data.
"""

from scipy.stats import pearsonr
import pp
import numpy
from progressbar import ProgressBar, Percentage, Bar, RotatingMarker, ETA, FileTransferSpeed
import copy

import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt

class feature:
	"""Class defines a cytological profiling feature. """
	def __init__(self, name):
		self.name = name
		self.val = dict()

	def __getitem__(self, run):
		return self.val[run]
	
	def __setitem__(self, run, value):
		self.val[run] = value
	
	def __delitem__(self, run):
		del self.val[run]
		
	def __repr__(self):
		return self.name
		
	def values(self):
		return self.val.values()
		
	def keys(self):
		return self.val.keys()
		
	def items(self):
		return self.val.items()
		
	def merge(self, feature):
		"""Function takes another feature object and merges it into this feature object if the
		names match."""
		if feature.name == self.name:
			for run, val in feature.items():
				self.val[run] = val

class fingerprint:
	"""Class defines a cytological fingerprint for one compound/well."""
	def __init__(self, name):
		self.name = name
		self.finger = dict()

	def __getitem__(self, feat):
		return self.finger[feat]
	
	def __setitem__(self, feat, value):
		self.finger[feat] = value
		
	def __delitem__(self, feat):
		del self.finger[feat]
		
	def __repr__(self):
		return self.name
		
	def values(self):
		return self.finger.values()
	
	def keys(self):
		return self.finger.keys()
		
	def items(self):
		return self.finger.items()
		
	def activity_score(self):
		return sum([val ** 2 for val in self.finger.values()])

class cp:
	"""Class defines a cytological profiling heatmap. """
	def __init__(self, cpfile):
		self.map = dict()
		if cpfile:
			feats, runs = self.parse_tab(open(cpfile))
			for run in runs:
				self.map[run.name] = run
			for feat in feats:
				self.map[feat.name] = feat
			self.pearson()

	def __getitem__(self, key):
		return self.map[key]
	
	def __setitem__(self, key, value):
		self.map[key] = value
	
	def __delitem__(self, key):
		sys.stderr.write("Cannot edit single fingerprints.\n")
		
	def keys(self):
		return self.fingerprints()
	
	def fingerprints(self):
		return [run for run in self.map.values() if isinstance(run, fingerprint)]
	
	def features(self):
		return [feat for feat in self.map.values() if isinstance(feat, feature)]

	def parse_tab(self, tab):
		"""Function takes an open tab file of CP data and inputs it into the internal
		dictionaries."""
		runs = []
		iter_num = int(tab.name.rsplit("_", 1)[1].split(".")[0])
		for line in tab:
			if "FEATURES" in line.strip().upper():
				features = [feature(feat) for feat in line.strip().upper().split("\t")[1:]]
			else:
				finger = line.strip().split("\t")
				runs.append(fingerprint("{}_{}".format(finger[0].split("_")[0], iter_num + 1)))
				for feat, value in zip(features, finger[1:]):
					runs[-1][feat.name] = feat[runs[-1].name] = float(value)
		return features, runs
		
	def remove_activity(self, score):
		"""Function accepts an activity score and returns a cp object with only fingerprints with
		an activity score greater than that given.
		"""
		runs = set([str(run) for run in self.fingerprints() if run.activity_score() > score])
		feats = [str(feat) for feat in self.features()]
		subcp = cp(None)
		sub = dict()
		for run in runs:
			sub[run] = fingerprint(run)
		for feat in feats:
			sub[feat] = feature(feat)
			for run in runs:
				sub[run][feat] = sub[feat][run] = self[run][feat]
		subcp.map = sub
		subcp.pearson()
		return subcp
		
	def restrict_runs(self, runs):
		"""Function accepts a list of fingerprint names and returns a new cp object with only
		those fingerprints."""
		runs = set(runs)
		feats = [str(feat) for feat in self.features()]
		subcp = cp(None)
		sub = dict()
		for run in runs:
			sub[run] = fingerprint(run)
		for feat in feats:
			sub[feat] = feature(feat)
			for run in runs:
				sub[run][feat] = sub[feat][run] = self[run][feat]
		subcp.map = sub
		subcp.pearson()
		return subcp
		
	def write_tab(self, tab):
		"""Function takes an open file object and writes the cp data matrix to that file as a 
		tab delimited text file.
		"""
		f = open(tab, 'w')
		runs = self.fingerprints()
		feats = self.features()
		f.write("Features\t{}\n".format("\t".join([str(feat) for feat in feats])))
		for run in runs:
			f.write("{}\t{}\n".format(str(run), "\t".join([str(run[str(feat)]) for feat in feats])))
		f.close()
		
	def write_specific(self, tab, specific=[]):
		f = open(tab, 'w')
		runs = [run for run in self.fingerprints() if run.name.split("_")[0] in specific or run.name.split("_")[0][:-1] in specific]
		feats = self.features()
		f.write("Features\t{}\n".format("\t".join([str(feat) for feat in feats])))
		for run in runs:
			f.write("{}\t{}\n".format(str(run), "\t".join([str(run[str(feat)]) for feat in feats])))
		f.close()
		
	def get_pearson(self, id1, id2):
		return self.nxn[(self.labels.index(id1), self.labels.index(id2))]
		
		
	def pearson(self):
		"""Function returns an nxn heatmap of pearson correlations with idruns as rows and columns
		as columns. The tuple that is returned is (nxn, labels)
		"""
		idruns = self.fingerprints()
		self.labels = [id.name for id in idruns]
		
		widgets = ['Pearson Correlations: ', Percentage(), ' ', Bar(marker=RotatingMarker()), ' ',
				ETA(), ' ', FileTransferSpeed()]
		pbar = ProgressBar(widgets=widgets, maxval=len(idruns)).start()
		
		self.nxn = numpy.zeros((len(idruns), len(idruns)))
		for ind1, id1 in enumerate(idruns):
			for ind2, id2 in enumerate(idruns):
				self.nxn[(ind1,ind2)] = pearsonr(id1.values(), id2.values())[0]
			pbar.update(ind1)
		pbar.finish()
		
	def cos(self, x, y):
		if not len(x) == len(y):
			print "lengths not same"
			return None
		magx = sum([xi**2 for xi in x]) ** (1.0/2.0)
		magy = sum([yi**2 for yi in y]) ** (1.0/2.0)
		
		dot = 0.0
		for xi, yi in zip(x, y):
			dot += xi * yi
		dot /= float(magx * magy)
		return dot
		
	def export_nxn(self):
		"""Function returns an nxn heatmap of retention times with idruns as rows and baskets
		as columns. The tuple that is returned is (nxn, row labels, column labels)
		"""
		runs = self.fingerprints()
		params = sorted(self.features(), key=lambda x: str(x))
		
		nxn = numpy.zeros((len(runs), len(params)))
		for rind, run in enumerate(runs):
			for pind, param in enumerate(params):
				if str(run) in self.map:
					if str(param) in self.map[str(run)].finger:
						nxn[(rind, pind)] = run[str(param)]
		return nxn, [str(run) for run in runs], [str(param) for param in params]
		
	def coscore(self):
		"""Function returns an nxn heatmap of pearson correlations with idruns as rows and columns
		as columns. The tuple that is returned is (nxn, labels)
		"""
		idruns = self.fingerprints()
		if not self.labels:
			self.labels = [id.name for id in idruns]
		
		widgets = ['CoScore: ', Percentage(), ' ', Bar(marker=RotatingMarker()), ' ',
				ETA(), ' ', FileTransferSpeed()]
		pbar = ProgressBar(widgets=widgets, maxval=len(idruns)).start()
		
		self.conxn = numpy.zeros((len(idruns), len(idruns)))
		for ind1, id1 in enumerate(idruns):
			for ind2, id2 in enumerate(idruns):
				self.conxn[(ind1,ind2)] = self.cos(id1.values(), id2.values())
			pbar.update(ind1)
		pbar.finish()
	
	def add_tab(self, tab):
		"""Function takes an open tab file of CP data and appends it to the internal
		dictionaries."""
		tab = open(tab)
		feats, runs = self.parse_tab(tab)
		for run in runs:
			self.map[run.name] = run
		for feat in feats:
			if str(feat) in self.map:
				self.map[str(feat)].merge(feat)
			else:
				self.map[str(feat)] = feat
				
	def anticluster(self, idrun, pmax=0, fraction=5, limit=None):
		if idrun not in self.labels:
			print idrun, ' not in cp map'
			return []
		aruns = []
		for ind, pear in enumerate(self.nxn[self.labels.index(idrun)]):
			if pear <= pmax and self.labels[ind][:-1] not in [id[:-1] for id in self.cruns]:
				aruns.append(self.labels[ind])
		if limit:
			return aruns[::int(len(aruns)/limit)]
		return aruns[::fraction]
				

	def cluster(self, idrun, min_tolerance=0.7, max_tolerance=0.4):
		if idrun not in self.labels:
			print idrun, ' not in cp map'
			return []
		if max_tolerance != "":
			self.limitruns = []
			self.cruns = []
		for ind, pear in enumerate(self.nxn[self.labels.index(idrun)]):
			#If it's the first iteration, find the limit runs and the runs to iterate over
			if max_tolerance != "":
				if pear > max_tolerance:
					self.limitruns.append(self.labels[ind])
				if pear > min_tolerance:
					self.cruns.append(self.labels[ind])
					self.cluster(self.labels[ind], min_tolerance, "")
			#If it's not the first iteration, find more runs to iterate over
			#Only iterate over runs that haven't been found and are within the limit
			else:
				if pear > min_tolerance and self.labels[ind] not in self.cruns and self.labels[ind] in self.limitruns:
					self.cruns.append(self.labels[ind])
					self.cluster(self.labels[ind], min_tolerance, "")
		return self.cruns
		
	def plot_parameters(self, param1x, param2y):
		x = []
		y = []
		for run in self.fingerprints():
			x.append(run[param1x])
			y.append(run[param2y])
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.scatter(x, y)
		ax.set_xlim(-1,1)
		ax.set_ylim(-1,1)
		plt.savefig("{}_{}.pdf".format(param1x, param2y))
		plt.close()
		
	def plot_activity(self, name):
		fingers = self.fingerprints()
		x = range(len(fingers))
		y = sorted([finger.activity_score() for finger in fingers])
		f, (ax1, ax2) = plt.subplots(1, 2, sharey=False)
		ax1.scatter(x, y)
		ax2.hist(y, 30)
		plt.savefig(name + "_activity.pdf")
		plt.close()
		
	def plot_parameter(self, param):
		y = sorted(self[param].values())
		x = range(len(x))
		f, (ax1, ax2) = plt.subplots(1, 2, sharey=False)
		ax1.scatter(x, y)
		ax2.hist(y, 30)
		plt.savefig(param + "_plot.pdf")
		plt.close()
		for pair in sorted(self[param].items(), key=lambda p:p[1]):
			print "{}, {:4f}".format(*pair)
import numpy as np
import networkx as nx
"""
For reading and comparing heatmaps of yeast data. Data can be appended together
and exported as a network or as a heatmap.

"""


class synleth:
	def __init__(self, heatmap):		
		self.data, self.strains, self.prefractions = self.parse_tab(heatmap)
				
	def __getitem__(self, key):
		return self.map[key]
	
	def __setitem__(self, key, val):
		self.map[key] = val
		
	def __delitem__(self, key):
		del self.map[key]
				
		
	def parse_tab(self, tabfile):
		strains = set()
		prefractions = set()
		data = dict()
		with open(tabfile, "U") as f:
			strains.update(set(f.readline().strip().split("\t")[1:]))
			for line in f:
				line = line.strip().split("\t")
				pref = line.pop(0)
				prefractions.add(pref)
				data.setdefault(pref, dict())
				
				for strn, val in zip(strains, line):
					data.setdefault(strn, dict())
					data[pref][strn] = data[strn][pref] = int(val)
		return data, strains, prefractions
		
	def append_data(self, data, strains, prefractions):
		self.strains.intersection_update(strains)
		self.prefractions.update(prefractions)
		
		new_data = dict()
		
		for strn in self.strains:
			new_data.setdefault(strn, dict())
			for pref in self.prefractions:
				new_data.setdefault(pref, dict())
				if pref in data:
					new_data[pref][strn] = new_data[strn][pref] = data[pref][strn]
				else:
					new_data[pref][strn] = new_data[strn][pref] = self.data[pref][strn]
		self.data = new_data
		
	def remove_zero_fingerprints(self):
		for pref in list(self.prefractions):
			fingerprint = self.data[pref].values()
			if fingerprint.count(0) == len(fingerprint):
				self.prefractions.discard(pref)
				del self.data[pref]
				for strn in self.strains:
					del self.data[strn][pref]

	def write_heatmap(self, filename):
		"""Function accepts a filename, a data dictionary and lists of prefractions and
		strains, and writes a file with prefraction rows and strain columns with the value
		in the corresponding row/column.
		"""
		with open(filename, 'w') as f:
			f.write("Prefractions\t{}\n".format("\t".join(self.strains)))
			for pref in self.prefractions:
				f.write("{}\t{}\n".format(pref, "\t".join([str(self.data[pref][strn]) for strn in self.strains])))
		return
		
	def write_network(self, filename, mismatches = 1):
		"""Function accepts a output network filename, a data dictionary, and a list of prefractions
		and outputs the network of prefraction nodes, with connections if their yeast fingerprint
		contains fewer than the percentage of mismatches specified in mismatches (default is half).
		"""
		g = nx.Graph()
		
		for pref in self.prefractions:
			g.add_node(pref, activity = sum(self.data[pref].values()))
			for pref2 in self.prefractions:
				if pref == pref2:
					continue
				if self.data[pref].values() == self.data[pref2].values() == [0] * len(self.data[pref].values()):
					continue
				mismatch = 1
				for pair in zip(self.data[pref].values(), self.data[pref2].values()):
					if pair.count(pair[0]) == 2:
						continue
					mismatch += 1
				
				if mismatch <= mismatches+1:
					g.add_edge(pref, pref2, weight = 1.0 / mismatch)
		nx.write_dot(g, filename)
		
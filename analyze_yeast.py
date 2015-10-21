
import numpy
import dbtables
from sqlalchemy.orm import sessionmaker, subqueryload_all
from progressbar import ProgressBar, Percentage, Bar, RotatingMarker, ETA, FileTransferSpeed
import yeast
import yeast_mzbin as mb
import btools
import sys

import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt

from scipy.stats import norm

import itertools
import argparse

import networkx as nx

cm = {'red':((0.0, 1.0, 1.0), (0.0000001, 1.0, 0.0), (1.0, 1.0, 1.0)),\
	'green':((0.0, 1.0, 1.0), (0.0000001, 0.0, 0.0), (1.0, 0.0, 0.0)),\
	'blue':((0.0, 1.0, 1.0), (0.0000001, 1.0, 1.0), (1.0, 0.0, 0.0))}
mycm = matplotlib.colors.LinearSegmentedColormap('mycm', cm)

class yeastbin():
	
	def __init__(self, db, tabfile, quant=True, remove=True, bqs=['00110', '00100', '10000', '10100', '10110', '11000', '11100', '11110','00111', '00101', '10001', '10101', '10111', '11001', '11101', '11111']):
		bqs = set(bqs)
		self.db = db
		self.engine = dbtables.connect(db)
		self.Session = sessionmaker(bind=self.engine)
		self.session = self.Session()
		
		self.synleth = self.get_synleth(tabfile)
		self.ad_heat = self.get_bins([str(pref) for pref in self.synleth.prefs], bqs)
		
		if remove:
			self.remove_systematic(0.1)
			self.remove_underrepresent()
		
		if quant:
			self.quant()
		
		
	def get_synleth(self, tabfile):
		"Function returns the cp object generated from the passed in tab file."
		return yeast.synleth(tabfile)

	def get_idruns(self, idruns, bqs):
		"""Function accepts a list of idruns of the format \d{4}[A-F]_\d{2} (4 digits followed by a 
		letter A through F *underscore* minute number). Returns a list of IDrun objects with loaded
		adducts and isotopes.
		"""
		return self.session.query(dbtables.Adduct).join(dbtables.IDrun).filter(dbtables.IDrun.name.\
				in_(idruns)).filter(dbtables.Adduct.bq.in_(bqs)).options(\
				subqueryload_all('isotopes')).all()
	
	def get_bins(self, idruns, bqs):
		"""Function takes a list of idruns of the format \d{4}[A-F] (4 digits followed by a
		letter A through F). Returns a mzbin heatmap object generated from the mass spec data
		in the loaded database.
		"""
		runs = self.get_idruns(idruns, bqs)

		return mb.heatmap(*[ad for ad in runs if ad.isotopes[0].rt > 0.35])
		
	def quant(self):		
		baskets = self.ad_heat.baskets()
		
		for ind, basket in enumerate(baskets):
			basket.subprint = basket.subfingerprint(*[self.synleth[pref].fingerprint for pref in basket.keys()])
			#for pref in basket.keys():
				#basket.fingerprints[pref] = basket.fingerprint_combos(self.synleth[pref].fingerprint)
				#basket.tally_fingerprints()
				
		
	def remove_systematic(self, max_length=0.25):
		"""Function takes a number represented the maximum fraction of idruns in the imported cp
		file(s) that can contain a compound/adduct in a basket (maximum length of a basket) before
		the compound/adduct is considered a systematic contaminant. Default is 0.25. Baskets with
		sizes over this are removed.
		"""
		idrun_count = len(self.synleth.prefs)
		n = 0
		initial = len(self.ad_heat.baskets())
		for bask in self.ad_heat.baskets():
			if len(bask.keys()) > float(idrun_count) * max_length:
				self.ad_heat.remove_basket(bask.name)
				n += 1
		print "Removed by Length Ceiling:", n, "of", initial
		
	def remove_low_fraction(self, min_effect=0.5):
		"""Function accepts a number representing the minimum fraction of features in the imported
		cp file(s) that must be affected (outside the standard deviation) for a compound/adduct to
		be considered 'active'. Inactive compounds/adducts are removed from the heatmap.
		"""
		n = 0
		initial = len(self.ad_heat.baskets())
		for basket in self.ad_heat.baskets():
			if basket.fraction_covered < min_effect:
				self.ad_heat.remove_basket(basket.name)
				n += 1
		print "Removed by Inactivity:", n, "of", initial
		
	def remove_underrepresent(self, min_num=1, min_percent=None):
		"""Function accepts a number representing the minimum number of peaks that must be in a
		basket for it to not be classified as underrepresented and removed as noise. Defualt is 1.
		"""
		if not min_percent == None:
			mp = min_percent * len(self.ad_heat.idruns())
			if mp > min_num:
				min_num = mp
		initial = len(self.ad_heat.baskets())
		
		self.ad_heat.remove_undercount(min_num)
		
		print "Removed by Length Floor:", initial - len(self.ad_heat.baskets()), "of", initial
		
	def remove_low_count(self, min_count=0.1):
		"""Function accepts a number representing the minimum specificity for the basket to not be
		considered noise. The specificity is determined by the cluster_score, a high cluster_score
		indicates a high average pearson correlation between all the idruns in the basket, meaning
		the compound/basket is fairly specific.
		"""
		n = 0
		initial = len(self.ad_heat.baskets())
		for basket in self.ad_heat.baskets():
			if basket.inverse_count < min_count:
				self.ad_heat.remove_basket(basket.name)
				n += 1
		print "Removed by Cluster Score:", n, "of", initial
		
	def remove_nonspecific(self, min_specificity=0.1):
		"""Function accepts a number representing the minimum specificity for the basket to not be
		considered noise. The specificity is determined by the cluster_score, a high cluster_score
		indicates a high average pearson correlation between all the idruns in the basket, meaning
		the compound/basket is fairly specific.
		"""
		n = 0
		initial = len(self.ad_heat.baskets())
		for basket in self.ad_heat.baskets():
			if basket.cluster_score < min_specificity:
				self.ad_heat.remove_basket(basket.name)
				n += 1
		print "Removed by Cluster Score:", n, "of", initial
		
	def remove_nonactive(self, min_effect=5):
		"""Function accepts a number representing the minimum fraction of features in the imported
		cp file(s) that must be affected (outside the standard deviation) for a compound/adduct to
		be considered 'active'. Inactive compounds/adducts are removed from the heatmap.
		"""
		n = 0
		initial = len(self.ad_heat.baskets())
		for basket in self.ad_heat.baskets():
			if basket.activity < min_effect:
				self.ad_heat.remove_basket(basket.name)
				n += 1
		print "Removed by Inactivity:", n, "of", initial

	def bask_network(self):
		"""Function exports a network of fingerprints as nodes and edges as
		compounds.
		"""
		h = [nx.Graph(), nx.Graph(), nx.Graph()]
		edges = dict()
		for id in self.synleth.prefs:
			h[0].add_node(id.name, weight=id.activity_score())
			h[1].add_node(id.name, weight=id.activity_score())
			h[2].add_node(id.name, weight=id.activity_score())
		for bask in self.ad_heat.baskets():
			for pair in itertools.combinations(sorted(bask.keys()), 2):
				edges.setdefault((pair[0][:-1], pair[1][:-1]), [])
				edges[(pair[0][:-1], pair[1][:-1])].append(bask.activity * bask.cluster_score)
		for edge in edges:
			h[0].add_edge(*edge, weight=sum(edges[edge]))
			h[1].add_edge(*edge, weight=numpy.mean(edges[edge]))
			h[2].add_edge(*edge, weight=max(edges[edge]))
		return h
		
	def co_express_network(self):
		"""Function generates a network of features, connecting them if they appear in a
		prefraction together.
		"""
		h = nx.Graph()
		edges = dict()
		for fin in self.ad_heat.baskets():
			h.add_node(fin, weight=fin.activity * fin.cluster_score, mass=fin.m, rt=fin.rt, bid=fin.name, num=len(fin.keys()))
		widgets = ['VectorMove: ', Percentage(), ' ', Bar(marker=RotatingMarker()), ' ', ETA(), ' ', FileTransferSpeed()]
		pbar = ProgressBar(widgets=widgets, maxval=len(self.ad_heat.idruns())).start()
		for num, idrun in enumerate(self.ad_heat.idruns()):
			pbar.update(num)
			for pair in itertools.combinations(sorted(idrun.keys()), 2):
				edges.setdefault((self.ad_heat[pair[0]], self.ad_heat[pair[1]]), [])
				edges[(self.ad_heat[pair[0]], self.ad_heat[pair[1]])].append(1)
		pbar.finish()
		for edge in edges:
			h.add_edge(*edge, weight=len(edges[edge]))
		return h
		
	def pref_mz_network(self):
		h = nx.Graph()
		for pref in self.synleth.prefs:
			h.add_node(str(pref), activity = pref.activity_score(), weight = pref.activity_score(), t=1)
		for mass in g.ad_heat.baskets():
			h.add_node(mass.m, mz=mass.m, rt=mass.rt, activity=mass.activity, cluster=mnass.cluster_score, weight=mass.activity * 0.5, t=2)
			for pref in mass.peaks:
				h.add_edge(pref, mass.m, weight = mass.cluster_score)
		return h
			
        
	def export_plot(self, idrun, name=None, pmax=None, maxt=0.4, mint=0.7, text=1, title=True, min_count=None, min_percent=None):
		"""Function accepts an idrun, finds the cluster and anticluster for that idrun, and saves
		a pdf of the heatmap for that cluster nxn matrix. Name is the idrun_heat.pdf.
		"""
		submap = self.ad_heat.subheat([idrun])
		nxn, rows, cols = submap.export_nxn()
		if len(rows) == 0 or len(cols) == 0:
			print "empty heatmap"
			return
		
		if title:
			title = "pmax={};maxt={};mint={}".format(pmax, maxt, mint)
		else:
			title = ""
		name = idrun+"_heatmap" if not name else name
		self.plot(name, nxn, cols, rows, title, text)
		
	def plot(self, name, matrix, bins, rows, title="", text=4):
		"""plots a heatmap of the supplied matrix, bins, and rows. title and text set
		title and text size of the plot. Graph is saved as name.pdf
		"""
		plt.matshow(matrix, cmap=mycm, vmin=0, vmax=6)
		plt.xticks(range(len(bins)), [str(bn) for bn in bins], rotation=90, size=text)
		plt.yticks(range(len(rows)), rows, size=text)
		ax = plt.gca()
		for posi in ax.spines:
			ax.spines[posi].set_color('none')
		ax.tick_params(labelcolor='k', top='off', bottom='off', left='off', right='off')
		plt.title(title)
		plt.savefig(name + ".pdf", bbox_inches='tight')
		plt.close()
		return
		
	def export_activity(self, idrun, maxt=0.5, mint=0.65, size=10, aline=8, cline=0.1, alpha=0.75, figsize=(4, 4)):
		"""Exports activity plot to idrun_activity_plot.pdf. maxt and mint are for cluster calling.
		size is for text size in the plot. aline and cline are the green lines visually
		showing the cutoffs. alpha specifies the transparency of the spots and figsize
		is the plot dimensions (in, in). idrun must be in the cpbin heatmap.
		"""
		#submap of the specified idrun
		submap = self.ad_heat.subheat([idrun])
		if submap == None:
			print "Empty Heatmap!"
			return
		basks = submap.baskets()
		plt.rc('xtick', labelsize=10)
		plt.rc('ytick', labelsize=10)
		fig = plt.figure(figsize=figsize, dpi=300)
		ax = fig.add_subplot(111)
		ax.set_ylim(0, 20)
		ax.set_xlim(0, 1)
		ax.axhline(aline, color='g')
		ax.axvline(cline, color='g')
		ax.set_xlabel("Cluster Score", size=10)
		ax.set_ylabel("Activity Score", size=10)

		x = [bask.cluster_score for bask in basks]
		y = [bask.activity for bask in basks]
		c = [bask.rt for bask in basks]

		plt.scatter(x, y, s=size, c=c, cmap=mycm, vmin=0, vmax=6, alpha=alpha)
		plt.savefig(idrun + "_activity_plot.pdf", bbox_inches='tight')
		plt.close()
		return	  
				
	
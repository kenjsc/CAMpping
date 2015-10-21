"""Analyze module for CAMpping. This is the module for running CAMpping analyses.
Module can be run as a command line utility for vector compression, or can be run
interactively from within python to generate heatmaps, activity plots, etc. If
running interactively, you need to first generate a cpbin object with the relevant
CP file and mz database of interest, and then run analyses through the cpbin structure.
"""

import numpy
import dbtables
from sqlalchemy.orm import sessionmaker, subqueryload_all
from progressbar import ProgressBar, Percentage, Bar, RotatingMarker, ETA, FileTransferSpeed
import cp
import mzbin as mb
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

class cpbin():
	
	def __init__(self, db, tabfile, quant=True, remove=True, bqs=['00110', '00100', '10000', '10100', '10110', '11000', '11100', '11110','00111', '00101', '10001', '10101', '10111', '11001', '11101', '11111']):
		bqs = set(bqs)
		self.db = db
		self.engine = dbtables.connect(db)
		self.Session = sessionmaker(bind=self.engine)
		self.session = self.Session()
		
		self.cp = self.get_cp(tabfile)
		self.ad_heat = self.get_bins([finger.name for finger in self.cp.fingerprints()], bqs)
		
		if remove:
			self.remove_systematic(0.1)
			self.remove_underrepresent()
		
		if quant:
			self.quant()
		
		
	def get_cp(self, tabfile):
		"Function returns the cp object generated from the passed in tab file."
		return cp.cp(tabfile)

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
		
	def remove_systematic(self, max_length=0.25):
		"""Function takes a number represented the maximum fraction of idruns in the imported cp
		file(s) that can contain a compound/adduct in a basket (maximum length of a basket) before
		the compound/adduct is considered a systematic contaminant. Default is 0.25. Baskets with
		sizes over this are removed.
		"""
		idrun_count = len(self.cp.fingerprints())
		n = 0
		initial = len(self.ad_heat.baskets())
		for bask in self.ad_heat.baskets():
			if len(bask.keys()) > float(idrun_count) * max_length:
				self.ad_heat.remove_basket(bask.name)
				n += 1
		print "Removed by Length Ceiling:", n, "of", initial
		
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
		
	def remove_inactive_runs(self, min_effect=5):
	        """Function removes all idruns from the mzbin heatmap that have an activity score below the
	        specified value"""
		for id in self.cp.fingerprints():
			if id.activity_score() <= min_effect and id.name + "_" in self.ad_heat.map:
				self.ad_heat.remove_run_baskets(id.name + "_")
				
		
	def get_finger(self, basket, features):
		"""Function accepts a basket, and a list of feature objects from the cp module.
		Using scores in the features, the function returns a fingerprint representative of the basket
		idruns by constructing the average fingerprint (mean of each feature between the idruns).
		"""
		idrun_names = basket.keys()
		activity = 0
		for feat in features:
			in_scores = [feat[run] for run in idrun_names]
			mu = numpy.mean(in_scores)
			std = numpy.std(in_scores)
			basket.param_score(feat.name, (mu, std))
			activity += mu ** 2
		basket.activity = activity
		
	def get_cluster(self, basket, power=3):
	        """Function calculates the cluster score for the passed in basket raised to the specified power.
	        Default is 3."""
		n = 0
		cluster_score = 0
		for pair in itertools.combinations(basket.keys(), 2):
			cluster_score += self.cp.get_pearson(pair[0], pair[1]) ** power
			n += 1
		if n == 0:
			basket.cluster_score = 1
		else:
			basket.cluster_score = cluster_score / n
		
	def quant(self):
		"""Function analyzes the mzbin heatmap and cp heatmap that are loaded in the object
		and quantifies the effects of a basket on each parameter from cp. Removes baskets
		that have no effects on cp features and baskets that contain over 25% of the total
		idruns in the CP heatmap. Calculates the mean fingerprint for each basket based on
		the mean of the features for the idruns that contain adducts in the basket.
		"""
		print "Quantifying adduct/compound effects..."
		
		#Store list just in case to conserve order of contents and to avoid repeated generation of
		# the lists from the dictionaries
		features = self.cp.features()
		baskets = self.ad_heat.baskets()
		#A list of the feature names is needed later in the function

		widgets = ['HisDiff: ', Percentage(), ' ', Bar(marker=RotatingMarker()), ' ',\
				ETA(), ' ', FileTransferSpeed()]
		pbar = ProgressBar(widgets=widgets, maxval=len(self.ad_heat.baskets())).start()
		
		for ind, basket in enumerate(baskets):
			self.get_finger(basket, features)
			self.get_cluster(basket)
			pbar.update(ind + 1)
		pbar.finish()
		
	def reset_cscore(self, power):
	       """Function accepts a new power value to reset the cluster score
	       of each basket with.
	       """
	       for basket in self.ad_heat.baskets():
	           self.get_cluster(basket, power)


	def bask_prob(self, bask, idrun):
		"""Function takes a basket object from the mzbin module and an fingerprint object (idrun)
		from the cp module. Returns a likelihood score that the basket gave the fingerprint of that
		idrun based on the fingerprints (synthetic and measured).
		"""
		probs = [[score, bask.fingerprint[param][0], bask.fingerprint[param][1]] for param, score in idrun.items()]
		probs = zip(*[[mu-abs(mu-x), mu, sigma] for x, mu, sigma in probs if sigma != 0])
		if len(probs) == 0:
		    return 0
                score = 2 * norm.cdf(probs[0], loc=probs[1], scale=probs[2]).prod()
		return score

	def pearson_plot(self, name):
		"""Function generates a histogram of pearson correlations between all non
		overlapping combinations of synthetic fingerprints.
		"""
		basks = self.ad_heat.baskets()
		pearsons = []
		for bask in basks:
			for ads in itertools.combinations(bask.values(), 2):
				pearsons.append(self.cp.get_pearson(ads[0].idrun.name[:-1], ads[1].idrun.name[:-1]))
		plt.hist(pearsons, bins=1000, range=(-1.0, 1.0))
		plt.savefig(name+".pdf")
		plt.close()

	def cp_submz(self, idrun, pmax=None, maxt=0.4, mint=0.7, min_count=None, min_percent=None):
		"""Function returns a sub mz heatmap with the specified idrun, subtracting
		features present in anticorrelated runs (specified with pmax pearson score). Also
		uses maxt and mint to specify prefractions and compounds present in a single
		cluster with specified tolerances=. min_count and min_percent specify the minimum
		number of prefractions a feature must be in to not be filtered out.
		"""
		runs = self.cp.cluster(idrun.split("_")[0], max_tolerance=maxt, min_tolerance=mint)
		print "RUNS: ", runs
		if not pmax == None:
			antiruns = self.cp.anticluster(idrun.split("_")[0], pmax=pmax, fraction=1)
			print "ANTIRUNS: ", antiruns
		
		if idrun not in self.ad_heat.map:
			return
		submap = self.ad_heat.subheat([idrun])
		
		if not pmax == None:
			submap = submap.subheat_restrict(runs, antiruns)
		else:
			submap = submap.subheat_restrict(runs)
			
		run_num = float(len(submap.idruns()))	
		if min_count and min_percent:
			if min_count > min_percent * run_num:
				submap.remove_undercount(min_count)
			else:
				submap.remove_undercount(min_percent * run_num)
		elif min_percent:
			submap.remove_undercount(min_percent * run_num)
		elif min_count:
			submap.remove_undercount(min_count)	
		return submap


	def bask_network(self):
		"""Function exports a network of fingerprints as nodes and edges as
		compounds.
		"""
		h = [nx.Graph(), nx.Graph(), nx.Graph()]
		edges = dict()
		for id in self.cp.fingerprints():
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
		for pref in self.cp.fingerprints():
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
		submap = self.cp_submz(idrun, pmax, maxt, mint, min_count, min_percent)
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
		submap = self.cp_submz(idrun, None, maxt, mint)
		if submap == None:
			print "Empty Heatmap!"
			return
		basks = submap.baskets()
		plt.rc('xtick', labelsize=10)
		plt.rc('ytick', labelsize=10)
		fig = plt.figure(figsize=figsize, dpi=300)
		ax = fig.add_subplot(111)
		ax.set_ylim(0, 60)
		ax.set_xlim(-0.5, 1)
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

		
def parse_args(args):
	"""
	"""
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
			description=__doc__)
	parser.add_argument("--cpfile", dest="cp", action="store", help="""The CP file to be
			analyzed along with the mass spectral data.""")
	parser.add_argument("--mzdb", "-m", dest="mzdb", action="store", help="""The Mass Spectral
			database to be analyzed along with the CP data.""")
	parser.add_argument("--recalculate", "-r", dest="recalculate", action="store_true", default=
			False, help="Should CP data be recalculated?")
	parser.add_argument("--start", "-s", dest="start", action="store", default=0, type=int,
	        help="CP iteration start value")
	parser.add_argument("--bqs", "-b", dest="bqs", action="store", default="10100,10110,11110",
			help="BQ values to be loaded from the database. Comma separated.")
	parser.add_argument("--ascore", "-a", dest="ascore", action="store", default=8,
			type=int, help="Activity Score cutoff for removing inactive features. Default is 8.")
	parser.add_argument("--cscore", "-c", dest="cscore", action="store", default=0.1,
			type=int, help="Cluster Score cutoff for removing inactive features. Default is 0.1.")

	args = parser.parse_args(args)
	args.bqs = args.bqs.split(",")
	return args
		
def main(args):
	options = parse_args(args)
	g = cpbin(options.mzdb ,options.cp, bqs=['00100', '10100', '10110', '11110'])
	t = open("meta.txt", "w")
	print "Running with Activity Score: {} and Cluster Score: {}".format(options.ascore, options.cscore)
	g.recalculate_cp(ascore=options.ascore, cscore=options.cscore, output=t, start = options.start)
	t.close()

if __name__ == "__main__":
	main(sys.argv[1:])

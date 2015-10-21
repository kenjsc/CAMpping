#Emerson Glassey
# 07-30-2012

#------------------------------------------------------------------------------#
# This module is for the analysis of peaks/compounds in a database. An optional
# kmeans nxn csvfile is supplied along with relevant tight and loose cutoffs.
# 

import numpy
import dbtables
import sqlalchemy
from sqlalchemy import create_engine, func
from sqlalchemy.orm import sessionmaker, subqueryload_all, deferred
from progressbar import ProgressBar, Percentage, Bar, RotatingMarker, ETA, FileTransferSpeed
import cp
import mzbin as mb
import pp
import btools
import sys

import matplotlib
matplotlib.use("pdf")
#from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
#import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA

from scipy.optimize import curve_fit
from scipy.stats import pearsonr, norm
from scipy.misc import logsumexp

import itertools
import argparse

import networkx as nx



class cpbin():
	
	def __init__(self, db, tabfile, quant=True, remove=True, bqs=['00110', '00100', '10000', '10100', '10110',
			'11000', '11100', '11110','00111', '00101', '10001', '10101', '10111', '11001', '11101',
			'11111']):
		bqs = set(bqs)
		self.db = db
		self.engine = dbtables.connect(db)
		self.Session = sessionmaker(bind=self.engine)
		self.session = self.Session()
		
		self.cp = self.get_cp(tabfile)
		self.ad_heat = self.get_bins([finger.name for finger in self.cp.fingerprints()], bqs)
		
		if remove:
			self.remove_systematic()
			self.remove_underrepresent()
		
		if quant:
			self.quant()
		
		
	def get_cp(self, tabfile):
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
		corrected_runs = []
		for idrun in idruns:
			if len(idrun) == 5:
				corrected_runs.append(idrun + "_")
			else:
				corrected_runs.append(idrun)
		runs = self.get_idruns(corrected_runs, bqs)

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
		for id in self.cp.fingerprints():
			if id.activity_score() <= min_effect and id.name + "_" in self.ad_heat.map:
				self.ad_heat.remove_run_baskets(id.name + "_")
				
		
	def get_finger(self, basket, features):
		"""Function accepts a list of idrun_names, and a list of feature objects from the cp module.
		Using scores in the features, the function returns a fingerprint representative of the list
		of idruns by constructing the average fingerprint (mean of each feature between the idruns).
		"""
		fingerprint = []
		idrun_names = basket.keys()
		activity = 0
		for feat in features:
			in_scores = [feat[run[:5]] for run in idrun_names]
			mu = numpy.mean(in_scores)
			std = numpy.std(in_scores)
			basket.param_score(feat.name, (mu, std))
			activity += mu ** 2
		basket.activity = activity
		
	def get_cluster(self, basket, power=3):
		n = 0
		cluster_score = 0
		for pair in itertools.combinations(basket.keys(), 2):
			cluster_score += self.cp.get_pearson(pair[0].split("_")[0], pair[1].split("_")[0]) ** power
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
		feat_names = [feat.name for feat in features]

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

	def recalculate_cp(self, ascore = 8, cscore = 0.1, output = sys.stdout, start=0):
		self.remove_nonactive(ascore)
		self.remove_nonspecific(cscore)
#		self.remove_inactive_runs(5)
		cp_heat = self.cp
#		cp_runs = [str(finger) for finger in cp_heat.fingerprints()]
		cp_heat.write_tab("cp_{}.tab".format(start))
		bnets = self.bask_network()
		nx.write_dot(bnets[0], "cp_{}_bask_sum_network.dot".format(start))
		nx.write_dot(bnets[1], "cp_{}_bask_mean_network.dot".format(start))
		nx.write_dot(bnets[2], "cp_{}_bask_max_network.dot".format(start))
		nx.write_dot(cp_heat.cp_network(0.9), "cp_{}_cp_network.dot".format(start))
#		next = self
		for n in range(start, 1000):
			output.write("iteration {}\n".format(n+1))
			output.flush()
			
			cp_heat = self.all_chem_adjust(self.ad_heat, cp_heat, output)
			cp_heat.write_tab("cp_{}.tab".format(n+1))
			cp_heat.pearson()
			nx.write_dot(cp_heat.cp_network(0.9), "cp_{}_cp_network.dot".format(n+1))
			
#			cp_heat = self.chem_adjust(next.ad_heat, next.cp, output)
#			cp_heat.write_tab("cp_{}.tab".format(n+1))
#			next = cpbin(self.db, "cp_{}.tab".format(n+1), bqs=['10100', '10110', '11110'])
#			next.remove_nonactive(8)
#			next.remove_nonspecific(0.1)
#			next.remove_inactive_runs(5)
#			nx.write_dot(next.cp_network(), "cp_{}_cp_network.dot".format(n+1))
#			nx.write_dot(next.bask_network(), "cp_{}_bask_network.dot".format(n+1))
		
		return
		
	def chem_adjust(self, mz_heat, cp_heat, output):
		new_cp = cp.cp(None)
		for feat in cp_heat.features():
			new_cp[str(feat)] = cp.feature(str(feat))
		
		runcount = 0
		n = 0
		av_dist = 0
			
		widgets = ['VectorMove: ', Percentage(), ' ', Bar(marker=RotatingMarker()), ' ',\
		ETA(), ' ', FileTransferSpeed()]
		
		pbar = ProgressBar(widgets=widgets, maxval=len(cp_heat.fingerprints())).start()
		
		
		
		for run in cp_heat.fingerprints():
			output.write("\t" + str(run) + "\n")
			runcount += 1
			pbar.update(runcount)

			all_basks = [bask for bask in mz_heat.grab_basks(str(run))]

			inruns = set([inrun + "_" for inrun in cp_heat.cluster(str(run), max_tolerance=0.5, min_tolerance=0.65)])
#			antiruns = set([antirun + "_" for antirun in cp_heat.anticluster(str(run), pmax=-0.2, fraction=1)])
			
			basks = []
			for bask in all_basks:
				bruns = set(bask.keys())
#				if len(bruns & antiruns) > 0:
#					continue
				if len(bruns & inruns) < 2:
					continue
				basks.append(bask)
			
			labels = run.keys()
			run_vec = numpy.array(run.values(), dtype=float)
			add_vector = numpy.zeros(len(labels))
			vector_num = 0
			largest_scaler = 0.0
			for bask in basks:
				av_scaler = 0.0
				for connect_run in bask.keys():
					cprun = connect_run.replace("_", "")
					if cprun == str(run):
						continue
					if cprun not in cp_heat.map:
						continue
					vec_dif = numpy.array([cp_heat[cprun][val] for val in labels]) - run_vec
					scaler = self.bask_prob(bask, cp_heat[cprun]) * self.bask_prob(bask, run)
					if scaler >= largest_scaler:
						largest_scaler = scaler
					av_scaler += scaler
					add_vector += vec_dif * scaler
					vector_num += 1
				if not vector_num == 0:
					av_scaler /= vector_num
					output.write("\t\t\t{}; {}\n".format(str(bask), av_scaler))
			if not vector_num == 0:
				add_vector /= vector_num
			if not largest_scaler == 0:
				add_vector /= 2 * largest_scaler

				run_vec += add_vector
			new_cp[str(run)] = cp.fingerprint(str(run))
			for param, value in zip(labels, run_vec):
				new_cp[param][str(run)] = new_cp[str(run)][param] = value
				
			output.write("\t\tRun Movement: {}\n".format(numpy.sqrt(add_vector.dot(add_vector))))
				
			av_dist += numpy.sqrt(add_vector.dot(add_vector))
			n += 1
		pbar.finish()

		output.write("Average Movement: " + str(av_dist / n) + "\n")			
		return new_cp
		
	def all_chem_adjust(self, mz_heat, cp_heat, output):
		"""Does vector compression, ignoring biological activity. Uses all baskets from a run, not
		just those that are within a biological cutoff range.
		"""
		new_cp = cp.cp(None)
		for feat in cp_heat.features():
			new_cp[str(feat)] = cp.feature(str(feat))
		
		runcount = 0
		n = 0
		av_dist = 0
			
		widgets = ['VectorMove: ', Percentage(), ' ', Bar(marker=RotatingMarker()), ' ',\
		ETA(), ' ', FileTransferSpeed()]
		
		pbar = ProgressBar(widgets=widgets, maxval=len(cp_heat.fingerprints())).start()
		
		largest_scaler = None
		#This is a dictionary of run (as string) - [run_vector, add_vector] pairs
		add_vectors = dict()
		for run in cp_heat.fingerprints():
			output.write("\t" + str(run) + "\n")
			runcount += 1
			pbar.update(runcount)
			
			labels = run.keys()
			#This is the original vector fingerprint in log scale
			run_vec = numpy.log(run.values(), dtype=float)
			print run_vec
			#This is going to be a list of vectors, one from each basket, in log scale
			bask_vectors = []
			#This is the number of vectors, one for each connection, with multiple connections per basket
			vector_num = 0
			for bask in mz_heat.grab_basks(str(run)):
			        #If there's only one run, then the vector has nothing to connect to
				if len(bask.keys()) <= 1:
					continue
				#This is the average value of the scaler, for use in line plots later
				av_scaler = 0.0
				for connect_run in bask.keys():
					cprun = connect_run.replace("_", "")
					#Don't connect the query run to itself, that's not useful
					if cprun == str(run):
						continue
					#If the run isn't in the cp_heatmap, just continue; that means it also wasn't used for creating synthetic fingerprints
					if cprun not in cp_heat.map:
						continue
					#Get the vector difference between the target and the source, but it's in log scale. Also, make sure label values are
					# in the same order between the vectors. Remember log scale!
					vec_dif = numpy.log(numpy.array([cp_heat[cprun][val] for val in labels]) - run_vec)
					scaler = self.bask_prob(bask, cp_heat[cprun]) + self.bask_prob(bask, run)
					print scaler
					if scaler >= largest_scaler:
						largest_scaler = scaler
					av_scaler = numpy.logaddexp(av_scaler, scaler)
					bask_vectors.append(vec_dif + scaler)
					print "bask_ind", vec_dif+scaler
					vector_num += 1
				if not vector_num == 0:
					output.write("\t\t\t{}; {}\n".format(str(bask), numpy.exp(av_scaler) / vector_num))
			#This is the sum of the basket vectors, still in log scale
			print "all", bask_vectors
			add_vector = logsumexp(bask_vectors, axis=0)
			print "summed", add_vector
			if not vector_num == 0:
				add_vector -= numpy.log(vector_num)
			print "averaged", add_vector
			print 'large', largest_scaler
			return
#			if not largest_scaler == 0:
#				add_vector /= 2 * largest_scaler
                        add_vectors[str(run)] = [run_vec, add_vector]
                
                for run, (run_vec, add_vector) in add_vectors.items():
			add_vector = numpy.exp(add_vector - (numpy.log(2) + largest_scaler))
			run_vec = numpy.exp(run_vec) + add_vector
			new_cp[str(run)] = cp.fingerprint(str(run))
			for param, value in zip(labels, run_vec):
				new_cp[param][str(run)] = new_cp[str(run)][param] = value
				
			output.write("\t\tRun Movement: {}\n".format(numpy.sqrt(add_vector.dot(add_vector))))
				
			av_dist += numpy.sqrt(add_vector.dot(add_vector))
			n += 1
		pbar.finish()

		output.write("Average Movement: " + str(av_dist / n) + "\n")			
		return new_cp
		
	def bask_prob(self, bask, idrun):
		"""Function takes a basket object from the mzbin module and an fingerprint object (idrun)
		from the cp module. Returns a likelihood score that the basket gave the fingerprint of that
		idrun based on the fingerprints (synthetic and measured).
		"""
#		likelihood = 1
#		for param, score in idrun.items():
#			mu, sigma = bask.fingerprint[param]
#			if sigma == 0:
#			     continue
#			likelihood *= self.prob_score(score, mu, sigma)
		probs = [[score, bask.fingerprint[param][0], bask.fingerprint[param][1]] for param, score in idrun.items()]
		probs = zip(*[[mu-abs(mu-x), mu, sigma] for x, mu, sigma in probs if sigma != 0])
		if len(probs) == 0:
		    return 0
		return numpy.log(2) + numpy.sum(numpy.log(norm.cdf(probs[0], loc=probs[1], scale=probs[2])))
#		return likelihood
		
	def prob_score(self, x, mu, sigma):
		"""Function takes a value of x and an average and standard deviation from a normal
		distribution. Returns the probability of the value x given the distribution.
		"""
		if sigma == 0:
			return 1.0
		return 2 * norm.cdf(mu - abs(mu - x), loc=mu, scale=sigma)

	def pearson_plot(self, name):
		basks = self.ad_heat.baskets()
		pearsons = []
		for bask in basks:
			for ads in itertools.combinations(bask.values(), 2):
				pearsons.append(self.cp.get_pearson(ads[0].idrun.name[:-1], ads[1].idrun.name[:-1]))
		plt.hist(pearsons, bins=1000, range=(-1.0, 1.0))
		plt.savefig(name+".pdf")
		plt.close()
		
	def isotope_plot(self, name):
		isos = []
		for b in self.ad_heat.baskets():
			for comb in itertools.combinations(b.values(), 2):
				isos.append(btools.iso_diff(comb[0].get_pat(), comb[1].get_pat()))
		plt.hist(isos, bins=1000, range=(0.0, 1.0))
		plt.savefig(name + ".pdf")

	def cp_submz(self, idrun, pmax=None, maxt=0.4, mint=0.7, min_count=None, min_percent=None):
		idrun = idrun + "_" if not idrun.endswith("_") else idrun
		runs = [run + "_" for run in self.cp.cluster(idrun.split("_")[0], max_tolerance=maxt,\
				min_tolerance=mint)]
		print "RUNS: ", runs
		if not pmax == None:
			antiruns = [run + "_" for run in self.cp.anticluster(idrun.split("_")[0], pmax=pmax,\
					fraction=1)]
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
#		features = self.cp.features()
#		for bask in submap.baskets():
#			self.get_finger(bask, features)
#			self.get_cluster(bask)
		for bask in submap.baskets():
			bask.activity = self.ad_heat[bask.name].activity
			bask.cluster_score = self.ad_heat[bask.name].cluster_score
			bask.fingerprint = self.ad_heat[bask.name].fingerprint
		
		return submap


	def bask_network(self):
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
	


		
#	def weighter(self, bask, id1, id2):
#		return self.bask_prob(bask, self.cp[id1]) * self.bask_prob(bask, self.cp[id2])
#		
#	def network(self, weight=weighter):
#		h = nx.MultiGraph()
#		for id in self.cp.idruns():                                 
#			h.add_node(id.name + "_", weight=id.activity_score())
#		
#		weights = dict()
#		nums = dict()
#		
#		for pair in [pair for bask in self.ad_heat.baskets() for pair in itertools.combinations(sorted(bask.keys()), 2)]:
#			h.add_edge(*pair, weight = 0.0, num_share = 0.0)
#				weights[pair] = 0
#				nums[pair] = 0
#			
#			weights[pair] *= weights[pair]
#			nums[pair] += 1
#			weights[pair] += weight(self, bask, pair[0][:-1], pair[1][:-1])
#			weights[pair] /= nums[pair]
#		nx.set_edge_attributes(h, 'weight', weights)
#		nx.set_edge_attributes(h, 'num_share', nums)
#
#		return h
#		
#	def mz_network(self):
#		h = nx.Graph()
#		for feat in self.ad_heat.baskets():
#			h.add_node(feat.name, m=feat.m, rt=feat.rt, ions=",".join(feat.ions), activity=feat.activity, cluster=feat.cluster_score,num=len(feat.peaks))		
#		
#		nums = dict()
#		for pair in [pair for id in self.ad_heat.idruns() for pair in itertools.combinations(sorted(id.keys()), 2)]:
#			if pair not in h.edges():
#				h.add_edge(*pair, num=0.0)
#				nums[pair] = 0
#			
#			nums[pair] += 1
#		nx.set_edge_attributes(h, 'num_share', nums)
#		return h
	
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
		name = idrun if not name else name
		self.plot(name, nxn, cols, rows, title, text)
		
	def export_polar(self, idrun, maxt=0.5, mint=0.65, text=1, scale=1, greyscale='0.5', ylim=1200):
                plt.close()
		submap = self.cp_submz(idrun, None, maxt, mint)
		basks = submap.baskets()
#		total = set([p for b in basks for p in b.peaks])
		if len(basks) == 0 or len(submap.idruns()) == 0:
			print "empty heatmap"
			return
#		plt.rc('grid', linewidth=0.5, color=greyscale, linestyle='-')
#		plt.rc('xtick', labelsize=10)
#		plt.rc('ytick', labelsize=10)
		fig = plt.figure()
		ax = fig.add_subplot(111)
#		ax.set_ylim(0, ylim)
#		ticklen = ax.set_yticks(numpy.arange(0, ylim, 200))
#		ax.set_yticklabels([''] * len(ticklen))
#		ticklen = ax.set_xticks(numpy.linspace(0.0, 2 * numpy.pi, 13))
#		ax.set_xticklabels([''] * len(ticklen))
#		maxi = max([b.activity*b.cluster_score for b in basks])
		for bask in sorted(basks, key=lambda b: b.activity * b.cluster_score, reverse=True):
			ax.plot([bask.rt * numpy.pi / 3], [bask.m], 'bo', ms = bask.activity * bask.cluster_score * scale)
		plt.savefig(idrun + "_polar.pdf", bbox_inches='tight')
		plt.close()
		
		
	def plot(self, name, matrix, bins, rows, title="", text=4):

		cm = {'red':((0.0, 1.0, 1.0), (0.0000000000000000000000000001, 1.0, 0.0), (1.0, 1.0, 1.0)),\
          		'green':((0.0, 1.0, 1.0), (0.0000000000000000000000000001, 0.0, 0.0), (1.0, 0.0, 0.0)),\
			'blue':((0.0, 1.0, 1.0), (0.00000000000000000000000000001, 1.0, 1.0), (1.0, 0.0, 0.0))}
		mycm = matplotlib.colors.LinearSegmentedColormap('mycm', cm)
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
		
	def export_activity(self, idrun, maxt=0.5, mint=0.65, size=10, aline=8, cline=0.1, alpha=0.75, figsize=(2, 2)):
	       cm = {'red':((0.0, 1.0, 1.0), (0.0000000000000000000000000001, 1.0, 0.0), (1.0, 1.0, 1.0)),\
          		'green':((0.0, 1.0, 1.0), (0.0000000000000000000000000001, 0.0, 0.0), (1.0, 0.0, 0.0)),\
			'blue':((0.0, 1.0, 1.0), (0.00000000000000000000000000001, 1.0, 1.0), (1.0, 0.0, 0.0))}
	       mycm = matplotlib.colors.LinearSegmentedColormap('mycm', cm)
	       submap = self.cp_submz(idrun, None, maxt, mint)
	       basks = submap.baskets()
	       if len(basks) == 0 or len(submap.idruns()) == 0:
	               print 'empty heatmap'
	               return
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
	
	       
        def plot_likelihood(self, idrun, maxt=0.5, mint=0.65, size=10, aline=8, lline=0.15, alpha=0.9, figsize=(5,5)):
	       cm = {'red':((0.0, 1.0, 1.0), (0.0000000000000000000000000001, 1.0, 0.0), (1.0, 1.0, 1.0)),\
          		'green':((0.0, 1.0, 1.0), (0.0000000000000000000000000001, 0.0, 0.0), (1.0, 0.0, 0.0)),\
			'blue':((0.0, 1.0, 1.0), (0.00000000000000000000000000001, 1.0, 1.0), (1.0, 0.0, 0.0))}
	       mycm = matplotlib.colors.LinearSegmentedColormap('mycm', cm)
	       submap = self.cp_submz(idrun, None, maxt, mint)
	       basks = submap.baskets()
	       if len(basks) == 0 or len(submap.idruns()) == 0:
	               print 'empty heatmap'
	               return
	       plt.rc('xtick', labelsize=10)
	       plt.rc('ytick', labelsize=10)
	       fig = plt.figure(figsize=figsize, dpi=300)
	       ax = fig.add_subplot(111)
	       ax.set_ylim(0, 60)
	       ax.set_xscale('log')
#	       ax.set_xlim(0, 1)
	       ax.axhline(aline, color='g')
	       ax.axvline(lline, color='g')
	       ax.set_xlabel("Likelihood Score", size=10)
	       ax.set_ylabel("Activity*Cluster Score", size=10)
	       
	       x = [self.bask_prob(bask, self.cp[idrun]) for bask in basks]
	       y = [bask.activity*bask.cluster_score for bask in basks]
	       c = [bask.rt for bask in basks]
	       
	       plt.scatter(x, y, s=size, c=c, cmap=mycm, vmin=0, vmax=6, alpha=alpha)
	       plt.savefig(idrun + "_likelihood_plot.pdf", bbox_inches='tight')
	       plt.close()
	       return
	       
	       
		
def parse_args(args):
	"""
	"""
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,\
			description=__doc__)
	parser.add_argument("--cpfile", "-c", dest="cp", action="store", help="""The CP file to be
			analyzed along with the mass spectral data.""")
	parser.add_argument("--mzdb", "-m", dest="mzdb", action="store", help="""The Mass Spectral
			database to be analyzed along with the CP data.""")
	parser.add_argument("--recalculate", "-r", dest="recalculate", action="store_true", default=\
			False, help="Should CP data be recalculated?")
	parser.add_argument("--start", "-s", dest="start", action="store", default=0, type=int, \
	                help="CP iteration start value")

	args = parser.parse_args(args)
	return args
		
def main(args):
	options = parse_args(args)
	g = cpbin(options.mzdb ,options.cp, bqs=['10100', '10110', '11110'])
	t = open("meta.txt", "w")
	ascore = 8
	cscore = 0.15
	print "Running with Activity Score: {} and Cluster Score: {}".format(ascore, cscore)
	g.recalculate_cp(ascore=ascore, cscore=cscore, output=t, start = options.start)
	t.close()




if __name__ == "__main__":
	main(sys.argv[1:])


"""Module contains classes defining objects for parsing and indexing cytological profiling data.
"""

#import scipy.stats
from scipy.stats import pearsonr
import numpy
from progressbar import ProgressBar, Percentage, Bar, RotatingMarker, ETA, FileTransferSpeed
import sys

import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
import networkx as nx
import itertools
import pp

import re


labels = "MICRO_AND_BINUCLEATED_CELLS_MICRONUCLEI_EDU	MICRO_AND_PROBEABPOSITIVE_CELLS_MICRONUCLEI_EDU	MICRO_AND_PROBEBPOSITIVE_CELLS_MICRONUCLEI_EDU	MICRO_AND_PROBEAPOSITIVE_CELLS_MICRONUCLEI_EDU	INTERPHASE_CELLS_MICRONUCLEI_EDU	MONONUCLEATED_CELLS_MICRONUCLEI_EDU	BINUCLEATED_CELLS_MICRONUCLEI_EDU	MICRONUCLEATED_CELLS_MICRONUCLEI_EDU	CELLS_WITH_ONE_MICRONUCLEUS_MICRONUCLEI_EDU	PCT_BINUCLEATED_CELLS_MICRONUCLEI_EDU	NUCLEAR_DIVISION_INDEX_MICRONUCLEI_EDU	PCT_CELLS_WITH_ONE_MICRONUCLEUS_MICRONUCLEI_CYTO	PCT_MICRONUCLEATED_CELLS_MICRONUCLEI_CYTO	MICRONUCLEI_PER_CELL_MICRONUCLEI_CYTO	TOTAL_MICRONUCLEI_MICRONUCLEI_CYTO	ALL_ACTIN_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_ACTIN_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_CYTO	ALL_ACTIN_MEAN_STAIN_INTEGR_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_ACTIN_MEAN_STAIN_INTEGR_INTENS_MULTIWAVESCORING_CYTO	ALL_ACTIN_MEAN_STAIN_AREA_MULTIWAVESCORING_CYTO	POSITIVE_ACTIN_MEAN_STAIN_AREA_MULTIWAVESCORING_CYTO	POSITIVE_TUBULIN_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_CYTO	ALL_TUBULIN_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_CYTO	CELL_TUBULIN_NUCLEUS_INTEGR_INTENSITY_MULTIWAVESCORING_CYTO	DNA_MEAN_AREA_MICRONUCLEI_EDU	BREADTH_IMA_SUMMARY_CYTO	INNER_RADIUS_IMA_SUMMARY_CYTO	OUTER_RADIUS_IMA_SUMMARY_CYTO	LENGTH_IMA_SUMMARY_CYTO	PERIMETER_IMA_SUMMARY_CYTO	MEAN_RADIUS_IMA_SUMMARY_CYTO	EQUIV_RADIUS_IMA_SUMMARY_CYTO	EQUIV_PROLATE_VOL_IMA_SUMMARY_CYTO	ALL_NUCLEI_MEAN_AREA_MULTIWAVESCORING_CYTO	AREA_IMA_SUMMARY_CYTO	EQUIV_SPHERE_SURFACE_AREA_IMA_SUMMARY_CYTO	TOTAL_AREA_IMA_SUMMARY_CYTO	PIXEL_AREA_IMA_SUMMARY_CYTO	EQUIV_SPHERE_VOL_IMA_SUMMARY_CYTO	EQUIV_OBLATE_VOL_IMA_SUMMARY_CYTO	NEGATIVE_EDU_MEAN_CELL_INTEGR_INTENS_MULTIWAVESCORING_EDU	NEGATIVE_EDU_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_EDU	OUTER_RADIUS_IMA_SUMMARY_EDU	LENGTH_IMA_SUMMARY_EDU	EQUIV_SPHERE_VOL_IMA_SUMMARY_EDU	EQUIV_OBLATE_VOL_IMA_SUMMARY_EDU	EQUIV_PROLATE_VOL_IMA_SUMMARY_EDU	DAPI_STAINED_AREA_MULTIWAVESCORING_EDU	CELL_TOTAL_AREA_MULTIWAVESCORING_EDU	WIDTH_IMA_SUMMARY_EDU	PERIMETER_IMA_SUMMARY_EDU	EQUIV_RADIUS_IMA_SUMMARY_EDU	MEAN_RADIUS_IMA_SUMMARY_EDU	ALL_CELLS_MEAN_AREA_MULTIWAVESCORING_EDU	ALL_NUCLEI_MEAN_AREA_MULTIWAVESCORING_EDU	PIXEL_AREA_IMA_SUMMARY_EDU	EQUIV_SPHERE_SURFACE_AREA_IMA_SUMMARY_EDU	AREA_IMA_SUMMARY_EDU	TOTAL_AREA_IMA_SUMMARY_EDU	NUCLEAR_AREA_PER_CELL_TRANSFLUOR_EDU	INNER_RADIUS_IMA_SUMMARY_EDU	BREADTH_IMA_SUMMARY_EDU	ALL_PH3_MEAN_CELL_AVER_INTENS_MULTIWAVESCORING_EDU	ALL_PH3_MEAN_NUCL_AVER_INTENS_MULTIWAVESCORING_EDU	ALL_PH3_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_EDU	ALL_PH3_MEAN_CELL_INTEGR_INTENS_MULTIWAVESCORING_EDU	NEGATIVE_PH3_MEAN_CELL_INTEGR_INTENS_MULTIWAVESCORING_EDU	NEGATIVE_PH3_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_EDU	PCT_CELLS_WITH_MULTI_MICRONUCLEI_MICRONUCLEI_EDU	ALL_EDU_MEAN_STAIN_AREA_MULTIWAVESCORING_EDU	PCT_POSITIVE_EDU_MULTIWAVESCORING_EDU	CELL_PH3_STAINED_AVERAGE_INTENSITY_MULTIWAVESCORING_CYTO	CELL_PH3_STAINED_AVERAGE_INTENSITY_MULTIWAVESCORING_EDU	NEGATIVE_PH3_MEAN_STAIN_AREA_MULTIWAVESCORING_EDU	NEGATIVE_PH3_MEAN_STAIN_INTEGR_INTENS_MULTIWAVESCORING_EDU	POSITIVE_ACTIN_MEAN_CELL_AVER_INTENS_MULTIWAVESCORING_CYTO	ALL_ACTIN_MEAN_CELL_AVER_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_ACTIN_MEAN_CYTO_AVER_INTENS_MULTIWAVESCORING_CYTO	ALL_ACTIN_MEAN_CYTO_AVER_INTENS_MULTIWAVESCORING_CYTO	PIT_INTEGRATED_INTENSITY_TRANSFLUOR_CYTO	PIT_COUNT_TRANSFLUOR_CYTO	PIT_TOTAL_AREA_TRANSFLUOR_CYTO	GRADIENT_INDEX_TRANSFLUOR_CYTO	LAPLACIAN_INDEX_TRANSFLUOR_CYTO	POSITIVE_TUBULIN_MEAN_CYTO_AVER_INTENS_MULTIWAVESCORING_CYTO	ALL_TUBULIN_MEAN_CYTO_AVER_INTENS_MULTIWAVESCORING_CYTO	ALL_TUBULIN_MEAN_STAIN_AVER_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_TUBULIN_MEAN_STAIN_AVER_INTENS_MULTIWAVESCORING_CYTO	ALL_TUBULIN_MEAN_CELL_AVER_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_TUBULIN_MEAN_CELL_AVER_INTENS_MULTIWAVESCORING_CYTO	NEGATIVE_ACTIN_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_CYTO	NEGATIVE_ACTIN_MEAN_STAIN_AVER_INTENS_MULTIWAVESCORING_CYTO	NEGATIVE_ACTIN_MEAN_CYTO_AVER_INTENS_MULTIWAVESCORING_CYTO	NEGATIVE_ACTIN_MEAN_CELL_AVER_INTENS_MULTIWAVESCORING_CYTO	NEGATIVE_ACTIN_MEAN_NUCL_AVER_INTENS_MULTIWAVESCORING_CYTO	MITOTIC_CELLS_MICRONUCLEI_CYTO	TOTAL_CELLS_MICRONUCLEI_CYTO	TOTAL_CELLS_MULTIWAVESCORING_CYTO	POSITIVE_TUBULIN_MULTIWAVESCORING_CYTO	POSITIVE_ACTIN_MULTIWAVESCORING_CYTO	DNA_TOTAL_AREA_MICRONUCLEI_CYTO	INTERPHASE_CELLS_MICRONUCLEI_CYTO	MONONUCLEATED_CELLS_MICRONUCLEI_CYTO	MITOTIC_CELLS_MICRONUCLEI_EDU	POSITIVE_EDU_MULTIWAVESCORING_EDU	PIT_TOTAL_AREA_TRANSFLUOR_EDU	PIT_COUNT_TRANSFLUOR_EDU	NUCLEAR_INTEGRATED_INTENSITY_TRANSFLUOR_EDU	NUCLEAR_TOTAL_AREA_TRANSFLUOR_EDU	DNA_TOTAL_AREA_MICRONUCLEI_EDU	TOTAL_CELLS_MICRONUCLEI_EDU	TOTAL_CELLS_MULTIWAVESCORING_EDU	NUCLEAR_COUNT_TRANSFLUOR_EDU	ELL_FORM_FACTOR_IMA_SUMMARY_CYTO	ELL_FORM_FACTOR_IMA_SUMMARY_EDU	NEGATIVE_PH3_MEAN_STAIN_AVER_INTENS_MULTIWAVESCORING_EDU	ALL_PH3_MEAN_STAIN_INTEGR_INTENS_MULTIWAVESCORING_EDU	ALL_PH3_MEAN_STAIN_AREA_MULTIWAVESCORING_EDU	PCT_POSITIVE_PH3_MULTIWAVESCORING_EDU	POSITIVE_PH3_MULTIWAVESCORING_EDU	POSITIVE_PH3_MULTIWAVESCORING_CYTO	PCT_POSITIVE_PH3_MULTIWAVESCORING_CYTO	ALL_PH3_MEAN_STAIN_AREA_MULTIWAVESCORING_CYTO	ALL_PH3_MEAN_STAIN_INTEGR_INTENS_MULTIWAVESCORING_CYTO	NEGATIVE_PH3_MEAN_STAIN_AREA_MULTIWAVESCORING_CYTO	NEGATIVE_PH3_MEAN_STAIN_INTEGR_INTENS_MULTIWAVESCORING_CYTO	PIT_AVERAGE_INTENSITY_TRANSFLUOR_EDU	NEGATIVE_EDU_MEAN_STAIN_AVER_INTENS_MULTIWAVESCORING_EDU	CELL_PROBEA_AVERAGE_INTENSITY_MICRONUCLEI_EDU	CELL_EDU_NUCLEUS_AVERAGE_INTENSITY_MULTIWAVESCORING_EDU	CELL_EDU_CELL_AVERAGE_INTENSITY_MULTIWAVESCORING_EDU	ALL_EDU_MEAN_CELL_INTEGR_INTENS_MULTIWAVESCORING_EDU	ALL_EDU_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_EDU	ALL_EDU_MEAN_STAIN_INTEGR_INTENS_MULTIWAVESCORING_EDU	POSITIVE_EDU_MEAN_CELL_AVER_INTENS_MULTIWAVESCORING_EDU	POSITIVE_EDU_MEAN_NUCL_AVER_INTENS_MULTIWAVESCORING_EDU	ALL_EDU_MEAN_STAIN_AVER_INTENS_MULTIWAVESCORING_EDU	POSITIVE_EDU_MEAN_STAIN_AVER_INTENS_MULTIWAVESCORING_EDU	POSITIVE_EDU_MEAN_STAIN_INTEGR_INTENS_MULTIWAVESCORING_EDU	POSITIVE_EDU_MEAN_CELL_INTEGR_INTENS_MULTIWAVESCORING_EDU	POSITIVE_EDU_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_EDU	CELLULAR_GRADIENT_INDEX_TRANSFLUOR_EDU	CELLULAR_LAPLACIAN_INDEX_TRANSFLUOR_EDU	PIT_INTEGRATED_INTENSITY_TRANSFLUOR_EDU	ALL_EDU_MEAN_NUCL_AVER_INTENS_MULTIWAVESCORING_EDU	ALL_EDU_MEAN_CELL_AVER_INTENS_MULTIWAVESCORING_EDU	CELLULAR_TEXTURE_INDEX_TRANSFLUOR_EDU	LAPLACIAN_INDEX_TRANSFLUOR_EDU	GRADIENT_INDEX_TRANSFLUOR_EDU	TEXTURE_INDEX_TRANSFLUOR_EDU	CELL_EDU_STAINED_AVERAGE_INTENSITY_MULTIWAVESCORING_EDU	POSITIVE_EDU_MEAN_STAIN_AREA_MULTIWAVESCORING_EDU	CELL_EDU_STAINED_AREA_MULTIWAVESCORING_EDU	POSITIVE_PH3_MEAN_STAIN_AREA_MULTIWAVESCORING_EDU	NEGATIVE_PH3_MEAN_STAIN_AVER_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_PH3_MEAN_CYTO_AVER_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_PH3_MEAN_CELL_AVER_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_PH3_MEAN_STAIN_AREA_MULTIWAVESCORING_CYTO	ALL_PH3_MEAN_STAIN_AVER_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_PH3_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_PH3_MEAN_STAIN_INTEGR_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_PH3_MEAN_STAIN_AVER_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_PH3_MEAN_NUCL_AVER_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_PH3_MEAN_CELL_INTEGR_INTENS_MULTIWAVESCORING_EDU	POSITIVE_PH3_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_EDU	POSITIVE_PH3_MEAN_STAIN_INTEGR_INTENS_MULTIWAVESCORING_EDU	POSITIVE_PH3_MEAN_STAIN_AVER_INTENS_MULTIWAVESCORING_EDU	POSITIVE_PH3_MEAN_NUCL_AVER_INTENS_MULTIWAVESCORING_EDU	POSITIVE_PH3_MEAN_CELL_AVER_INTENS_MULTIWAVESCORING_EDU	ALL_PH3_MEAN_STAIN_AVER_INTENS_MULTIWAVESCORING_EDU	NEGATIVE_TUBULIN_MEAN_NUCL_AVER_INTENS_MULTIWAVESCORING_CYTO	NEGATIVE_TUBULIN_MEAN_CELL_AVER_INTENS_MULTIWAVESCORING_CYTO	NEGATIVE_TUBULIN_MEAN_CYTO_AVER_INTENS_MULTIWAVESCORING_CYTO	NEGATIVE_TUBULIN_MEAN_STAIN_INTEGR_INTENS_MULTIWAVESCORING_CYTO	NEGATIVE_TUBULIN_MEAN_STAIN_AREA_MULTIWAVESCORING_CYTO	NEGATIVE_TUBULIN_MEAN_STAIN_AVER_INTENS_MULTIWAVESCORING_CYTO	NEGATIVE_TUBULIN_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_CYTO	NEGATIVE_TUBULIN_MEAN_CELL_INTEGR_INTENS_MULTIWAVESCORING_CYTO	PCT_POSITIVE_TUBULIN_MULTIWAVESCORING_CYTO	BINUCLEATED_CELLS_MICRONUCLEI_CYTO	MULTINUCLEATED_CELLS_MICRONUCLEI_CYTO	MULTIMONO_CELL_RATIO_MICRONUCLEI_CYTO	MULTIDUALMONO_CELL_RATIO_MICRONUCLEI_CYTO	DUALMONO_CELL_RATIO_MICRONUCLEI_CYTO	PCT_BINUCLEATED_CELLS_MICRONUCLEI_CYTO	MULTIDUALMONO_CELL_RATIO_MICRONUCLEI_EDU	DUALMONO_CELL_RATIO_MICRONUCLEI_EDU	DNA_MEAN_AREA_MICRONUCLEI_CYTO	ALL_PH3_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_CYTO	CELL_ACTIN_CYTOPLASM_AVERAGE_INTENSITY_MULTIWAVESCORING_CYTO	CELL_ACTIN_CELL_AVERAGE_INTENSITY_MULTIWAVESCORING_CYTO	NEGATIVE_ACTIN_MEAN_STAIN_INTEGR_INTENS_MULTIWAVESCORING_CYTO	NEGATIVE_ACTIN_MEAN_STAIN_AREA_MULTIWAVESCORING_CYTO	PCT_CELLS_WITH_MULTI_MICRONUCLEI_MICRONUCLEI_CYTO	SHAPE_FACTOR_IMA_SUMMARY_CYTO	PCT_PROBEA_POSITIVE_CELLS_MICRONUCLEI_EDU	PCT_PROBEB_POSITIVE_CELLS_MICRONUCLEI_EDU	PCT_PROBEAB_POSITIVE_CELLS_MICRONUCLEI_EDU	PCT_CELLS_WITH_ONE_MICRONUCLEUS_MICRONUCLEI_EDU	PCT_MICRONUCLEATED_CELLS_MICRONUCLEI_EDU	MICRONUCLEI_PER_CELL_MICRONUCLEI_EDU	MICRONUCLEI_PER_HEALTHY_CELL_MICRONUCLEI_EDU	MICRONUCLEI_PER_MONONUCLEATED_CELL_MICRONUCLEI_CYTO	MICRONUCLEI_PER_PROBEBPOSITIVE_CELL_MICRONUCLEI_EDU	MICRONUCLEI_PER_PROBEABPOSITIVE_CELL_MICRONUCLEI_EDU	MICRONUCLEI_PER_PROBEAPOSITIVE_CELL_MICRONUCLEI_EDU	MICRONUCLEI_PER_MONONUCLEATED_CELL_MICRONUCLEI_EDU	NEGATIVE_EDU_MEAN_STAIN_AREA_MULTIWAVESCORING_EDU	NEGATIVE_EDU_MEAN_STAIN_INTEGR_INTENS_MULTIWAVESCORING_EDU	CELLS_WITH_MULTI_MICRONUCLEI_MICRONUCLEI_CYTO	CELL_ACTIN_NUCLEUS_AVERAGE_INTENSITY_MULTIWAVESCORING_CYTO	CELL_ACTIN_STAINED_AVERAGE_INTENSITY_MULTIWAVESCORING_CYTO	TEXTURE_INDEX_TRANSFLUOR_CYTO	PIT_AVERAGE_INTENSITY_TRANSFLUOR_CYTO	CELL_TUBULIN_CELL_AVERAGE_INTENSITY_MULTIWAVESCORING_CYTO	ALL_TUBULIN_MEAN_NUCL_AVER_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_TUBULIN_MEAN_NUCL_AVER_INTENS_MULTIWAVESCORING_CYTO	CELLULAR_GRADIENT_INDEX_TRANSFLUOR_CYTO	CELL_GRADIENT_INDEX_TRANSFLUOR_CYTO	CELL_LAPLACIAN_INDEX_TRANSFLUOR_CYTO	CELL_TEXTURE_INDEX_TRANSFLUOR_CYTO	CELL_TUBULIN_NUCLEUS_AVERAGE_INTENSITY_MULTIWAVESCORING_CYTO	CELL_PIT_AVERAGE_INTENSITY_TRANSFLUOR_CYTO	PIT_COUNT_PER_CELL_TRANSFLUOR_CYTO	PIT_AREA_PER_CELL_TRANSFLUOR_CYTO	POSITIVE_ACTIN_MEAN_NUCL_AVER_INTENS_MULTIWAVESCORING_CYTO	MULTINUCLEATED_CELLS_MICRONUCLEI_EDU	MULTIMONO_CELL_RATIO_MICRONUCLEI_EDU	PCT_MULTINUCLEATED_CELLS_MICRONUCLEI_EDU	CELLS_WITH_MULTI_MICRONUCLEI_MICRONUCLEI_EDU	ALL_NUCLEI_MEAN_INTEGR_ITENSITY_MULTIWAVESCORING_EDU	TOTAL_INTENSITY_IMA_SUMMARY_EDU	ALL_NUCLEI_MEAN_INTEGR_ITENSITY_MULTIWAVESCORING_CYTO	TOTAL_INTENSITY_IMA_SUMMARY_CYTO	PCT_MITOTIC_CELLS_MICRONUCLEI_EDU	ALL_NUCLEI_MEAN_AVERAGE_ITENSITY_MULTIWAVESCORING_CYTO	AVERAGE_INTENSITY_IMA_SUMMARY_CYTO	CELL_DNA_AVERAGE_INTENSITY_MICRONUCLEI_CYTO	CELL_NUCLEAR_AVERAGE_INTENSITY_MICRONUCLEI_CYTO	DAPI_AVERAGE_INTENSITY_MULTIWAVESCORING_CYTO	CELL_NUCLEAR_AVERAGE_INTENSITY_TRANSFLUOR_CYTO	CELL_DNA_AVERAGE_INTENSITY_MICRONUCLEI_EDU	CELL_NUCLEAR_AVERAGE_INTENSITY_MICRONUCLEI_EDU	CELL_NUCLEAR_AVERAGE_INTENSITY_TRANSFLUOR_EDU	DAPI_AVERAGE_INTENSITY_MULTIWAVESCORING_EDU	NUCLEAR_AVERAGE_INTENSITY_TRANSFLUOR_EDU	ALL_NUCLEI_MEAN_AVERAGE_ITENSITY_MULTIWAVESCORING_EDU	AVERAGE_INTENSITY_IMA_SUMMARY_EDU".split("\t")

#def pearsonr(v1, v2):
#	return scipy.stats.pearsonr(v1, v2)[0]

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
		self.parent = None
		self.children = None
		self.pearson = dict()

	def __getitem__(self, feat):
		return self.finger[feat]
	
	def __setitem__(self, feat, value):
		self.finger[feat] = value
		
	def __delitem__(self, feat):
		del self.finger[feat]
		
	def __repr__(self):
		return self.name
		
	def setp(self, id, p):
		self.pearson[id] = p
		
	def getp(self, id):
		return self.pearson[id]
		
	def delp(self, id):
		del self.pearson[id]
		
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
	def __init__(self, cpfile, pearson=True):
		self.map = dict()
		if cpfile:
			self.name = cpfile.split("/")[-1].split(".")[0]
			feats, runs = self.parse_tab(open(cpfile, "U"))
			for run in runs:
				self.map[run.name] = run
			for feat in feats:
				self.map[feat.name] = feat
			if pearson:
			        self.pearson()

	def __getitem__(self, key):
		return self.map[key]
	
	def __setitem__(self, key, value):
		self.map[key] = value
	
	def __delitem__(self, key):
		sys.stderr.write("Cannot edit single fingerprints.\n")
		
	def __repr__(self):
		return self.name
		
	def keys(self):
		return self.fingerprints()
	
	def fingerprints(self):
		return [run for run in self.map.values() if isinstance(run, fingerprint)]
	
	def features(self):
		try:
		    return [self.map[feat] for feat in labels]
		except KeyError:
		    return [feat for feat in self.map.values() if isinstance(feat, feature)]

	def parse_tab(self, tab):
		"""Function takes an open tab file of CP data and inputs it into the internal
		dictionaries."""
		runs = []
		for line in tab:
			if "DMSO" in line.strip().upper():
				continue
			elif "FEATURES" in line.strip().upper():
				features = [feature(feat) for feat in line.strip().upper().split("\t")[1:]]
			else:
				finger = line.strip().split("\t")
				runs.append(fingerprint(finger[0].rsplit("_", 1)[0]))
#                                runs.append(fingerprint(finger[0]))
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
		"""Function takes a file object name and writes the cp data matrix to that open file as a 
		tab delimited text file.
		"""
		f = open(tab, 'w')
		runs = sorted(self.fingerprints(), key=lambda x: str(x))
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
		
# 	def parallel_pearson(self, cpus=None):
# 		"""Function returns an nxn heatmap of pearson correlations with idruns as rows and columns
# 		as columns. The tuple that is returned is (nxn, labels)
# 		"""
# 		idruns = self.fingerprints()
# 		self.labels = [id.name for id in idruns]
# 		
# 		ppservers = ()
# 		job_server = pp.Server(cpus, ppservers=ppservers, secret='acetone') if cpus else pp.Server(ppservers=ppservers, secret='acetone')
# 		
# 		self.nxn = numpy.zeros((len(idruns), len(idruns)))
# 		jobs = numpy.empty((len(idruns), len(idruns)), dtype=object)
# 		
# 		for ind1, id1 in enumerate(idruns):
# 			for ind2, id2 in enumerate(idruns):
# 				jobs[(ind1, ind2)] = job_server.submit(pearsonr, (id1.values(), id2.values()), (), ("scipy.stats",))
# 		
# 		idruns = self.fingerprints()
# 		self.labels = [id.name for id in idruns]
# 		
# 		widgets = ['Pearson Correlations: ', Percentage(), ' ', Bar(marker=RotatingMarker()), ' ',
# 				ETA(), ' ', FileTransferSpeed()]
# 		pbar = ProgressBar(widgets=widgets, maxval=len(idruns)).start()
# 		
# 		self.nxn = numpy.zeros((len(idruns), len(idruns)))
# 		for ind1, id1 in enumerate(idruns):
# 			for ind2, id2 in enumerate(idruns):
# 				p = jobs[(ind1, ind2)]()
# 				self.nxn[(ind1,ind2)] = p
# 				self[id1.name].setp(id2.name, p)
# 			pbar.update(ind1)
# 		pbar.finish()
# 		self.nxn = self.nxn.astype(numpy.float)

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
				
	def get_ancestor(self, f1, f2, name):
		"""Function accepts two fingerprint objects and the name of the ancestor, returns the
		average fingerprint of the two (the ancestor), and then sets the parent of both fingerprint
		objects to the ancestor."""
		a = fingerprint(name)
		for feat in set(f1.keys()).union(set(f2.keys())):
			a[feat] = numpy.mean([f1[feat], f2[feat]])
		a.children = (f1, f2)
		f1.parent = a
		f2.parent = a
		return a
		
	def h_cluster(self, ids, number=0):
		if len(ids) == 1:
			return ids.values()[0]
		pair = []
		largest = -1.0
		for id in ids:
			for ppair in ids[id].pearson.items():
				if ppair[1] > largest and not ppair[0] == id:
					largest = ppair[1]
					pair = [id, ppair[0]]
				else:
					continue
		a = self.get_ancestor(ids[pair[0]], ids[pair[1]], "ancestor_{}".format(number))
		ids = dict([(id, ids[id]) for id in ids if id not in pair])
		for id in ids.values():
			id.delp(pair[0])
			id.delp(pair[1])
			p = pearsonr(id.values(), a.values())[0]
			a.setp(id.name, p)
			id.setp(a.name, p)
		ids[a.name] = a
		return self.h_cluster(ids, number + 1)
				
	def hierarchical(self):
		#A list of idrun/fingerprint objects
		idruns = self.fingerprints()
		idruns = dict([(id.name, id) for id in idruns])
		a0 = self.h_cluster(idruns)
		net = nx.Graph()
		self.network(a0, net)
		return a0, net
				
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
		
		
	def cp_network(self, min_sim = 0.9):
		h = nx.Graph()
		for id in self.fingerprints():
			h.add_node(id.name, weight=id.activity_score())
		for id1, id2 in itertools.combinations(self.fingerprints(), 2):
			if self.get_pearson(id1.name, id2.name) >= min_sim:
				h.add_edge(id1.name, id2.name, weight=self.get_pearson(id1.name, id2.name))
		return h
		
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
		x = range(len(y))
		f, (ax1, ax2) = plt.subplots(1, 2, sharey=False)
		ax1.scatter(x, y)
		ax2.hist(y, 30)
		plt.savefig(param + "_plot.pdf")
		plt.close()
		for pair in sorted(self[param].items(), key=lambda p:p[1]):
			print "{}, {:4f}".format(*pair)
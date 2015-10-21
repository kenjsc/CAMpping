
from progressbar import ProgressBar, Percentage, Bar, RotatingMarker, ETA, FileTransferSpeed
import itertools
import copy
import btools
import numpy

import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt

import sys

"""Module contains classes defining objects for parsing and indexing mass spectral basketing 
data.
"""



class idrun:
	def __init__(self, name):
		self.name = name
		self.peaks = dict()

	def __getitem__(self, baskname):
		return self.peaks[baskname]
	
	def __setitem__(self, baskname, peak):
		self.peaks[baskname] = peak
		
	def __delitem__(self, baskname):
		del self.peaks[baskname]

	def values(self):
		return self.peaks.values()
	
	def keys(self):
		return self.peaks.keys()
		
	def items(self):
		return self.peaks.items()
		
	def __repr__(self):
		return self.name


class basket:
	"""Class defines a basket for binning adducts."""
	def __init__(self, name):
		self.mode = None
		self.name = name
		self.peaks = dict()
		self.fingerprint = dict()
		self.m = 0
		self.rt = 0
		self.activity = None
		self.cluster_score = None
		self.ions = set()

	def __getitem__(self, run):
		return self.peaks[str(run)]
	
	def __setitem__(self, run, peak):
		self.m = ((self.m * len(self.peaks)) + peak.isotopes[0].mz) / (len(self.peaks) + 1)
		self.rt = ((self.rt * len(self.peaks)) + peak.isotopes[0].rt) / (len(self.peaks) + 1)
		self.ions.add(peak.isotopes[0].ion)
		self.peaks[str(run)] = peak
		if not self.mode:
			self.mode = str(peak.mode)
	
	def __delitem__(self, run):
		self.m *= len(self.peaks)
		self.rt *= len(self.peaks)
		run = self.peaks.pop(str(run))
		self.m -= run.isotopes[0].mz
		self.rt -= run.isotopes[0].rt
		self.ions = set([ad.isotopes[0].ion for ad in self.peaks.values()])
		if len(self.peaks) > 0:
			self.m /= len(self.peaks)
			self.rt /= len(self.peaks)
		else:
			self.m = 0
			self.rt = 0
			self.mode = None
			
	def __repr__(self):
#		return "{} : m/z={:4f}, rt={:2f}, AAS={:1f}".format(",".join(self.ions), self.m, self.rt,
#				self.activity)
		return "basket={}, m/z={:4f}, rt={:2f}; ions={}; scores={}_{}".format(self.name, self.m, self.rt, ",".join(self.ions), self.activity, self.cluster_score)
		
	def values(self):
		return self.peaks.values()
		
	def keys(self):
		return self.peaks.keys()
		
	def items(self):
		return self.peaks.items()
		
	def iso_pat(self):
		pat = []
		for group in itertools.izip_longest(*[ad.get_pat() for ad in self.values()]):
			mzs = [t[0] for t in group if t]
			absos = [t[1] for t in group if t]
			pat.append((numpy.mean(mzs), numpy.mean(absos)))
		return pat

		
	def param_score(self, parameter, scores):
		self.fingerprint[parameter] = scores
		
	def finprint(self, order = None):
		if not order:
			order = self.fingerprint.keys()
		return [self.fingerprint[param][0] for param in order]

class heatmap:
	"""Class defines a mass spectral binning heatmap."""
	def __init__(self, *peaks):
		self.map = dict()
		if not len(peaks):
			return
		for bask in self.basketer(7, 0.2, 0.5, *peaks):
			self.map[bask.name] = bask
			for idname, peak in bask.items():
				if idname not in self.map:
					self.map[idname] = idrun(idname)
				self.map[idname][bask.name] = peak

	def __getitem__(self, key):
		return self.map[key]
	
	def __setitem__(self, key, value):
		sys.stderr.write("Cannot insert single runs/baskets\n")
	
	def __delitem__(self, key):
		del self.map[key]
		
	def keys(self):
		return self.idruns()
	
	def idruns(self):
		return [run for run in self.map.values() if isinstance(run, idrun)]
	
	def baskets(self):
		return [bask for bask in self.map.values() if isinstance(bask, basket)]
		
	def basketer(self, ppm, drt, diso, *peaks):
		"""Function takes a delta rt and ppm, and any number of adducts and bins adducts that are
		within the delta rt and ppm differences. Returns a list of basket objects."""

		peaks = sorted(peaks, key=lambda ad: ad.isotopes[0].mz)

		widgets = ['BASKETING: ', Percentage(), ' ', Bar(marker=RotatingMarker()), ' ', ETA(),\
				' ', FileTransferSpeed()]
		pbar = ProgressBar(widgets=widgets, maxval=len(peaks)).start()
	
		baskets = []
		inpro = []
		n = 0
		basknum = 0
		for peak in peaks:
			match = False
			for indbin, in_bask in enumerate(inpro):
				ppmdif = abs(peak.isotopes[0].mz - in_bask.m) / in_bask.m * 1000000
				rtdif = abs(peak.isotopes[0].rt - in_bask.rt)
				if ppmdif <= ppm and rtdif <= drt and str(peak.mode) == in_bask.mode and \
						btools.iso_diff(peak.get_pat(), in_bask.iso_pat()) >= diso:
					in_bask[peak.idrun] = peak
					match = True
					break
				elif ppmdif > ppm:
					baskets.append(in_bask)
					inpro[indbin] = None

			k = 0
			while k < len(inpro) and len(inpro):
				if inpro[k]:
					k += 1
					continue
				else:
					inpro.pop(k)

			if not match:
				inpro.append(basket(basknum))
				inpro[-1][peak.idrun] = peak
				basknum += 1
			n += 1
			pbar.update(n)
		baskets.extend(inpro)
		pbar.finish()
		return baskets
		
	def remove_basket(self, bname):
		"""Function accepts the name of a basket to be removed from the heatmap. Removes all links
		to that basket object."""
		if bname in self.map:
			for run in self.map[bname].keys():
				del self.map[run][bname]
			del self.map[bname]
	
	def remove_run_complete(self, runname):
		"""Function accepts the name of an idrun to be removed from the heatmap. Removes the idrun
		from the heatmap and removes all baskets that were contained in that idrun.
		"""
		if runname in self.map:
			for bask in self.map[runname].keys():
				self.remove_basket(bask)
			del self.map[runname]
				
	def remove_run_baskets(self, runname):
		"""Function accepts the name of an idrun to have its baskets removed. Removes all baskets
		that include that idrun, but does not remove the idrun itself. Run remove_empty_runs to do
		that.
		"""
		if runname in self.map:
			for bask in self.map[runname].keys():
				self.remove_basket(bask)
			
	def remove_run(self, runname):
		"""Function accepts the name of an idrun to be removed from the heatmap. Removes all links
		to the idrun from the heatmap, but does not remove the baskets that the idrun contributes
		to.
		"""
		if runname in self.map:
			for bask in self.map[runname].keys():
				del self.map[bask][runname]
			del self.map[runname]
				
	def remove_empty_baskets(self):
		"""Function checks through the heatmap and removes baskets that don't have links to any
		idruns.
		"""
		for bask in self.baskets():
			if len(bask.keys()) <= 0:
				del self.map[bask.name]
				
	def remove_empty_runs(self):
		"""Function checks through the heatmap and removes idruns that don't have links to any
		baskets.
		"""
		for run in self.idruns():
			if len(run.keys()) <= 0:
				del self.map[run.name]
		
	def find_mz(self, mz, ppm=50):
		"""Function accepts a mass and a ppm value and prints a list of baskets in order of mz,
		which are present in the heatmap with mass of mz +- ppm.
		"""
		basks = []
		for bask in self.baskets():
			if btools.ppm(bask.m, mz) < ppm:
				basks.append(bask)
		return sorted(basks, key=lambda x: x.m)
		
	def output_tab(self, file, restrict=False):
		"""Function takes an open file and outputs the baskets to a .tab file in the same format as
		the tab CP file that was loaded in. Rows are baskets instead of idruns and columns are the
		same features. File can be clustered with Cluster3.0 and viewed in Java TreeView.
		"""
		basks = self.baskets()
                labels = "MICRO_AND_BINUCLEATED_CELLS_MICRONUCLEI_EDU	MICRO_AND_PROBEABPOSITIVE_CELLS_MICRONUCLEI_EDU	MICRO_AND_PROBEBPOSITIVE_CELLS_MICRONUCLEI_EDU	MICRO_AND_PROBEAPOSITIVE_CELLS_MICRONUCLEI_EDU	INTERPHASE_CELLS_MICRONUCLEI_EDU	MONONUCLEATED_CELLS_MICRONUCLEI_EDU	BINUCLEATED_CELLS_MICRONUCLEI_EDU	MICRONUCLEATED_CELLS_MICRONUCLEI_EDU	CELLS_WITH_ONE_MICRONUCLEUS_MICRONUCLEI_EDU	PCT_BINUCLEATED_CELLS_MICRONUCLEI_EDU	NUCLEAR_DIVISION_INDEX_MICRONUCLEI_EDU	PCT_CELLS_WITH_ONE_MICRONUCLEUS_MICRONUCLEI_CYTO	PCT_MICRONUCLEATED_CELLS_MICRONUCLEI_CYTO	MICRONUCLEI_PER_CELL_MICRONUCLEI_CYTO	TOTAL_MICRONUCLEI_MICRONUCLEI_CYTO	ALL_ACTIN_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_ACTIN_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_CYTO	ALL_ACTIN_MEAN_STAIN_INTEGR_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_ACTIN_MEAN_STAIN_INTEGR_INTENS_MULTIWAVESCORING_CYTO	ALL_ACTIN_MEAN_STAIN_AREA_MULTIWAVESCORING_CYTO	POSITIVE_ACTIN_MEAN_STAIN_AREA_MULTIWAVESCORING_CYTO	POSITIVE_TUBULIN_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_CYTO	ALL_TUBULIN_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_CYTO	CELL_TUBULIN_NUCLEUS_INTEGR_INTENSITY_MULTIWAVESCORING_CYTO	DNA_MEAN_AREA_MICRONUCLEI_EDU	BREADTH_IMA_SUMMARY_CYTO	INNER_RADIUS_IMA_SUMMARY_CYTO	OUTER_RADIUS_IMA_SUMMARY_CYTO	LENGTH_IMA_SUMMARY_CYTO	PERIMETER_IMA_SUMMARY_CYTO	MEAN_RADIUS_IMA_SUMMARY_CYTO	EQUIV_RADIUS_IMA_SUMMARY_CYTO	EQUIV_PROLATE_VOL_IMA_SUMMARY_CYTO	ALL_NUCLEI_MEAN_AREA_MULTIWAVESCORING_CYTO	AREA_IMA_SUMMARY_CYTO	EQUIV_SPHERE_SURFACE_AREA_IMA_SUMMARY_CYTO	TOTAL_AREA_IMA_SUMMARY_CYTO	PIXEL_AREA_IMA_SUMMARY_CYTO	EQUIV_SPHERE_VOL_IMA_SUMMARY_CYTO	EQUIV_OBLATE_VOL_IMA_SUMMARY_CYTO	NEGATIVE_EDU_MEAN_CELL_INTEGR_INTENS_MULTIWAVESCORING_EDU	NEGATIVE_EDU_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_EDU	OUTER_RADIUS_IMA_SUMMARY_EDU	LENGTH_IMA_SUMMARY_EDU	EQUIV_SPHERE_VOL_IMA_SUMMARY_EDU	EQUIV_OBLATE_VOL_IMA_SUMMARY_EDU	EQUIV_PROLATE_VOL_IMA_SUMMARY_EDU	DAPI_STAINED_AREA_MULTIWAVESCORING_EDU	CELL_TOTAL_AREA_MULTIWAVESCORING_EDU	WIDTH_IMA_SUMMARY_EDU	PERIMETER_IMA_SUMMARY_EDU	EQUIV_RADIUS_IMA_SUMMARY_EDU	MEAN_RADIUS_IMA_SUMMARY_EDU	ALL_CELLS_MEAN_AREA_MULTIWAVESCORING_EDU	ALL_NUCLEI_MEAN_AREA_MULTIWAVESCORING_EDU	PIXEL_AREA_IMA_SUMMARY_EDU	EQUIV_SPHERE_SURFACE_AREA_IMA_SUMMARY_EDU	AREA_IMA_SUMMARY_EDU	TOTAL_AREA_IMA_SUMMARY_EDU	NUCLEAR_AREA_PER_CELL_TRANSFLUOR_EDU	INNER_RADIUS_IMA_SUMMARY_EDU	BREADTH_IMA_SUMMARY_EDU	ALL_PH3_MEAN_CELL_AVER_INTENS_MULTIWAVESCORING_EDU	ALL_PH3_MEAN_NUCL_AVER_INTENS_MULTIWAVESCORING_EDU	ALL_PH3_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_EDU	ALL_PH3_MEAN_CELL_INTEGR_INTENS_MULTIWAVESCORING_EDU	NEGATIVE_PH3_MEAN_CELL_INTEGR_INTENS_MULTIWAVESCORING_EDU	NEGATIVE_PH3_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_EDU	PCT_CELLS_WITH_MULTI_MICRONUCLEI_MICRONUCLEI_EDU	ALL_EDU_MEAN_STAIN_AREA_MULTIWAVESCORING_EDU	PCT_POSITIVE_EDU_MULTIWAVESCORING_EDU	CELL_PH3_STAINED_AVERAGE_INTENSITY_MULTIWAVESCORING_CYTO	CELL_PH3_STAINED_AVERAGE_INTENSITY_MULTIWAVESCORING_EDU	NEGATIVE_PH3_MEAN_STAIN_AREA_MULTIWAVESCORING_EDU	NEGATIVE_PH3_MEAN_STAIN_INTEGR_INTENS_MULTIWAVESCORING_EDU	POSITIVE_ACTIN_MEAN_CELL_AVER_INTENS_MULTIWAVESCORING_CYTO	ALL_ACTIN_MEAN_CELL_AVER_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_ACTIN_MEAN_CYTO_AVER_INTENS_MULTIWAVESCORING_CYTO	ALL_ACTIN_MEAN_CYTO_AVER_INTENS_MULTIWAVESCORING_CYTO	PIT_INTEGRATED_INTENSITY_TRANSFLUOR_CYTO	PIT_COUNT_TRANSFLUOR_CYTO	PIT_TOTAL_AREA_TRANSFLUOR_CYTO	GRADIENT_INDEX_TRANSFLUOR_CYTO	LAPLACIAN_INDEX_TRANSFLUOR_CYTO	POSITIVE_TUBULIN_MEAN_CYTO_AVER_INTENS_MULTIWAVESCORING_CYTO	ALL_TUBULIN_MEAN_CYTO_AVER_INTENS_MULTIWAVESCORING_CYTO	ALL_TUBULIN_MEAN_STAIN_AVER_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_TUBULIN_MEAN_STAIN_AVER_INTENS_MULTIWAVESCORING_CYTO	ALL_TUBULIN_MEAN_CELL_AVER_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_TUBULIN_MEAN_CELL_AVER_INTENS_MULTIWAVESCORING_CYTO	NEGATIVE_ACTIN_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_CYTO	NEGATIVE_ACTIN_MEAN_STAIN_AVER_INTENS_MULTIWAVESCORING_CYTO	NEGATIVE_ACTIN_MEAN_CYTO_AVER_INTENS_MULTIWAVESCORING_CYTO	NEGATIVE_ACTIN_MEAN_CELL_AVER_INTENS_MULTIWAVESCORING_CYTO	NEGATIVE_ACTIN_MEAN_NUCL_AVER_INTENS_MULTIWAVESCORING_CYTO	MITOTIC_CELLS_MICRONUCLEI_CYTO	TOTAL_CELLS_MICRONUCLEI_CYTO	TOTAL_CELLS_MULTIWAVESCORING_CYTO	POSITIVE_TUBULIN_MULTIWAVESCORING_CYTO	POSITIVE_ACTIN_MULTIWAVESCORING_CYTO	DNA_TOTAL_AREA_MICRONUCLEI_CYTO	INTERPHASE_CELLS_MICRONUCLEI_CYTO	MONONUCLEATED_CELLS_MICRONUCLEI_CYTO	MITOTIC_CELLS_MICRONUCLEI_EDU	POSITIVE_EDU_MULTIWAVESCORING_EDU	PIT_TOTAL_AREA_TRANSFLUOR_EDU	PIT_COUNT_TRANSFLUOR_EDU	NUCLEAR_INTEGRATED_INTENSITY_TRANSFLUOR_EDU	NUCLEAR_TOTAL_AREA_TRANSFLUOR_EDU	DNA_TOTAL_AREA_MICRONUCLEI_EDU	TOTAL_CELLS_MICRONUCLEI_EDU	TOTAL_CELLS_MULTIWAVESCORING_EDU	NUCLEAR_COUNT_TRANSFLUOR_EDU	ELL_FORM_FACTOR_IMA_SUMMARY_CYTO	ELL_FORM_FACTOR_IMA_SUMMARY_EDU	NEGATIVE_PH3_MEAN_STAIN_AVER_INTENS_MULTIWAVESCORING_EDU	ALL_PH3_MEAN_STAIN_INTEGR_INTENS_MULTIWAVESCORING_EDU	ALL_PH3_MEAN_STAIN_AREA_MULTIWAVESCORING_EDU	PCT_POSITIVE_PH3_MULTIWAVESCORING_EDU	POSITIVE_PH3_MULTIWAVESCORING_EDU	POSITIVE_PH3_MULTIWAVESCORING_CYTO	PCT_POSITIVE_PH3_MULTIWAVESCORING_CYTO	ALL_PH3_MEAN_STAIN_AREA_MULTIWAVESCORING_CYTO	ALL_PH3_MEAN_STAIN_INTEGR_INTENS_MULTIWAVESCORING_CYTO	NEGATIVE_PH3_MEAN_STAIN_AREA_MULTIWAVESCORING_CYTO	NEGATIVE_PH3_MEAN_STAIN_INTEGR_INTENS_MULTIWAVESCORING_CYTO	PIT_AVERAGE_INTENSITY_TRANSFLUOR_EDU	NEGATIVE_EDU_MEAN_STAIN_AVER_INTENS_MULTIWAVESCORING_EDU	CELL_PROBEA_AVERAGE_INTENSITY_MICRONUCLEI_EDU	CELL_EDU_NUCLEUS_AVERAGE_INTENSITY_MULTIWAVESCORING_EDU	CELL_EDU_CELL_AVERAGE_INTENSITY_MULTIWAVESCORING_EDU	ALL_EDU_MEAN_CELL_INTEGR_INTENS_MULTIWAVESCORING_EDU	ALL_EDU_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_EDU	ALL_EDU_MEAN_STAIN_INTEGR_INTENS_MULTIWAVESCORING_EDU	POSITIVE_EDU_MEAN_CELL_AVER_INTENS_MULTIWAVESCORING_EDU	POSITIVE_EDU_MEAN_NUCL_AVER_INTENS_MULTIWAVESCORING_EDU	ALL_EDU_MEAN_STAIN_AVER_INTENS_MULTIWAVESCORING_EDU	POSITIVE_EDU_MEAN_STAIN_AVER_INTENS_MULTIWAVESCORING_EDU	POSITIVE_EDU_MEAN_STAIN_INTEGR_INTENS_MULTIWAVESCORING_EDU	POSITIVE_EDU_MEAN_CELL_INTEGR_INTENS_MULTIWAVESCORING_EDU	POSITIVE_EDU_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_EDU	CELLULAR_GRADIENT_INDEX_TRANSFLUOR_EDU	CELLULAR_LAPLACIAN_INDEX_TRANSFLUOR_EDU	PIT_INTEGRATED_INTENSITY_TRANSFLUOR_EDU	ALL_EDU_MEAN_NUCL_AVER_INTENS_MULTIWAVESCORING_EDU	ALL_EDU_MEAN_CELL_AVER_INTENS_MULTIWAVESCORING_EDU	CELLULAR_TEXTURE_INDEX_TRANSFLUOR_EDU	LAPLACIAN_INDEX_TRANSFLUOR_EDU	GRADIENT_INDEX_TRANSFLUOR_EDU	TEXTURE_INDEX_TRANSFLUOR_EDU	CELL_EDU_STAINED_AVERAGE_INTENSITY_MULTIWAVESCORING_EDU	POSITIVE_EDU_MEAN_STAIN_AREA_MULTIWAVESCORING_EDU	CELL_EDU_STAINED_AREA_MULTIWAVESCORING_EDU	POSITIVE_PH3_MEAN_STAIN_AREA_MULTIWAVESCORING_EDU	NEGATIVE_PH3_MEAN_STAIN_AVER_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_PH3_MEAN_CYTO_AVER_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_PH3_MEAN_CELL_AVER_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_PH3_MEAN_STAIN_AREA_MULTIWAVESCORING_CYTO	ALL_PH3_MEAN_STAIN_AVER_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_PH3_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_PH3_MEAN_STAIN_INTEGR_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_PH3_MEAN_STAIN_AVER_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_PH3_MEAN_NUCL_AVER_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_PH3_MEAN_CELL_INTEGR_INTENS_MULTIWAVESCORING_EDU	POSITIVE_PH3_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_EDU	POSITIVE_PH3_MEAN_STAIN_INTEGR_INTENS_MULTIWAVESCORING_EDU	POSITIVE_PH3_MEAN_STAIN_AVER_INTENS_MULTIWAVESCORING_EDU	POSITIVE_PH3_MEAN_NUCL_AVER_INTENS_MULTIWAVESCORING_EDU	POSITIVE_PH3_MEAN_CELL_AVER_INTENS_MULTIWAVESCORING_EDU	ALL_PH3_MEAN_STAIN_AVER_INTENS_MULTIWAVESCORING_EDU	NEGATIVE_TUBULIN_MEAN_NUCL_AVER_INTENS_MULTIWAVESCORING_CYTO	NEGATIVE_TUBULIN_MEAN_CELL_AVER_INTENS_MULTIWAVESCORING_CYTO	NEGATIVE_TUBULIN_MEAN_CYTO_AVER_INTENS_MULTIWAVESCORING_CYTO	NEGATIVE_TUBULIN_MEAN_STAIN_INTEGR_INTENS_MULTIWAVESCORING_CYTO	NEGATIVE_TUBULIN_MEAN_STAIN_AREA_MULTIWAVESCORING_CYTO	NEGATIVE_TUBULIN_MEAN_STAIN_AVER_INTENS_MULTIWAVESCORING_CYTO	NEGATIVE_TUBULIN_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_CYTO	NEGATIVE_TUBULIN_MEAN_CELL_INTEGR_INTENS_MULTIWAVESCORING_CYTO	PCT_POSITIVE_TUBULIN_MULTIWAVESCORING_CYTO	BINUCLEATED_CELLS_MICRONUCLEI_CYTO	MULTINUCLEATED_CELLS_MICRONUCLEI_CYTO	MULTIMONO_CELL_RATIO_MICRONUCLEI_CYTO	MULTIDUALMONO_CELL_RATIO_MICRONUCLEI_CYTO	DUALMONO_CELL_RATIO_MICRONUCLEI_CYTO	PCT_BINUCLEATED_CELLS_MICRONUCLEI_CYTO	MULTIDUALMONO_CELL_RATIO_MICRONUCLEI_EDU	DUALMONO_CELL_RATIO_MICRONUCLEI_EDU	DNA_MEAN_AREA_MICRONUCLEI_CYTO	ALL_PH3_MEAN_NUCL_INTEGR_INTENS_MULTIWAVESCORING_CYTO	CELL_ACTIN_CYTOPLASM_AVERAGE_INTENSITY_MULTIWAVESCORING_CYTO	CELL_ACTIN_CELL_AVERAGE_INTENSITY_MULTIWAVESCORING_CYTO	NEGATIVE_ACTIN_MEAN_STAIN_INTEGR_INTENS_MULTIWAVESCORING_CYTO	NEGATIVE_ACTIN_MEAN_STAIN_AREA_MULTIWAVESCORING_CYTO	PCT_CELLS_WITH_MULTI_MICRONUCLEI_MICRONUCLEI_CYTO	SHAPE_FACTOR_IMA_SUMMARY_CYTO	PCT_PROBEA_POSITIVE_CELLS_MICRONUCLEI_EDU	PCT_PROBEB_POSITIVE_CELLS_MICRONUCLEI_EDU	PCT_PROBEAB_POSITIVE_CELLS_MICRONUCLEI_EDU	PCT_CELLS_WITH_ONE_MICRONUCLEUS_MICRONUCLEI_EDU	PCT_MICRONUCLEATED_CELLS_MICRONUCLEI_EDU	MICRONUCLEI_PER_CELL_MICRONUCLEI_EDU	MICRONUCLEI_PER_HEALTHY_CELL_MICRONUCLEI_EDU	MICRONUCLEI_PER_MONONUCLEATED_CELL_MICRONUCLEI_CYTO	MICRONUCLEI_PER_PROBEBPOSITIVE_CELL_MICRONUCLEI_EDU	MICRONUCLEI_PER_PROBEABPOSITIVE_CELL_MICRONUCLEI_EDU	MICRONUCLEI_PER_PROBEAPOSITIVE_CELL_MICRONUCLEI_EDU	MICRONUCLEI_PER_MONONUCLEATED_CELL_MICRONUCLEI_EDU	NEGATIVE_EDU_MEAN_STAIN_AREA_MULTIWAVESCORING_EDU	NEGATIVE_EDU_MEAN_STAIN_INTEGR_INTENS_MULTIWAVESCORING_EDU	CELLS_WITH_MULTI_MICRONUCLEI_MICRONUCLEI_CYTO	CELL_ACTIN_NUCLEUS_AVERAGE_INTENSITY_MULTIWAVESCORING_CYTO	CELL_ACTIN_STAINED_AVERAGE_INTENSITY_MULTIWAVESCORING_CYTO	TEXTURE_INDEX_TRANSFLUOR_CYTO	PIT_AVERAGE_INTENSITY_TRANSFLUOR_CYTO	CELL_TUBULIN_CELL_AVERAGE_INTENSITY_MULTIWAVESCORING_CYTO	ALL_TUBULIN_MEAN_NUCL_AVER_INTENS_MULTIWAVESCORING_CYTO	POSITIVE_TUBULIN_MEAN_NUCL_AVER_INTENS_MULTIWAVESCORING_CYTO	CELLULAR_GRADIENT_INDEX_TRANSFLUOR_CYTO	CELL_GRADIENT_INDEX_TRANSFLUOR_CYTO	CELL_LAPLACIAN_INDEX_TRANSFLUOR_CYTO	CELL_TEXTURE_INDEX_TRANSFLUOR_CYTO	CELL_TUBULIN_NUCLEUS_AVERAGE_INTENSITY_MULTIWAVESCORING_CYTO	CELL_PIT_AVERAGE_INTENSITY_TRANSFLUOR_CYTO	PIT_COUNT_PER_CELL_TRANSFLUOR_CYTO	PIT_AREA_PER_CELL_TRANSFLUOR_CYTO	POSITIVE_ACTIN_MEAN_NUCL_AVER_INTENS_MULTIWAVESCORING_CYTO	MULTINUCLEATED_CELLS_MICRONUCLEI_EDU	MULTIMONO_CELL_RATIO_MICRONUCLEI_EDU	PCT_MULTINUCLEATED_CELLS_MICRONUCLEI_EDU	CELLS_WITH_MULTI_MICRONUCLEI_MICRONUCLEI_EDU	ALL_NUCLEI_MEAN_INTEGR_ITENSITY_MULTIWAVESCORING_EDU	TOTAL_INTENSITY_IMA_SUMMARY_EDU	ALL_NUCLEI_MEAN_INTEGR_ITENSITY_MULTIWAVESCORING_CYTO	TOTAL_INTENSITY_IMA_SUMMARY_CYTO	PCT_MITOTIC_CELLS_MICRONUCLEI_EDU	ALL_NUCLEI_MEAN_AVERAGE_ITENSITY_MULTIWAVESCORING_CYTO	AVERAGE_INTENSITY_IMA_SUMMARY_CYTO	CELL_DNA_AVERAGE_INTENSITY_MICRONUCLEI_CYTO	CELL_NUCLEAR_AVERAGE_INTENSITY_MICRONUCLEI_CYTO	DAPI_AVERAGE_INTENSITY_MULTIWAVESCORING_CYTO	CELL_NUCLEAR_AVERAGE_INTENSITY_TRANSFLUOR_CYTO	CELL_DNA_AVERAGE_INTENSITY_MICRONUCLEI_EDU	CELL_NUCLEAR_AVERAGE_INTENSITY_MICRONUCLEI_EDU	CELL_NUCLEAR_AVERAGE_INTENSITY_TRANSFLUOR_EDU	DAPI_AVERAGE_INTENSITY_MULTIWAVESCORING_EDU	NUCLEAR_AVERAGE_INTENSITY_TRANSFLUOR_EDU	ALL_NUCLEI_MEAN_AVERAGE_ITENSITY_MULTIWAVESCORING_EDU	AVERAGE_INTENSITY_IMA_SUMMARY_EDU".split("\t")
		file.write("Features\t{}\n".format("\t".join(labels)))
		for bask in basks:
			if restrict and bask.fingerprint["NUCLEAR_COUNT_TRANSFLUOR_EDU"][0] > -0.6:
				file.write("{}_{}\t{}\n".format(bask.m, bask.rt, "\t".join([str(bask.fingerprint[param][0])
						for param in labels])))
			elif restrict:
				continue
			else:	
				file.write("{}\t{}\n".format(str(bask), "\t".join([str(bask.fingerprint[param][0])
						for param in labels])))
					
	def grab_basks(self, run):
		if run not in self.map:
			return []
		return [self.map[bask] for bask in self.map[run].keys()]

					
	def subheat(self, runs, antiruns=[]):
		"""Function accepts an idrun, extracts all the baskets that are in that idrun, and creates
		a sub heatmap composed of only those baskets (and the idruns that compose them).
		"""
		antiruns = set([run for run in antiruns if run in self.map])
		baskets = set([bask for run in runs for bask in self.map[run].keys() if not len(set(self.map[bask].keys()) & antiruns)])
		subheatmap = heatmap()
		sub = dict()
		for bask in baskets:
			sub[bask] = copy.copy(self.map[bask])
			sub[bask].peaks = copy.copy(self.map[bask].peaks)
			for inner_run in self.map[bask].keys():
				if inner_run not in sub:
					sub[inner_run] = idrun(inner_run)
				sub[inner_run][bask] = self.map[bask][inner_run]
		subheatmap.map = sub
		return subheatmap
		
	def subheat_restrict(self, runs, antiruns=[]):
		"""
		"""
		antiruns = set([run for run in antiruns if run in self.map])
		antibaskets = set([bask for run in antiruns for bask in self.map[run].keys()])
		subheatmap = heatmap()
		sub = dict()
		
		for run in runs:
			if run not in self.map:
				continue
			sub[run] = copy.copy(self.map[run])
			sub[run].peaks = copy.copy(self.map[run].peaks)
			for bask in sub[run].keys():
				if bask in antibaskets:
					del sub[run][bask]
					continue
				if bask not in sub:
					sub[bask] = basket(bask)
				sub[bask][run] = self.map[run][bask]
		subheatmap.map = sub
		subheatmap.remove_empty_runs()
		for bask in subheatmap.baskets():
			bask.activity = self.map[bask.name].activity
			bask.cluster_score = self.map[bask.name].cluster_score
			bask.fingerprint = self.map[bask.name].fingerprint
		return subheatmap
		
	def remove_undercount(self, count):
		"""Function accepts a minimum count value, and removes all baskets from the heatmap
		that contain fewer than that number of adducts.
		"""
		for bask in self.baskets():
			if len(bask.keys()) <= count:
				self.remove_basket(bask.name)
		
	def export_nxn(self, idrun_list=[]):
		"""Function returns an nxn heatmap of retention times with idruns as rows and baskets
		as columns. The tuple that is returned is (nxn, row labels, column labels)
		"""
		if idrun_list:
			idruns = idrun_list
		else:
			idruns = [id.name for id in self.idruns()]
		baskets = sorted(self.baskets(), key=lambda bask: bask.m)
		
		nxn = numpy.zeros((len(idruns), len(baskets)))
		for iind, idrun_name in enumerate(idruns):
			for bind, bask in enumerate(baskets):
				if idrun_name in self.map:
					if bask.name in self.map[idrun_name].peaks:
						nxn[(iind, bind)] = bask.rt
		return nxn, idruns, [str(bask) for bask in baskets]
		
	def isotope_plot(self, name):
		isos = []
		for b in self.baskets():
			for comb in itertools.combinations(b.values(), 2):
				isos.append(btools.iso_diff(comb[0].get_pat(), comb[1].get_pat()))
		plt.hist(isos, bins=1000, range=(0.0, 1.0))
		plt.savefig(name + "_isotope_score.pdf")
		plt.close()
		
	def plot_cluster_score(self, name):
		basks = self.baskets()
		x = range(len(basks))
		y = sorted([bask.cluster_score for bask in basks])
		f, (ax1, ax2) = plt.subplots(1, 2, sharey=False)
		ax1.scatter(x, y)
		ax2.hist(y, 30)
		plt.savefig(name + "_cluster_score.pdf")
		plt.close()
		
	def plot_activity(self, name):
		basks = self.baskets()
		x = range(len(basks))
		y = sorted([bask.activity for bask in basks])
		f, (ax1, ax2) = plt.subplots(1, 2, sharey=False)
		ax1.scatter(x, y)
		ax2.hist(y, 30)
		plt.savefig(name + "_activity.pdf")
		plt.close()

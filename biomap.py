"""Function takes a text file of dilution data from HiTS to perform the
BioMAP algorithm and assign BioMAP profiles to individual samples. The input data
is tab delimited text file where the order of the columns does not matter, but the 
column headings must be consistent. In particular, these columns are absolutely necessary:

	HAW_COORIDNATES
	IC_ALIAS_ID
	HA_SDESC
	IW_FINAL_DILUTION
	HAR_TIMEPOINT_MS
	HAR_VALUE_DOU

The HAR values must be comma delimeted. Furthermore, the negative control (no drug)
must be in the 24th column. All data to be analyzed must be between columns 3 and 22 
(inclusive). Both the concentration and the sample name must be included for each sample
with the exception of the controls. The arrangement of the data in between columns 3 and 
22 does not matter, but it is typically a dilution series down the column.
"""

import argparse
import pyeq2

import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt

import math
import numpy as np
import sys

import os.path



def score(absorbances):
	return max(absorbances) - min(absorbances)

def four_param(x, a, b, c, d):
	return ((a - d) / (1 + (x / c) ** b)) + d

def fp_fit(x, y):
	"""Given a list of x and y values, function fits a four parameter logistic 
	distribution to the data. The distribution being fit has bounds set on several of
	the parameters to keep the distribution in the proper orientation. These are:
	
	a	-0.25	0.25
	b	-inf	-0.1
	c	0		inf
	d	0.75	1.25
	
	The function returns a tuple of the parameters and the covariance of the fit:
	
	((a, b, c, d), cov)
	
	"""
	equation = pyeq2.Models_2D.Sigmoidal.FourParameterLogistic()

	data = "\n".join("{} {}".format(x1, y1) for x1, y1 in zip(x, y))
	
	equation.upperCoefficientBounds = [0.25, -0.1, None, 1.25]
	equation.lowerCoefficientBounds = [-0.25, None, 0, 0.75]
	
	pyeq2.dataConvertorService().ConvertAndSortColumnarASCII(data, equation, False)
	equation.Solve()
	
	return equation.solvedCoefficients, equation.CalculateAllDataFittingTarget(equation.solvedCoefficients)
				
class dilution:
	"""Class defines the dilution data for a given prefraction with a given bacteria"""
	def __init__(self):
		self.dilutions = dict()
		self.ic = None
		self.mic = None
		self.nmic = None
		
	def __getitem__(self, dil):
		"""Function returns the fold change value for a given dilution"""
		return self.dilutions[dil]
		
	def __setitem__(self, dil, fold):
		"""Function assigns a fold change value to the given dilution"""
		self.dilutions[dil] = fold
		
	def __delitem(self, dil):
		"""Function deletes the specified dilution value from the object"""
		del self.dilutions[dil]
	
	def values(self):
		"""function returns the values of the underlying dilutions dictionary"""
		return self.dilutions.values()
	
	def keys(self):
		"""Function returns the keys of the underlying dilutions dictionary"""
		return self.dilutions.values()
		
	def items(self):
		"""Function returns the items of the underlying dilutions dictionary"""
		return self.dilutions.items()
		
	def scale(self, value):
		"""Function scales all of the absorbances in the dilution dictionary by the 
		specified value (usually the average of the controls"""
		for val in self.dilutions:
			self.dilutions[val] /= value
		
	def fit(self, pdf=None, redo=False):
		"""Function fits a four parameter logistic curve to the dilution series and
		saves the fit parameters, the ic50, and the covariance. Optional pdf argument
		specifies to output the graph to the pdf file specified. The output plots
		the observed data as blue circles, the ic50 as a green pentagon, and the fit
		data as red crosses.
		"""
		x, y = zip(*self.items())
		if not self.ic or redo:
			vals, cov = fp_fit(x, y)
			self.ic = vals[2]
			self.vals = vals
			self.cov = cov
		if not pdf == None:
			yn = [four_param(xi, *self.vals) for xi in x]
			plt.scatter(x, yn, c='r', marker="+")
			plt.scatter(x, y)
			plt.plot([self.ic], [four_param(self.ic, *self.vals)], c='g', marker="p")
			plt.xscale('log')
			try:
				plt.savefig(pdf + "dilutions.pdf")
			except ValueError:
				pass
			plt.close()
			
	def calc_mic(self, pdf=None):
		"""Function calculates the mic of the dilution series using the curve fit. Calls
		the curve fit function if it has not already been called.
		"""
		if not self.ic:
			self.fit(pdf)
		conc = sorted(self.dilutions.keys(), reverse=True)
		if self.cov > 0.5:
			self.mic = -1
		elif len([ab for ab in self.dilutions.values() if ab < 0.6]) == 0:
			self.mic = -2
		elif len([ab for ab in self.dilutions.values() if ab > 0.3]) == 0:
			self.mic = -3
		else:
			n = 1
			while (conc[n] > self.ic) and (n < len(conc)):
				n += 1
			self.mic = conc[n-1]
		return self.mic
	

class bacterium:
	"""Class defines a bacterium in BioMAP. """
	def __init__(self, name):
		self.name = name
		self.prefractions = dict()

	def __getitem__(self, pref):
		return self.prefractions[pref]
	
	def __setitem__(self, pref, dil):
		"""Function assigns a dilution object to a prefraction in the prefractions dictionary"""
		self.prefractions[pref] = dil
	
	def __delitem__(self, pref):
		del self.prefractions[pref]
		
	def __repr__(self):
		return self.name
		
	def values(self):
		return self.prefractions.values()
		
	def keys(self):
		return self.prefractions.keys()
		
	def items(self):
		return self.prefractions.items()
		
		
		
class prefraction:
	"""Class defines a BioMAP fingerprint for one compound/well."""
	def __init__(self, name):
		self.name = name
		self.bacteria = dict()

	def __getitem__(self, bact):
		return self.bacteria[bact]
	
	def __setitem__(self, bact, dil):
		self.bacteria[bact] = dil
		
	def __delitem__(self, bact):
		del self.bacteria[bact]
		
	def __repr__(self):
		return self.name
		
	def values(self):
		return self.bacteria.values()
	
	def keys(self):
		return self.bacteria.keys()
		
	def items(self):
		return self.bacteria.items()
	
	def calc_profile(self, err=sys.stderr, pdf=None):
		"""Function calculates the biomap profile for a prefraction based on all of the
		bacteria that have been treated with this compound dilution series. The 
		normalized MIC value is written to nmic in each of the bacteria. To access:
		
		prefraction[bacteria].nmic
		
		where you replace bacteria with the bacterium of interest
		"""
		largest = 0
		for bac in self.bacteria:
			if self.bacteria[bac].calc_mic("{}_{}_{}_".format(pdf, str(self), str(bac))) > largest:
				largest = self.bacteria[bac].mic
		largestnmic = 0
		for bac in self.bacteria:
			dil = self.bacteria[bac]
			
			if dil.mic == -1:
				err.write("Sample: {}; Organism: {} - Poor Fit.\n".format(str(self), str(bac)))
			elif dil.mic == -2:
				err.write("Sample: {}; Organism: {} - Too dilute. Never reached MIC.\n".format(str(self), str(bac)))
			elif dil.mic == -3:
				err.write("Sample: {}; Organism: {} - Too concentrated. Dilute and re-screen.\n".format(str(self), str(bac)))
			
			if dil.mic < 0:
				dil.nmic = 0
			else:
				dil.nmic = math.log10(10 * ((dil.mic / largest) ** (-1)))
			
			if dil.nmic > largestnmic:
				largestnmic = dil.nmic
				
		for bac in self.bacteria:
			if largestnmic == 0:
				break
			self.bacteria[bac].nmic /= largestnmic
	
	def plot_profile(self, dir="./"):
		fig = plt.figure(figsize=(8, len(self.bacteria) * 0.5))
		ax = fig.add_subplot(111)

		pos = np.arange(len(self.bacteria)) + 0.5    # Center bars on the Y-axis ticks
	
		ax.barh(pos, [dil.nmic for dil in self.bacteria.values()], align='center', height=0.5,)
		
		plt.yticks(pos, self.bacteria.keys())
		plt.xlabel('Normalized MIC')
		plt.title(str(self))
		
		fig.tight_layout(rect=[0, 0, 1, 1])
		
		plt.savefig("{}_{}_profile.pdf".format(dir, str(self)))
		plt.close()
		
class biomap:
	"""Class defines a BioMAP heatmap. Because data is normalized to controls, this should
	only be used on a per plate basis."""
	def __init__(self, biofile):
		self.headings = dict()
		self.map = dict()
		self.parse_txt(biofile)
		self.control_scale()
	
	def __getitem__(self, key):
		return self.map[key]
	
	def __setitem__(self, key, value):
		self.map[key] = value
	
	def __delitem__(self, key):
		sys.stderr.write("Cannot edit single fingerprints.\n")
		
	def prefractions(self):
		return [pref for pref in self.map.values() if isinstance(pref, prefraction) and str(pref) != 'control']
	
	def bacteria(self):
		return [bac for bac in self.map.values() if isinstance(bac, bacterium)]

	def keys(self):
		return self.prefractions()
		
	def calc_profiles(self, err=sys.stderr, pdf=None):
		"""Function calculates all of the biomap profiles for each of the compounds/
		prefractions in the heatmap.
		"""
		for pref in self.prefractions():
			pref.calc_profile(err, pdf)
		
	def write_tab(self, output=sys.stdout):
		"""Function accepts an open file object (default is stdout) and writes the
		biomap profiles of all of the prefractions/compounds to that file.
		"""
		bacs = map(str, self.bacteria())
		output.write("Name\t{}\n".format("\t".join(bacs)))
		for pref in self.prefractions():
			output.write("{}\t{}\n".format(str(pref), "\t".join(map(str, [pref[bac].nmic for bac in bacs]))))
		output.close()
    
	def parse_txt(self, file, delimeter="\t"):
		"""Reads a tab delimited text file which has been read in as a string from HiTS Biomap. 
		Creates the dictionary grid of bacterium/prefraction dictionaries.
		"""
		file = file.split("\n")
		for ind, field in enumerate(file[0].strip("\n").split(delimeter)):
			self.headings[field] = ind
		for line in file[1:]:
			if line == "":
				continue
			line = line.strip("\n").split(delimeter)
			well = line[self.headings["HAW_COORDINATES"]]
			if "02" in well or "23" in well or "01" in well:
				continue
			pref = line[self.headings["IC_ALIAS_ID"]]
			if "24" in well:
				pref = "control"
			bac = line[self.headings["HA_SDESC"]]
			try:
				dil = float(line[self.headings["IW_FINAL_DILUTION"]])
			except ValueError:
				dil = 0
			if "24" in well:
				dil = 2 ** ord(well[0]) - 65
			timepoints = map(float, line[self.headings["HAR_TIMEPOINT_MS"]].strip("\"").split(","))
			abs = map(float, line[self.headings["HAR_VALUE_DOU"]].strip("\"").split(","))
			
			self.map.setdefault(bac, bacterium(bac))
			self.map.setdefault(pref, prefraction(pref))
			
			if pref not in self.map[bac].prefractions:
				self.map[bac][pref] = self.map[pref][bac] = dilution()
			
			self.map[bac][pref][dil] = score(abs)
		
	def control_scale(self):
		"""Scales each bacterial dataset by its control values.
		"""
		for bac in self.bacteria():
			val = np.mean(self.map["control"][str(bac)].values())
			for pref in self.map[str(bac)].keys():
				self.map[str(bac)][pref].scale(val)
			
	def export_profiles(self, dir="./"):
		for pref in self.prefractions():
			pref.plot_profile(dir)

def parse_args(args):
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
			description=__doc__)
	
	parser.add_argument('--infile', '-f', dest='infile', type=argparse.FileType('U'),
			default=sys.stdin, help='plate reader output from the necessary bacteria. Pulled from database.')
	parser.add_argument('--out', '-o', dest='out', type=argparse.FileType('w'),
			default=sys.stdout, help='output filename for generated tab file. Default is stdout')
	parser.add_argument('--outic', '-i', dest='ic', default="",
			help='output directory for generating IC50 curvefits. Default is None')
	parser.add_argument('--outp', '-p', dest='prof', default="",
			help='output directory for generating profile bar charts. Default is None')
	parser.add_argument('--log', '-l', dest='log', type=argparse.FileType('w'), 
			default=sys.stderr, help='output filename for any errors. Please read the\
			log file after running the program!')
		
	args = parser.parse_args(args)
	return args
	
def write_plates(output, *plates):
	"""
	"""
	bacs = set(map(str, plates[0].bacteria()))
	for plate in plates[1:]:
		bacs = bacs & set(map(str, plate.bacteria()))
	bacs = list(bacs)
	output.write("Name\t{}\n".format("\t".join(bacs)))
	for plate in plates:
		for pref in plate.prefractions():
			output.write("{}\t{}\n".format(str(pref), "\t".join(map(str, [pref[bac].nmic for bac in bacs]))))
	output.close()
	
def split_plates(txt):
	labels = txt.readline().strip("\n")
	label_dict = dict()
	for lind, label in enumerate(labels.split("\t")):
		label_dict[label] = lind
	plates = dict()
	for line in txt:
		line.strip("\n")
		name = line.split("\t")[label_dict["HAP_BARCODE"]]
		plates.setdefault(name, [labels])
		plates[name].append(line)
	for plate in plates:
		plates[plate] = "\n".join(plates[plate])
	print plates
	return plates

		
def main(args):
	"""Main function. Exports cdt, gradient, and binary heatmaps."""
	options = parse_args(args)
	processed = []
	for plate, text in split_plates(options.infile).items():
		processed.append(biomap(text))
		if options.ic:
			options.ic = "{}/{}".format(*os.path.split(options.ic))
			processed[-1].calc_profiles(options.log, options.ic)
		else:
			processed[-1].calc_profiles(options.log)
		if options.prof:
			options.prof = "{}/{}".format(*os.path.split(options.prof))
			processed[-1].export_profiles(options.prof)
	write_plates(options.out, *processed)
	options.infile.close()
	options.log.close()

if __name__ == "__main__":
	try:
		sys.exit(main(sys.argv[1:]))
	except EnvironmentError as (errno,strerr):
		sys.stderr.write("ERROR: " + strerr + "\n")
		sys.exit(errno)
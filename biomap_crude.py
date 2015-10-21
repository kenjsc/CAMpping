"""
Command line utility for generating heatmap of hit/no-hit data for each sample against
several bacterial strains. The bacteria are columns along the x-axis in the heatmaps,
and the drug samples are rows along the y-axis. The hit/no-hit is denoted by color,
where red indicates hit.

The input data is a tab delimeted text file, with comma delimeted absorbance data as
one of the columns. While the order of the columns does not matter, the column
headings are very important. The important headings are:

IC_ALIAS_ID
HA_SDESC
HAR_TIMEPOINT_MS
HAR_VALUE_DOU

These are ascribed in the HiTS database system for QB3.
"""


import argparse
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt




def parse_file(ofile, delimeter='\t'):
	"""Function iterates through lines in a file, stripping whitespace and splitting
	by a defined delimeter (default is tab). Returns an iterator per line.
	"""
	for line in ofile:
		yield line.strip().split(delimeter)
		
def score(timepoints, absorbances):
	"""Function returns a score or fold-change value that describes how much the bacteria
	has grown. Currently defined as max absorbance - min absorbance. This means a higher
	value corresponds to more growth, and a lower value corresponds to more inhibition.
	"""
	return max(absorbances) - min(absorbances)
	# / (timepoints[-1] - timepoints[0]))

def get_data(ofile, delimeter='\t'):
	"""Pulls data from the provided open file object and returns it as a dictionary
	matrix with wells and bacteria as keys.
	"""
	data = dict()
	plate = ""
	labels = dict()
	for ind, label in enumerate(ofile.readline().strip().split(delimeter)):
		labels[label] = ind
	scores = dict()
	for line in parse_file(ofile, delimeter):
		well, bacteria, timepoints, absorbances = (line[labels['IC_ALIAS_ID']],
			line[labels['HA_SDESC']], map(float, line[labels['HAR_TIMEPOINT_MS']].strip('"').split(',')),
			map(float, line[labels['HAR_VALUE_DOU']].strip('"').split(',')))
		sc = score(timepoints, absorbances)
		data.setdefault(well, dict())[bacteria] = sc
		scores.setdefault(bacteria, [1, 0])
		if sc > scores[bacteria][1]:
			scores[bacteria][1] = sc
		if sc < scores[bacteria][0]:
			scores[bacteria][0] = sc
	for well in data:
		for bac in data[well]:
			data[well][bac] = (data[well][bac]) / (scores[bac][1])
	return data
		
def export_cdt(data, output):
	"""Exports the same data observed in the heatmap to a cdt file that can be opened
	with java TreeView.
	"""
	bacteria = sorted(data.values()[0].keys())
	output.write("GID\tFeatures\tNAME\tGWEIGHT\t{}\n".format("\t".join(bacteria)))
	output.write("AID\t\t\t\t{}\n".format("\t".join(["ARRY{}X".format(i) for i in range(len(bacteria))])))
	output.write("EWEIGHT\t\t\t\t{}\n".format("\t".join(['1.0'] * len(bacteria))))
	n = 0
	for well in sorted(data.keys(), key=lambda w: sum([d**11 for d in data[w].values()])):
		output.write("GENE{0}X\t{1}_{2}\t{1}_{2}\t1.0\t{3}\n".format(n, well,
			sum([d**3 for d in data[well].values()]), "\t".join([str(data[well][bac]) for bac in bacteria])))
	output.close()
	

def export_pdf(data, output, gradient = True):
	"""Exports the data as a heatmap to the file specified by output (.pdf). The gradient
	marker specifies if the color in the heatmap should be a gradient that shows the
	fold change score, or as binary red/white colors with a cutoff for maximum fold-change
	score to be considered a hit.
	"""
	bacteria = sorted(data.values()[0].keys())
	matrix = np.zeros((len(data), len(bacteria)))
	rows = sorted(data.keys(), key=lambda w: np.sum([d**3 for d in data[w].values()]) / (len([d for d in data[w].values() if d < 0.3]) + 1))
	for wind, well in enumerate(rows):
		for bind, bacterium in enumerate(bacteria):
			if gradient:
				matrix[(wind, bind)] = data[well][bacterium]
			else:
				matrix[(wind, bind)] = 0 if data[well][bacterium] < 0.3 else 0.5
	plt.matshow(matrix, cmap=plt.get_cmap('RdYlGn'), vmin=0, vmax=1)
	plt.xticks(range(len(bacteria)), bacteria, rotation=90)
	plt.yticks(range(len(rows)), rows)
	ax = plt.gca()
	for posi in ax.spines:
		ax.spines[posi].set_color('none')
	ax.tick_params(labelcolor='k', top='off', bottom='off', left='off', right='off')
	fig = plt.gcf()
	fig.set_size_inches(10, 100)
	plt.savefig(output + ".pdf", bbox_inches='tight', dpi=200)
	plt.close()


def parse_args(args):
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
			description=__doc__)
	
	parser.add_argument('--infile', '-f', dest='infile', type=argparse.FileType('U'),
			default=sys.stdin, help='plate reader output from the necessary bacteria. Pulled from database.')
	parser.add_argument('--outcdt', '-o', dest='out', type=argparse.FileType('w'),
			default=sys.stdout, help='output filename for generated cdt. Default is stdout')
	parser.add_argument('--outpdf', '-p', dest='pdf', default='output',
			help='output filename for generated pdf. Default is output')
	parser.add_argument('--delimeter', '-d', dest='delimeter', default='\t', 
			help="""delimeter used in file. Default is tab. Delimeter cannot be comma
			because that is what is the delimeter for the absorbance/timestamp data.""")
		
	args = parser.parse_args(args)
	return args
		
def main(args):
	"""Main function. Exports cdt, gradient, and binary heatmaps."""
	options = parse_args(args)
	data = get_data(options.infile, options.delimeter)
	export_cdt(data, options.out)
	export_pdf(data, options.pdf + "_gradient")
	export_pdf(data, options.pdf + "_binary", False)
	options.infile.close()
	options.out.close()


if __name__ == "__main__":
	try:
		sys.exit(main(sys.argv[1:]))
	except EnvironmentError as (errno,strerr):
		sys.stderr.write("ERROR: " + strerr + "\n")
		sys.exit(errno)

#This is a command line utility written in python to take in a tab delimeted table of values
#with column headings and x values as the first column. 

import argparse
import sys
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
import numpy as np


def parse_args(args):
	parser = argparse.ArgumentParser()
	
	parser.add_argument("--file", "-f", dest="file", nargs="?", type=
			argparse.FileType('r'), default=sys.stdin, help=
			"Input file for plotting. Default is stdin.")
	parser.add_argument("--l_title", "-l", action="store", dest="l", default="", help=
			"Legend title for the plot. Default is no title.")
	parser.add_argument("--title", "-t", action="store", dest="t", default="", help=
			"Plot title. Default is not title.")
	parser.add_argument("--l_position", "-p", action="store", dest="p", type=int, default=4, help=
			"Location for legend. Options are: 1 (upper right), 2 (upper left), 3 (lower left), \
			4 (lower right), 5 (right), 6 (center left), 7 (center right), 8 (lower center), \
			9 (upper center), and 10 (center).")
	parser.add_argument("--output", "-o", action="store", dest="out", default="plot",
			help="File output")
			
	
	return parser.parse_args(args)
	
	
def plot(file, output, title=""):
	
	data = []
	
	for line in file:
		if line.startswith("#"):
			continue
		data.append(line.split(","))

	for ind, curve in enumerate(data):
		print ind, curve
		data[ind] = (curve[0], map(np.float, curve[1:]))
			
	fig = plt.figure(figsize=(4, 3))
	ax = fig.add_subplot(111)
	ax.set_xlabel("Iteration Number")
	ax.set_ylabel("Likelihood Score")
	for ind, curve in enumerate(sorted(data, key=lambda tup: tup[0])):
		ax.plot(range(1, len(curve[1]) + 1), curve[1], label=curve[0], lw=1,\
				linestyle="-", color = plt.get_cmap('spectral')(float(ind)/(len(data) + 1)))
	ax.set_yscale('log')
#	leg = ax.legend(title="", loc=1, fancybox=True)
	#leg.get_frame().set_alpha(0.3)
	ax.set_title(title)
	plt.savefig(output + ".pdf", bbox_inches='tight')
	
def main(args):
	options = parse_args(args)
	plot(options.file, options.out, options.t)


if __name__ == "__main__":
	main(sys.argv[1:])
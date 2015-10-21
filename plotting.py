
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
	parser.add_argument("--yaxis", "-y", action="store", dest="y", default="", help=
			"Y-axis label for the plot. Default is no label. X-axis label is from file.")
	parser.add_argument("--l_title", "-l", action="store", dest="l", default="", help=
			"Legend title for the plot. Default is no title.")
	parser.add_argument("--title", "-t", action="store", dest="t", default="", help=
			"Plot title. Default is not title.")
	parser.add_argument("--l_position", "-p", action="store", dest="p", type=int, default=4, help=
			"Location for legend. Options are: 1 (upper right), 2 (upper left), 3 (lower left), \
			4 (lower right), 5 (right), 6 (center left), 7 (center right), 8 (lower center), \
			9 (upper center), and 10 (center).")
			
	
	return parser.parse_args(args)
	
	
def main(args):
	options = parse_args(args)
	
	data = []
	
	for line in options.file:
		if line.startswith("#"):
			continue
		data.append(line.split("\t"))
		
	av_data = dict()
	er_data = dict()
	x_label = None
	x = []
	for curve in zip(*data):
		if not x:
			x_label = curve[0]
			x = map(float, curve[1:])
			continue
		if " ".join(curve[0].split(" ")[:-1]) not in av_data:
			av_data[" ".join(curve[0].split(" ")[:-1])] = map(float, curve[1:])
		else:
			er_data[" ".join(curve[0].split(" ")[:-1])] = [np.std([one, two]) for one, two in zip(av_data[" ".join(curve[0].split(" ")[:-1])], map(float, curve[1:]))]
			av_data[" ".join(curve[0].split(" ")[:-1])] = [np.mean([one, two]) for one, two in zip(av_data[" ".join(curve[0].split(" ")[:-1])], map(float, curve[1:]))]
	for curve in av_data:
		if curve in er_data:
			continue
		else:
			er_data[curve] = [0.0] * len(av_data[curve])
			
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_xlabel(x_label)
	ax.set_ylabel(options.y)
	for curve in sorted(av_data.keys(), key=lambda name: name.split(" ")[::-1]):
		if curve.startswith("X"):
			continue
		ax.errorbar(x, av_data[curve], yerr=er_data[curve], label=curve, marker="o", lw=2, linestyle="-")
	leg = ax.legend(title=options.l, loc=options.p, fancybox=True)
	leg.get_frame().set_alpha(0.3)
	ax.set_title(options.t)
	plt.savefig(options.file.name.split(".")[0] + ".pdf")


if __name__ == "__main__":
	main(sys.argv[1:])
"""Command line tool for importing cef files to a ms SQLite database. The
syntax is straightforward, use the help associated with the argparse module.
"""

from sqlalchemy.orm import sessionmaker
import argparse
import sys
import pp
import dbtables
import val
import cefparse
import btools
from progressbar import ProgressBar, Percentage, Bar, RotatingMarker, ETA, FileTransferSpeed
import re
import itertools

def run_name(cef):
	"""Function accepts a filename for a cef file and picks out the unique country code, extract,
	fraction, minute motif. Function returns a string of this motif (e.g. RLUS-1647E-74 or
	RLUS-1674E).
	"""
	meta = re.search(r"RL[A-Z]{2}\-\d{4}[A-F]\-\d{1,2}", cef)
	if meta:
		return meta.group()
	meta = re.search(r"RL[A-Z]{2}\-\d{4}[A-F]", cef)
	if meta:
		return meta.group()
	return None

def group_cefs(*cefs):
	"""Function accepts a number of .cef files and organizes them by positive/negative mode and
	4GHz/2GHz mode. If the number of .cef files supplied is not divisable by 4 (i.e. there are not
	the same number of POS_2G, POS_4G, NEG_2G, and NEG_4G files) an error is raised. Function
	returns a dictionary with keys: pos_2g, pos_4g, neg_2g, and neg_4g and values: list of cef
	file names in order of 
	"""
	#dictionary with the four modes as keys and a list of cef file names
	# as values in its respective mode
	ceflists = dict([('pos_2g', []), ('pos_4g', []), ('neg_2g', []), ('neg_4g',[])])
	
	for cef in cefs:
		upcef = cef.upper()
		if "POS" in upcef and "2G" in upcef:
			ceflists['pos_2g'].append(cef)
		elif "POS" in upcef and "4G" in upcef:
			ceflists['pos_4g'].append(cef)
		elif "NEG" in upcef and "2G" in upcef:
			ceflists['neg_2g'].append(cef)
		elif "NEG" in upcef and "4G" in upcef:
			ceflists['neg_4g'].append(cef)
	
	length = len(ceflists.values()[0])
	for ceflist in ceflists:
		ceflists[ceflist] = sorted(ceflists[ceflist], key=run_name)
		if not len(ceflists[ceflist]) == length:
			print "Unequal lengths of modes"
			for group in itertools.izip_longest(*ceflists.values()):
				print group
			sys.exit(1)
	return ceflists

def process(run_cefs, blank_cefs, cutoff, rtcut):
	"""Function accepts 2 dictionaries. One of actual run cef files which are mode-file NAME
	pairs, while the other is blank (noise) cef files which are mode-IDrun object pairs. Function
	subtracts the blanks, removes saturated peaks, validates the peaks between 4GHz and 2GHz mode,
	combines positive and negative data, and writes information to disk. Note the modes have to
	match between the two cef file dictionaries with pos_2g, pos_4g, neg_2g, and neg_4g as keys.
	"""
	for mode in run_cefs:
		run_cefs[mode] = cefparse.read_cef(run_cefs[mode])
		btools.rmblank(run_cefs[mode], blank_cefs[mode])
	return val.validate(cutoff, rtcut, *run_cefs.values())

def parse_args(args):
	"""Function takes a list of arguments (sys.args passed in via command line) and
	parses them. This parser accepts flags:
		--files, -f
		--blanks, -b
		--database, -d
	"""
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,\
			description=__doc__)
	parser.add_argument("--files", "-f", dest="cefs", nargs="+", help="""The .cef files to import
			into the database. There must be positive/negative mode and 4GHz/2GHz mode data for each
			idrun, and the location code, extract, fraction, and minute motif must be the same for
			all datasets for one idrun.""")
	parser.add_argument("--blanks", "-b", dest="blanks", nargs="+", help="""The .cef files to
			subtract as blanks from the other cef files. The same mode restrictions apply.""")
	parser.add_argument("--db", "-d", dest="db", action="store", help="The database to write the\
			mass spectral data")
	parser.add_argument("--cpus", "-n", dest="cpus", action="store", type=int, default=None, help=\
			"The number of CPU cores to use parallely when processing data. Default is all that are\
			available.")
	parser.add_argument("--cutoff", "-c", dest="cutoff", action="store", type=int, default=1000,\
			help="Default minimum abundance cutoff for compounds, default is 1000.")
	parser.add_argument("--minrt", "-rt", dest="rt", action="store", type=float, default=0.35,\
			help="Minimum retention time for compounds to avoid solvent front. Default is 0.35")
	args = parser.parse_args(args)
	return args


def main(args):
	
	options = parse_args(args)
	
	#Open up connection to database
	engine = dbtables.connect(options.db)
	Session = sessionmaker(bind=engine)
	session = Session()

        #Get all of the blanks in groups, as a dictionary by mode. Then compress to a single run
        # for each mode. Note: the compound information in compressed runs is not very reliable.
        # the adduct information is what is preserved. 
	blanks = group_cefs(*options.blanks)
	for mode in blanks:
		for ind, blank in enumerate(blanks[mode]):
			blanks[mode][ind] = cefparse.read_cef(blank)
		blanks[mode] = btools.comb_runs(20, 0.4, *blanks[mode])
	runs = []
	ceflists = group_cefs(*options.cefs)	
	modes = ceflists.keys()
	
	#This is for the parallelization, this says only use this computer, no network cluster
	ppservers = ()
	if options.cpus:
		job_server = pp.Server(options.cpus, ppservers=ppservers, secret='acetone')
	else:
		job_server = pp.Server(ppservers=ppservers, secret='acetone')
	print "Running with ", job_server.get_ncpus(), " CPU's"
	
	widgets = ['Submitting Jobs: ', Percentage(), ' ', Bar(marker=RotatingMarker()), ' ', ETA(),\
			' ', FileTransferSpeed()]
	pbar = ProgressBar(widgets=widgets, maxval=len(ceflists.values()[0])).start()
	for ind, group in enumerate(zip(*ceflists.values())):
		if len(group) < 4:
			print "Skipping: ", group
			continue
		small = dict([(mode, cef) for mode, cef in zip(modes, group)])
		#Submit the group as a job to parallel process
		
		print group

		runs.append(job_server.submit(process, (small, blanks, options.cutoff, options.rt,), (),
				("btools", "sys", "re", "cefparse", "sqlalchemy.orm", "sqlalchemy.ext.declarative",
				"dbtables", "val",)))
		pbar.update(ind + 1)
	pbar.finish()

	#This portion is for the progress bar.
	widgets = ['COMMITING: ', Percentage(), ' ', Bar(marker=RotatingMarker()), ' ', ETA(),\
			' ', FileTransferSpeed()]
	pbar = ProgressBar(widgets=widgets, maxval=len(runs)).start()

	for ind, run in enumerate(runs):
		session.add(run())
		session.commit()
		pbar.update(ind + 1)
	pbar.finish()
	
if __name__ == "__main__":
	main(sys.argv[1:])
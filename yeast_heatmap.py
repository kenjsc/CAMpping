"""
For converting .csv files downloaded from HiTS into heatmaps for cluster 3.0

Cannot include the absorption column, because that is a time series with more
commas.

"""



import itertools
from scipy.stats import johnsonsu, kstest

def parse_data(filename):
	"""Function accepts a filename that represents yeast data downloaded from hits
	and outputs a dictionary of dictionary of dictionary that can be indexed by
	prefraction, strain, and metric (in any order) to get the raw value.
	
	ex. data[pref][strain][metric] = data[strain][pref][metric] = etc... = value
	
	Function also returns the prefractions, strains, and sets that make up the
	keys in the data dictionary
	
	"""
	#This will be the dictionary of data that can be indexed by strain, metric, or 
	# prefraction. It is three levels of dictionary, with 6 different routes to the
	# data:
	# prefraction --> metric --> strain --> value
	# prefraction --> strain --> metric --> value
	# metric --> prefraction --> strain --> value
	# metric --> strain --> prefraction --> value
	# strain --> prefraction --> metric --> value
	# strain --> metric --> prefraction --> value
	# ex. data_dict[prefraction][metric][strain] will give the value for that
	# prefraction, metric and strain (as will the other 5 combinations of that query)
	data_dict = dict()
	
	#because all of these keys will be intermixed and we won't be able to tell
	# strain from prefraction from metric, we need to keep track:
	prefractions = set()
	strains = set()
	metrics = set()
	
	with open(filename, "U") as f:
		strain_line = f.readline().strip().split(",")
		metric_line = f.readline().strip().split(",")
		#This is a list of labels, which is a merge between row 1 and row 2 in the hits
		# data, it will be strain_metric, in the order of appearance in the data file.
		labels = []
		#strain keeps track of the strain in row 1 for the metrics in line 2 until the
		# next strain in line 1 is reached
		strain = ""
		#keeps track of the number of columns before the data is actually reached, this
		# way we can skip to that position when we're looking at the data
		start = 0
		for strn, met in zip(strain_line, metric_line):
			if not strain:
				start += 1
			if strn:
				strain = strn
			labels.append("{}_{}".format(strain, met))
	
		data = [[column.strip() for column in line.strip().split(",")] for line in f]
		
		#index in each row of the prefraction (the prefraction column number)
		pref_index = labels.index("_Catalog Number")
		

		
		for row in data:
			#get the prefraction name from the data row
			pref = row[pref_index]
			prefractions.add(pref)
			data_dict.setdefault(pref, dict())
			for label, value in zip(labels[start:], row[start:]):
				strain, metric = label.split("_")
				strains.add(strain)
				metrics.add(metric)
				data_dict[pref].setdefault(strain, dict())
				data_dict[pref].setdefault(metric, dict())
				data_dict.setdefault(strain, dict())
				data_dict[strain].setdefault(pref, dict())
				data_dict[strain].setdefault(metric, dict())
				data_dict.setdefault(metric, dict())
				data_dict[metric].setdefault(pref, dict())
				data_dict[metric].setdefault(strain, dict())
				if value == "":
					value = 0
				else:
					value = float(value)
					
				data_dict[pref][metric][strain] = data_dict[pref][strain][metric]\
						= data_dict[strain][pref][metric]\
						= data_dict[strain][metric][pref]\
						= data_dict[metric][strain][pref]\
						= data_dict[metric][pref][strain] = value

	return data_dict, prefractions, strains, metrics
		
def remove_dimension(data, metric, prefractions, strains):
	"""Function accepts a data dictionary that is three nested dictionaries, with each
	nesting being a dimension (explained more in parse_data function). Function also takes
	the name of the metric to flatten around, and the lists of keys for the prefractions
	and strains
	ex. 
	
	one dimension is strains
	another is prefractions
	another is metrics
	
	we are only interested in the metric absorptiondiff, so we include absorptiondiff
	as the metric and the lists of strains and prefractions.
	
	The dimension of metrics is lost, and all values become absorptiondiff, so:
	
	data[strain][prefraction][metric] = value becomes data[strain][prefraction] = absorptiondiff value
	"""
	
	new_data = dict()
	for pref, strn in itertools.product(prefractions, strains):
		new_data.setdefault(pref, dict())
		new_data.setdefault(strn, dict())
		new_data[pref][strn] = new_data[strn][pref] = data[pref][strn][metric]
	return new_data
			
def curvefit(series):
	"""Function accepts an iterable series of datapoints representing a johnson su
	distribution. The series is fit, tested for goodness of fit (p value greater than 0.1
	and data mean value greater than 0.5) and the mean - 1/2 the std deviation is
	returned. None is returned if fit is bad.
	"""
	parameters = johnsonsu.fit(series)
	func = johnsonsu(*parameters)
	D, p = kstest(series, 'johnsonsu', N = len(series), args=parameters)
	if p < 0.1 or func.mean() < 0.5:
		return None
	return func.mean() - func.std() / 2
	
def binarify_dataset(data, strains):
	"""Function accepts a flattened data dictionary generated by the remove_dimension and
	parse functions, a list of strains of interest. Function fits the data series for each
	strain to a johnsonsu distribution, identifies a live/dead cutoff value and the
	mean - 1/2 the std deviation, and replaces all 'live' values with 0 and all 'dead'
	values with 1.
	"""
	binary_data = dict()
	good_strains = []
	for strn in strains:
		cutoff = curvefit(data[strn].values())
		if not cutoff:
			continue
		binary_data.setdefault(strn, dict())
		for pref in data[strn].keys():
			binary_data[strn].setdefault(pref, dict())
			binary_data.setdefault(pref, dict())
			binary_data[pref].setdefault(strn, dict())
			if data[strn][pref] >= cutoff:
				binary_data[strn][pref] = binary_data[pref][strn] = 0
			else:
				binary_data[strn][pref] = binary_data[pref][strn] = 1
		good_strains.append(strn)
	return binary_data, good_strains


def write_heatmap(filename, data, prefractions, strains):
	"""Function accepts a filename, a data dictionary and lists of prefractions and
	strains, and writes a file with prefraction rows and strain columns with the value
	in the corresponding row/column.
	"""
	with open(filename, 'w') as f:
		f.write("Prefractions\t{}\n".format("\t".join(strains)))
		for pref in prefractions:
			f.write("{}\t{}\n".format(pref, "\t".join([str(data[pref][strn]) for strn in strains])))
	return
	
def main(infile, outfile, metric):
	"""Function accepts a filename with yeast data and a metric of interest from the data,
	parses the data, and writes a heatmap of analyzed data for that metric.
	"""
	
	raw_data, prefractions, strains, metrics = parse_data(infile)
	data = remove_dimension(raw_data, metric, prefractions, strains)
	binary_data, good_strains = binarify_dataset(data, strains)
	write_heatmap(outfile, binary_data, prefractions, good_strains)

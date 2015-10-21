
import numpy as np
import scipy.stats.stats as st

def get_dmsos(cp):
	"""function accepts a cp object and returns a list of fingerprint objects with
	'DMSO' in their title.
	"""
	return [finger for finger in cp.fingerprints() if "DMSO" in str(finger)]
	
def dmso_means(cp):
	"""Function accepts a cp object and returns a mean and standard deviation of all 
	parameters for the dmso fingerprints in this form:
	
	list of (feature, mu, std) for each feature
	"""
	dmsos = [dmso.items() for dmso in get_dmsos(cp)]
	data = []
	for parameter_pair in zip(*dmsos):
		parameter, values = zip(*parameter_pair)
		data.append((parameter[0], st.nanmean(values), st.nanstd(values)))
	return data
	
def all_means(cp):
	"""Function accepts a cp object and returns a mean and standard deviation of all
	parameters for all fingerprints in this form:
	list of (feature, mu, std) for each feature
	"""
	data = []
	for parameter_pair in zip(*[finger.items() for finger in cp.fingerprints()]):
		parameter, values = zip(*parameter_pair)
		print parameter, values
		print st.nanmean(values), st.nanstd(values)
		data.append([parameter[0], st.nanmean(values), st.nanstd(values)])
	return data
	
def batch_correct(cp, std=0.2):
	"""Function accepts a cp object and a standard deviation value. Function edits the
	cp object values so that the DMSO fingerprint values have a mean of 0 and standard
	deviation of supplied value. All other fingerprints are altered in an identical 
	fashion. (Editing is actually done on a per parameter basis and includes subtracting
	the offset of the DMSO mean values, and scaling by the supplied std value over the 
	observed std value).
	"""
	for parameter, mu, obs_std in dmso_means(cp):
		for finger in cp.fingerprints():
#			if obs_std == 0:
#				new = finger[parameter] - mu
#			else:
#				new = (finger[parameter] - mu) * (std/obs_std)
			new = finger[parameter] - mu
			finger[parameter] = cp[parameter][str(finger)] = new
			
def plate_batch(*cps):
	"""Function accepts a broken out list of cp objects (any number of comma separated
	cp objects), calculates their mean and standard deviation for each parameter, and 
	normalizes the objects together. Writes a file of all cp objects combined.
	"""
	params = dict()
	plate_params = dict()
	for cp in cps:
		plate_params.setdefault(str(cp), [])
		for parameter, mu, obs_std in all_means(cp):
			params.setdefault(parameter, [])
			plate_params[str(cp)].append((parameter, mu, obs_std))
			params[parameter].append((mu, obs_std))
	for param in params:
		params[param] = np.mean(zip(*params[param]), axis=1)
	for cp in cps:
		for parameter, mu, obs_std in plate_params[str(cp)]:
			print parameter, 
			for finger in cp.fingerprints():
				if obs_std == 0 or params[parameter][1] == 0:
					new = finger[parameter] - mu + params[parameter][0]
				else:
					new = (finger[parameter] - mu + params[parameter][0]) * \
							(params[parameter][1]/obs_std)
				finger[parameter] = cp[parameter][str(finger)] = new
				
def print_cps(filename, *cps):
	"""Function takes a file object name and writes the cp data matrix to that open file
	as a tab delimited text file.
	"""
	f = open(filename, 'w')
	runs = [finger for cp in cps for finger in cp.fingerprints()]
	feats = cps[0].features()
	f.write("Features\t{}\n".format("\t".join([str(feat) for feat in feats])))
	for run in runs:
		f.write("{}\t{}\n".format(str(run), "\t".join([str(run[str(feat)]) for feat in feats])))
	f.close()
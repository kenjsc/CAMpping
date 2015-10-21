"""Module contains several tools for analyzing mass spec data including basketing compounds and
adducts, comparing isotope patterns, and merging idruns.
"""

import sys

def ppm(m1, m2):
        """Function accepts two exact masses and returns the absolute value of the ppm
        difference between the two masses.
        """
	return abs(float(m1) - m2) / m2 * 1000000
	
def iso_diff(pat1, pat2):
        """Function accepts to isotope patterns (represented as a list of (m/z, abundance) pairs)
        and returns the isotope difference product score (adapted from Pluskal, et al).
        For two isotopes to be considered the same exact mass, they have to be within 20 ppm.
        """
	p = 0
	n = 0
	diffs = []
	while n < len(pat1):
		while p < len(pat2):
			if ppm(pat1[n][0], pat2[p][0]) < 20:
				diffs.append(1.0 - abs(pat1[n][1] - pat2[p][1]))
				p += 1
				n += 1
				break
			elif pat1[n][0] < pat2[p][0]:
				diffs.append(1.0 - pat1[n][1])
				n += 1
				break
			elif pat1[n][0] > pat2[p][0]:
				diffs.append(1.0 - pat2[p][1])
				p += 1
				break
			elif n >= len(pat1) or p >= len(pat2):
				break
		if p >= len(pat2):
			while n < len(pat1):
				diffs.append(1.0 - pat1[n][1])
				n += 1
		continue
	while p < len(pat2):
		diffs.append(1.0 - pat2[p][1])
		p += 1
	dif = 1
	for d in diffs:
		dif *= d
	return dif

def add_match(id1_adds, id2_adds, dppm, drt, diso=0.5, pair=False):
	"""Function bins adducts by a ppm and retention time difference. Input is a list of adducts from
	idrun1, a list of adducts from idrun2, a ppm difference, a retention time difference, and a pair
	flag. Function returns a list of tuples where index 0 is the index of each adduct in the idrun1
	adduct list, and index 1 is the index of the adduct match in the idrun2 adduct list (value is
	None if no match is found). If the pair flag is set to true, the function instead returns tuples
	where index 0 is the adduct from idrun1 list of adducts, and index 1 is the matched adduct from
	idrun2 list of adducts (or None if there was no match). Note, the isotope list in each adduct
	needs to be sorted by volume in reverse order for function to work properly."""
	#List of the same length as id1 adducts and each value is None if no match is found in the
	# other adduct list, or the index of the match in the other adduct list
	matches = [None] * len(id1_adds)
	for id1_ind, id1_add in enumerate(id1_adds):
		for id2_ind, id2_add in enumerate(id2_adds):
			ppmdif = ppm(id1_add.isotopes[0].mz, id2_add.isotopes[0].mz)
			rtdif = abs(id1_add.isotopes[0].rt - id2_add.isotopes[0].rt)
			isodif = iso_diff(id1_add.get_pat(), id2_add.get_pat())
			#remember to put back 'and isodif>=diso'
			if ppmdif <= dppm and rtdif <= drt and isodif >= diso:
				matches[id1_ind] = id2_ind
				break
	if pair:
		#pair adducts together and return the baskets
		baskets = []
		for id1_ind, id2_ind in enumerate(matches):
			if not id2_ind == None:
				baskets.append((id1_adds[id1_ind], id2_adds[id2_ind]))
			else:
				baskets.append((id1_adds[id1_ind], None))
		return baskets
	return matches
	
def comp_match(id1_comps, id2_comps, ppm, drt, pair=False):
	"""Function bins compounds by a ppm and retention time difference. Input is a list of compounds
	from idrun1, a list of compounds from idrun2, a ppm difference, a retention time difference, and
	a pair flag. Function returns a list of tuples where index 0 is the index of each compound in
	the idrun1 compound list, and index 1 is the index of the compound match in the idrun2 compound
	list (value is None if no match is found). If the pair flag is set to true, the function instead
	returns tuples where index 0 is the compound from idrun1 list of compounds, and index 1 is the
	matched compound from idrun2 list of compounds (or None if there was no match). Note, the
	isotope list in each adduct of the compounds needs to be sorted by volume in reverse order for
	function to work properly."""
	#List of is the same length as id1 compounds and each value is None if no match is found in the
	# other compound list, or the index of the match in the other compound list
	matches = [None] * len(id1_comps)
	for id1_ind, id1_comp in enumerate(id1_comps):
		for id2_ind, id2_comp in enumerate(id2_comps):
			ppmdif = abs(id1_comp.mass - id2_comp.mass) / id1_comp.mass * 1000000
			rtdif = abs(id1_comp.rt - id2_comp.rt)
			if ppmdif <= ppm and rtdif <= drt:
				matches[id1_ind] = id2_ind
				break
	if pair:
		#pair compounds together and return the baskets
		baskets = []
		for id1_ind, id2_ind in enumerate(matches):
			if id2_ind:
				baskets.append((id1_comps[id1_ind], id2_comps[id2_ind]))
			else:
				baskets.append((id1_comps[id1_ind], None))
		return baskets
	return matches
	
def comb_runs(dppm, drt, *runs):
	"""Function takes any number of IDrun objects and combines all the adducts into one run by
	merging those with less than the specified retention time and ppm difference. Note: this 
	function payes no attention to compounds, only adducts are merged.
	"""
	if len(runs) == 1:
		#Recursive function, break recursion when there are no more blanks to combine
		return runs[0]
	if not runs[0] and not runs[1]:
		#If both blanks are None type, skip to next level of recursion
		return comb_runs(dppm, drt, None, *list(runs)[2:])
	#If either blank is None type, return it (or enter next level of recursion without it)
	elif not runs[0]:
		return comb_runs(dppm, drt, runs[1], *list(runs)[2:])
	elif not runs[1]:
		return comb_runs(dppm, drt, runs[0], *list(runs)[2:])
		
	#Get bins of adducts
	del_list = add_match(runs[0].adducts, runs[1].adducts, 10, 0.4)
	for add, rm in zip(runs[0].adducts, del_list):
		#Add adduct to the second blank from the first if it wasn't a match between the two.
		if rm:
			continue
		else:
			add.compound = None
			add.idrun = runs[1]
	#Continue recursion with all but the first blank (because it was merged into the second)
	return comb_runs(dppm, drt, *runs[1:])
	
			
def rmblank(idrun, blank, drt=0.4, dppm=10):
	"""Function takes two IDrun objects: first, one from an idrun, the second from a blank run and
	subtracts peaks from the idrun that are in the blank within specified tolerance (default is 0.4
	minutes and 10 ppm).
	"""
	#Catch function if idrun or blank are None type objects
	if not idrun:
		return None
	if not blank:
		sys.stderr.write("ERROR: Couldn't Find Blank.\n")
		exit()
	del_list = add_match(idrun.adducts, blank.adducts, dppm, drt)
	
	
	for add, rm in zip(idrun.adducts, del_list):
		#Remove the adducts that are in the blanks
		if not rm:
			continue
		else:
			add.compound = None
			add.idrun = None


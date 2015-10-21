"""Program validates peaks reported by MassHunter between 2GHz and 4GHz mode, matches up 
positive and negative data when provided, and cleans up saturated peaks (which MassHunter
will pick several times). All functions use IDrun objects for comparison and editing.
"""

#Algorithm for validating compounds/peaks in prefractions.

import sys
import btools

def gigmatch(id1, id2):
	"""Function that takes a particular mode (positive or negative) of mass spec data (2 GHz
	and 4 Ghz) and validates the data by comparing 2 GHz to 4 GHz. BQ scores are given based on
	the reproducibility of the data between 2 and 4 GHz, and the bq values are generated based
	on which mass spec data (rt, m/z, intensity, abundance, etc.) is stored for the sets of data
	determined to be of the same adduct. 4 GHz data is favored unless it is saturated, when 2
	GHz becomes favored. BQ scores are generated with binary digits as follows: 
		X				X				X				X				X
		peak			peak			peak			peak			2GHz to 4GHz
		identified		saturated		identified		saturated		match within
		in 2 GHz		in 2 GHz		in 4 GHz		in 4 GHz		10 ppm (0) else
		(1), else (0)	(1), else (0)	(1), else (0)	(1), else (0)	within 20 ppm (0)
	."""
	
	#four and two are supposed to become an IDrun object that is 4GHz and 2GHz, respectively
	four = None
	two = None

	for idrun in [id1, id2]:
		if len(idrun.adducts) == 0:
			continue
		if all(add.hertz == 4 for add in idrun.adducts):
			four = idrun
		elif all(add.hertz == 2 for add in idrun.adducts):
			two = idrun
		else:
			sys.stderr.write("ERROR: Couldn't determine 2GHz from 4GHz data.\n")
			sys.exit()

	if four == None and two == None:
		return id1
	elif four == None:
		if two == id1:
			four = id2
		elif two == id2:
			four = id1
	elif two == None:
		if four == id1:
			two = id2
		elif four == id2:
			two = id1
	
#	print "BIGGENING, BIGEFSKLDJGKS:DJKG:LSJDKF:LJSDK:LFJSKLD:FJKSL:D"
	bins = btools.add_match(sorted(four.adducts, key=lambda ad:ad.isotopes[0].mz), two.adducts, 10, 0.4, pair=True)
	for add4, add2 in bins:
#		print "BINNING: ", add4.isotopes[0].mz, add4.isotopes[0].rt, add2
		if add4.isotopes[0].sat and add2:
			if add2.isotopes[0].sat:
				#If the 4 and 2 are saturated, pull the 2, and replace the 4 with the 2
				#Also assign the bq of "11110"
				add2.bq = "11110"
				add2.compound = None
				add2.idrun = None
				add2.compound = add4.compound
				add2.idrun = add4.idrun
				add4.compound = None
				add4.idrun = None
			else:
				#If the 4 is saturaged, but 2 isn't, pull the 2, and replace the 4 with the 2
				#Also assign the bq of "10110"
				add2.bq = "10110"
				add2.compound = None
				add2.idrun = None
				add2.compound = add4.compound
				add2.idrun = add4.idrun
				add4.compound = None
				add4.idrun = None
		elif not add4.isotopes[0].sat:
			if add2:
				if add2.isotopes[0].sat:
					#If the 4 is unsaturated but the 2 is, delete both of the adducts.
					#This is not a possible scenario, since 4 is more sensitive than 2.
					add2.idrun = None
					add4.idrun = None
					add2.compound = None
					add4.compound = None
				else:
					#If the 4 and the 2 are unsaturated this is ideal. Keep the 4,
					# pull the 2 so it isn't in the 20ppm binning
					#Also assign the bq of "10100"
					add2.idrun = None
					add2.compound = None
					add4.bq = "10100"
			else:
				#If the 4 is unsaturated, but no 2 match was found, just store the 4
				#Also assign the bq of "00100"
				add4.bq = "00100"
#		if add2:
#			print add4.isotopes[0].mz, add4.bq, add2.isotopes[0].mz, add2.bq
#		else:
#			print add4.isotopes[0].mz, add4.bq
				
	#Repeat for 20 ppm binning window
	bins = btools.add_match([add for add in four.adducts if not add.bq], two.adducts, 20, 0.4,
			pair=True)
	for add4, add2 in bins:
		if add4.isotopes[0].sat:
			if add2:
				if add2.isotopes[0].sat:
					#If the 4 and the 2 are saturated, pull the 2, and replace the 4 with 2
					#Also assign the and bq of "11111"
					add2.bq = "11111"
					add2.compound = None
					add2.idrun = None
					add2.compound = add4.compound
					add2.idrun = add4.idrun
					add4.compound = None
					add4.idrun = None
				else:
					#If the 4 is saturated but not the 2, pull the 2, and replace 4 with 2
					#Also assign the bq of "10111"
					add2.bq = "10111"
					add2.compound = None
					add2.idrun = None
					add2.compound = add4.compound
					add2.idrun = add4.idrun
					add4.compound = None
					add4.idrun = None
			else:
				#If the 4 is saturated but no 2 match is found, store the 4
				#Also assign the bq of "00111"
				add4.bq = "00111"
	#All unsaturated 4 GHz peaks should have allready been removed, so only deal with remaining
	# 2GHz peaks
	while len(two.adducts):
		add2 = two.adducts.pop(0)
		if add2.isotopes[0].sat:
			add2.idrun = None
			add2.compound = None
			continue
		else:
			add2.compound.adducts = [add for add in add2.compound.adducts if not\
					add.isotopes[0].sat]
			add2.compound.idrun = None
			add2.compound.idrun = four
			for add in add2.compound.adducts:
				add.idrun = None
				add.idrun = four
				add.bq = "10001"	
	return four

def addmatch(pos, neg):
	"""Takes positive and negative mode mass data and combines adducts into compounds based on
	the molecular mass that is given to the positive mode and negative mode compounds. If the
	ppm between the masses is less than 10 and the retention time less than 0.4 min (24
	seconds)they are called the same compound."""
	if pos and neg:
		#Get compounds that bin from the negative to the positive
		bins = btools.comp_match(pos.compounds, neg.compounds, 10, 0.4, pair=True)
	else:
		return pos if pos else neg
	for pcomp, ncomp in bins:
		if ncomp:
			while len(ncomp.adducts):
				#Transfer the adducts from the negative compound to the positive and reset mode
				ncomp.adducts[0].idrun = pcomp.idrun
				ncomp.adducts[0].compound = pcomp
				pcomp.mode = "B"
	#For compounds left in negative mode (weren't matched to positive), transfer them over to
	#the positive IDrun
	for ncomp in neg.compounds:
		if not len(ncomp.adducts):
			continue
		ncomp.idrun = pos
		for add in ncomp.adducts:
			add.idrun = pos
	return pos

def rmsat(idrun):
	"""Function accepts an IDrun object and edits the IDrun in place by finding saturated 
	adducts and removing any other adducts within 1 dalton and 0.4 minutes as they are most
	likely improperly picked by MassHunter.
	"""
	#Sort adducts by adbundance so saturated peaks are first
	idrun.adducts = sorted(idrun.adducts, key=lambda add: add.isotopes[0].abundance,\
			reverse=True)
	n = 0
	while n < len(idrun.adducts):
		#Only remove ghost peaks if the largest is saturated
		if idrun.adducts[n].isotopes[0].sat == 1:
			p = n + 1
			while p < len(idrun.adducts):
				mzdif = abs(idrun.adducts[p].isotopes[0].mz - idrun.adducts[n].isotopes[0].mz)
				rtdif = abs(idrun.adducts[p].isotopes[0].rt - idrun.adducts[n].isotopes[0].rt)
				if mzdif < 1 and rtdif < 0.4:
					#Remove all peaks within a retention time and m/z tolerance of the largest
					#peak that was saturated
					ad = idrun.adducts.pop(p)
					ad.idrun = None
					ad.compound = None
					continue
				else:
					p += 1
			n += 1
		else:
			break

def abund_test(idrun, cutoff=1000, rtcut=0.35):
	"""Function takes an idrun object and edits the IDrun in place by updating the abundance of each
	compound (in case it has gained or lost adducts during validationg) and removing it if the
	abundance is below the specified cutoff (Default is 1000).
	"""
	n = 0
	while n < len(idrun.compounds):
		comp = idrun.compounds[n]
		#Reset the abundance of the compounds
		comp.get_abund()
		if comp.abundance < cutoff or comp.rt < rtcut:
			#Remove adducts and compound
			for add in comp.adducts:
				add.idrun = None
			comp.idrun = None
			continue
		n += 1

def validate(cutoff, rtcut, *idruns):
	"""Function takes either 2 or 4 IDrun objects (2 of 2GHz and 4GHz of either positive or negative
	mode; 4 of 2GHz, 4GHz, positive, and negative modes) and merges 2GHz and 4GHz files according to
	gigmatch rules and matches positive and negative features if applicable.
	"""
	if len(idruns) != 2 and len(idruns) != 4:
		sys.stderr.write("""Incorrect file input. Software takes either 2 or 4 files. There should
be both 2Ghz and 4Ghz files for either or both positive and negative data.""")
		sys.exit(1)
	for idrun in idruns:
		rmsat(idrun)
	#pos and neg are lists of positive and negative IDrun objects
	for idrun in idruns:
		for ad in idrun.adducts:
			if ad.isotopes[0].mz == 647.1189:
				print idrun.cef_file, ad.isotopes
	pos = []
	neg = []
	for idrun in idruns:	# added 20140224
		if not len([compound.mode for compound in idrun.compounds]):
			if "POS" in idrun.cef_file.upper():
				pos.append(idrun)
			elif "NEG" in idrun.cef_file.upper():
				neg.append(idrun)   #####
		elif all([compound.mode == "+" for compound in idrun.compounds]):
			pos.append(idrun)
		elif all([compound.mode == "-" for compound in idrun.compounds]):
			neg.append(idrun)
	if len(pos):
		pos = gigmatch(*pos)
	else:
		pos = None
	if len(neg):
		neg = gigmatch(*neg)
	else:
		neg = None
	if pos and neg:
		pos = addmatch(pos, neg)
	abund_test(pos, cutoff=cutoff, rtcut=rtcut)
	return pos

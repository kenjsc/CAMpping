"""Module for importing data from *.cef files from MassHunter. Main function is cefimport which
takes a *.cef filename as input and returns an IDrun object with it's compounds, adducts, and
isotopes.


"""
#Written by Emerson Glassey
#around 12-16-2012

import sys
import dbtables
import re
import os


def pfind(line, *preceeding):
	"""Function takes a string which is a line from a *.cef file and any number of strings that are
	representative of data in the *.cef files and returns the value of the data (in the same order).
	Example is line: <rt=5.0, m/z=300.0987>
	       with preceeding: rt, m/z
	       will return: (5.0, 300.0987)
	"""
	values = []
	for term in preceeding:
		#find the term value, which there should only be one (index 0). Pull the value out
		#from the search term.
		t = re.search(r"{}=\"[\d\.A-Za-z\-\+]+\"".format(term), line)
		values.append(t.group().split("=")[-1][1:-1])
	return values
			
def add_adducts(openfile, idrun, compound, mode, hertz):
	"""Function takes an open *.cef file object, an idrun object, a compound object, and the hertz
	mode. Parses the file until the next compound is reached, populating the adducts and isotope
	tables with relevant data from the *.cef file.
	"""
	#Read the next line in the file, continuing until '<p' (adduct/isotope line) is reached
	line = openfile.readline().strip()
	while "<p" not in line:
		line = openfile.readline().strip()
	#Create a new adduct object and set the parent ion. new lines will be added to this adduct until
	#a new ion species is reached, then a new adduct will be created and parent will be reset.
	new_add = dbtables.Adduct(mode, hertz=hertz)
	#ion/adduct name, ex: M+H
	adduct = pfind(line, "s")[0]
	parent = adduct
	while True:
		if "/MSPeaks" in line:
			#break as soon as the end of the isotope/adduct info section is passed.
			break
		if "<p" not in line:
			#Skip lines that don't have any information (aren't '<p')
			line = openfile.readline().strip()
			continue
		#Get the isotope data and convert it to floats
		rt, mz, abundance, volume = map(float, pfind(line, "rt", "x", "y", "v"))
		#also get the charge and adduct information
		charge, adduct = pfind(line, "z", "s")
		#'sat="True"' is only in there if the peak is saturated, so store saturation tag
		sat = 1 if "sat" in line else 0
		if adduct[:len(parent)] == parent:
			#If the read isotope is the same ion as the parent, then create a new isotope
			# object and append it to the parent adduct
			new_add.isotopes.append(dbtables.Isotope(rt, mz, abundance, volume, sat, charge,
					adduct))
		else:
			#If the read isotope is a new ion, sort the set of isotopes from the last adduct
			#by m/z
			new_add.isotopes = sorted(new_add.isotopes, key=lambda isotope: isotope.mz)
			if new_add.calc_pat():
				#append the last adduct to the compound object adduct list and idrun object list
				compound.adducts.append(new_add)
				idrun.adducts.append(new_add)
			#create a new adduct and add this new isotope to it, also reset the new parent ion
			new_add = dbtables.Adduct(mode, hertz=hertz)
			new_add.isotopes.append(dbtables.Isotope(rt, mz, abundance, volume, sat, charge,
					adduct))
			parent = adduct
		#continue to the next line in the file
		line = openfile.readline().strip()
	#as a last step, make sure to sort the isotopes by m/z and append the last adduct
	new_add.isotopes = sorted(new_add.isotopes, key=lambda isotope: isotope.mz)
	if new_add.calc_pat():
		compound.adducts.append(new_add)
		idrun.adducts.append(new_add)

def add_compound(openfile, idrun, hertz):
	"""Function takes an open *.cef file object, an idrun object and the hertz mode. Parses the
	file, populating the compound table, and linking to adducts/isotopes with data from the *.cef
	file.
	"""
	#Read the next line in the file, continuing until the beginning of the first compound
	#at the "Location" flag.
	line = openfile.readline().strip()
	#Set rt and mode to None to keep track of when their values get assigned
	rt = None
	mode = None
	while "</CEF>" not in line:
		if "Location" not in line and not rt:
			line = openfile.readline().strip()
			continue
		if "Location" in line:
			#Get the compound information and convert to floats
			mass, rt, abund, volume =  map(float, pfind(line, "m", "rt", "y", "v"))
		#Skip ahead to the mode information
		if "MSDetails" not in line and not mode:
			line = openfile.readline().strip()
			continue
		if "MSDetails" in line:
			#Get the mode (pos/neg)
			mode = pfind(line, "p")[0]
		#Create the new compound object and add its object
		new_comp = dbtables.Compound(rt, mass, mode, abund, volume)
		add_adducts(openfile, idrun, new_comp, mode, hertz)
		#After adding the adducts, append the compound to the idrun
		idrun.compounds.append(new_comp)
		#Reset the rt and mode to none to find the next compound entry
		rt = None
		mode = None
		line = openfile.readline().strip()

def read_cef(cef):
	"""Function takes a *.cef filename and imports the data from the *.cef file into an IDrun object
	with compounds/adducts/isotopes. Returns the IDrun object.
	
	Notes about format of file names:
		The filename needs to include the extract number, prefraction, mode (pos/neg), mode
			(4GHz/2GHz), location code, and minute (if peak library).
		The format of this data in the filename is slightly specific. The easiest way to do things
			is follow the template:
			
				RLUS-1478C-45_POS_4G.cef
				
			The important aspects that MUST be followed:
				
				RLUS-1478C		-		4 letter country code *dash* extract fraction
				
				or 
				
				RLUS-1487C-45	-		*dash* minute in peak library
				
		If these aspects are followed, the pos/neg and 4G/2G can be anywhere in the filename.
		It also doesn't matter if they are separated by underscores, as long as the important
			motif detailed above can be picked out. But in order to avoid errors (in case my error
			checking hasn't been perfect), it's easiest to simply follow the template given.
		Names are not case sensitive as all the data is converted to uppercase, but the case does
			need to match the actual filename
	"""
	#Convert to upper case and add .cef if it's not in the filename
	cef_file = os.path.splitext(os.path.basename(cef))[0]

	#Get Hertz mode, which is 2 or 4 (or None, if the filename doesn't have it)
	
	sample_name, mode, hertz = cef_file.rsplit("_", 2)
	if "2G" in hertz:
		hertz = 2
	elif "4G" in hertz:
		hertz = 4
	else:
		sys.stderr.write("Hertz mode not supplied!\n")
		sys.exit(1)

	#Open new IDrun
	run = dbtables.IDrun(cef)
	with open(cef) as f:
		add_compound(f, run, hertz)
	return run

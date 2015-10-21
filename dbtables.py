"""Module includes classes and functions useful in the PrinCIS platform. Class objects are
mapped through sqlalchemy to a database that includes UHPLC-TOF-HRMS spectral data. Functions
are particularly useful for binning adducts and compounds and opening databases.
"""

#Emerson Glassey
#12-16-2012
#This is a definition of all of the classes for PrinCIS

import sqlalchemy
import sqlalchemy.orm
import sqlalchemy.ext.declarative
import os
import re

#In order to parallelize with pp, you cannot have from ... import clauses, so I import sqlalchemy
# and then reassign different functions within it:

declarative_base = sqlalchemy.ext.declarative.declarative_base
create_engine = sqlalchemy.create_engine
ForeignKey = sqlalchemy.ForeignKey
Table = sqlalchemy.Table
Column = sqlalchemy.Column
Integer = sqlalchemy.Integer
String = sqlalchemy.String
MetaData = sqlalchemy.MetaData
Float = sqlalchemy.Float
sessionmaker = sqlalchemy.orm.sessionmaker
relationship = sqlalchemy.orm.relationship
backref = sqlalchemy.orm.backref
mapper = sqlalchemy.orm.mapper

#Open declarative base builder to build object-db mappers
Base = declarative_base()

class IDrun(Base):
	"""IDrun class stores metadata about the run:
		rid			-		run id (assigned by the database upon import)
		file		-		run name that was saved (without .cef)
		location	-		location code (RLUS)
		extract		-		bacterial extract number
		fraction	-		prefraction letter
		name		-		combination of extract with fraction and minute (e.g. 1478D_47)
		minute		-		minute number (for peak libraries)
		
Notes about format of file names:
	The filename needs to include the extract number, prefraction, mode (pos/neg), mode
		(4GHz/2GHz), location code, and minute (if peak library).
	The format of this data in the filename is slightly specific. The easiest way to do
		things is follow the template:
		
			SampleName_POS_4G.cef
			
		The Pos/Neg and 4G/2G must be at the end in that order, and separated by '_' as
			shown above. The SampleName must be formated as follows for proper parsing:
			
			RLUS-1478C	-	4 letter country code *dash* extract fraction

			
			Import will still function without this, but the prefraction, extract,
				and fraction columns in the database will not be properly
				populated. If you are also analyzing with yeast and/or CP data, the
				SampleName must match the name in the CP or yeast files perfectly.
				
	Names are not case sensitive as all the data is converted to uppercase.
	"""
	__tablename__ = 'idruns'
	
	rid = Column(Integer, primary_key=True)
	
	#name of file run was saved as without .cef
	cef_file = Column(String)
	#Location code in 4 letter format (RLUS)
	location = Column(String)
	#4 integer extract code
	extract = Column(Integer)
	#A-F prefraction code
	fraction = Column(String)
	#concatonation of extract and prefraction codes
	prefraction = Column(String)
	#SampleName
	name = Column(String)
	
	#There is a link between IDrun and Compounds (one-to-many)
	#There is a link between IDrun and Adducts (one-to-many)

	def __init__(self, filename):
	
		self.cef_file = os.path.splitext(os.path.basename(filename))[0]	
		self.name, mode, hertz = self.cef_file.rsplit("_", 2)
		
		#Get other metadata: filename, location, extract, fraction, and minute
		
		#This is metadata for peak library or regular, respectively
		meta = re.search(r"RL[A-Z]{2}\-\d{4}[A-F]", self.name)
		if meta:
			meta = meta.group()
			self.location, prefraction = meta.split("-")
			self.extract = prefraction[:-1]
			self.fraction = prefraction[-1]
			meta = re.search(r"SYP\d[A-F]", self.name)
		elif meta:
			meta = meta.group()
			self.location = None
			self.extract = 9000 + int(meta[3])
			self.fraction = meta[4]
			self.prefraction = "{}{}".format(extract, fraction)
		else:
			self.location = None
			self.extract = None
			self.fraction = None
			self.prefraction = None
		
	def __repr__(self):
		return self.name

class Compound(Base):
	"""Compound class stores data about each compound:
		cid			-		compound id (assigned by the database upon import)
		rid			-		run id (assigned to match compound to IDrun)
		rt			-		retention time of compound
		abundance	-		height of compound peak (sum of adducts)
		volume		-		volume of compound peak (sum of adducts)
		mass		-		calculated mass of compound (based on adducts m/z)
		mode		-		polarity mode of compound ('+', '-', or 'B' if pos and neg)
	"""
	__tablename__ = 'compounds'

	cid = Column(Integer, primary_key=True)
	rid = Column(Integer, ForeignKey('idruns.rid'))
	#Retention time in minutes
	rt = Column(Float)
	#Abundance in counts
	abundance = Column(Float)
	#Volume in counts
	volume = Column(Float)
	#Mass of compound in daltons
	mass = Column(Float)
	#Mode of compound data '+', '-', or 'B' for both
	mode = Column(String)

	#relationship of compounds with IDrun (many-to-one)
	idrun = relationship(IDrun, backref=backref('compounds'))
	#There is also a link between Compound and Adducts (one-to-many)
	
	def __init__(self, rt, mass, mode, abundance, volume):
		self.rt = rt
		self.mass = mass
		self.mode = mode
		self.abundance = abundance
		self.volume = volume

	def get_abund(self):
		"""Function to calculate abundance and volume of compound peak based on its adducts. This is
		in case some adducts are removed/added, the abundance and volume must be recalculated."""
		self.abundance = sum([iso.abundance for add in self.adducts for iso in add.isotopes])
		self.volume = sum([iso.volume for add in self.adducts for iso in add.isotopes])

	def __repr__(self):
		"Representation for printing Compound objects."
		return "<COMPOUND({},{};{})>".format(self.rt, self.mass, len(self.adducts))

class Adduct(Base):
	"""Adduct class stores data about each adduct:
		aid			-		adduct id (assigned by the database upon import)
		cid			-		compound id (assigned to match adduct to compound)
		rid			-		run id (assigned to match adduct to idrun)
		bq			-		bq score (4 binary digits to denote saturation in 2Ghz and 4GHz)
		hertz		-		hertz mode the adduct was found in (4GHz or 2GHz)
	"""
	__tablename__ = 'adducts'
	
	aid = Column(Integer, primary_key=True)
	cid = Column(Integer, ForeignKey('compounds.cid'))
	rid = Column(Integer, ForeignKey('idruns.rid'))
	#bq score is a string of 4 binary digits (e.g. "1010") which denotes which modes were
	#matched up and saturated (e.ge "1010" is found in 2GHz and 4GHz, but saturated in neither)
	bq = Column(String)
	#Gigahertz mode data was acquired in (4 or 2)
	hertz = Column(Integer)
	mode = Column(String)

	#relationship of adducts with compound (many-to-one)
	compound = relationship(Compound, backref=backref('adducts'))
	#relationship of adducts with IDrun (many-to-one)
	idrun = relationship(IDrun, backref=backref('adducts'))
	#There is also a link between Adduct and Isotopes (one-to-many)

	def __init__(self, mode, hertz=None, bq=None):
		self.hertz = hertz
		self.bq = bq
		self.mode = mode
		self.calc_pat()
		
	@sqlalchemy.orm.reconstructor
	def calc_pat(self):
		if len(self.isotopes):
			if self.isotopes[0].abundance == 0:
				return False
			max = self.isotopes[0].abundance
			for iso in self.isotopes:
				iso.n_abund = float(iso.abundance) / max
			return True
		else:
			return True
			
	def get_pat(self):
		return [(iso.mz, iso.n_abund) for iso in self.isotopes]
			
		
	def __repr__(self):
		"Representation for printing Adduct objects."
		return "<ADDUCT({},{};{})>".format(self.bq, self.hertz, len(self.isotopes))

class Isotope(Base):
	"""Isotope class stores data about each isotope:
		iid			-		isotope id (assigned by the database upon import)
		aid			-		compound id (assigned to match isotope to adduct)
		rt			-		retention time of isotope
		mz			-		m/z ratio of isotope
		abundance	-		height of isotope peak
		volume		-		area of isotope peak
		sat			-		saturation flag (1 if peak was saturated, 0 otherwise)
		ion			-		adduct ion string (e.g. "M+H+1" is second isotope of M+H)
		charge		-		absolute value of charge z (must determine polarity from ion)
	"""
	__tablename__ = 'isotopes'

	iid = Column(Integer, primary_key=True)
	aid = Column(Integer, ForeignKey('adducts.aid'))
	rt = Column(Float)
	mz = Column(Float)
	abundance = Column(Float)
	volume = Column(Float)
	sat = Column(Integer)
	ion = Column(String)
	charge = Column(String)

	#relationship of isotopes with adduct (many-to-one), sort by m/z ratio.
	adduct = relationship(Adduct, backref=backref('isotopes'), order_by=mz)

	def __init__(self, rt, mz, abundance, volume, sat, charge, ion):
		self.rt = rt
		self.mz = mz
		self.abundance = abundance
		self.volume = volume
		self.sat = sat
		self.charge = charge
		self.ion = ion
		#normalized abundance (normalized within the adduct), to be calculated later.
		self.n_abund = None

	@sqlalchemy.orm.reconstructor
	def norm_abundance(self):
		self.n_abund = None
	
	def __repr__(self):
		"Representation of printing Isotope objects."
		return "<ISOTOPE({},{},{},{},{},{},{})>".format(self.rt, self.mz, self.abundance,\
				self.volume, self.sat, self.charge, self.ion)

def connect(dbdir, verb=False):
	"""Function takes a sqlite database filename and an optional verbose flag (Default is False).
	Returns an open sqlalchemy database connection engine which can be used to open a session."""
	db = dbdir
	engine = create_engine('sqlite:///'+db, echo=verb)
	metadata = Base.metadata
	metadata.create_all(engine)
	return engine
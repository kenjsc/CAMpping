README

I have all of the code in the same folder for two reasons:
It makes adding the path to python easier
I use many of these modules between projects

Breakdown of Files

CAMpping Code:

analyze.py - main analysis python file. Imports many other modules. Interactive Python tool.
btools.py - misc. functions for data analysis. things like peak alignment, etc.
mzbin.py - defines an object for analyzing the MS data (baskets and prefractions that contain peaks).




cefimport.py - imports four cef files per prefraction to db. Command line tool.
cefparse.py - imported by cefimport. Parses the files and reads cef files into python

val.py - functions for performing validation between 4 modes.

dbtables.py - defines all of the tables in the sqlite database, and how those map to python objects.

dynamic_gexf.py - python tool for compressing .dot files from networks into a massive dynamic network in .gexf form.
dynamic_network.py - another attempt at the same thing. Never got it to work. Use dynamic_gexf instead.

example_cefimport.txt - these are the commands I used to generate the mass database used for all the data analysis we did on plate 4

likeplots.py - command line utility. generates line plots for the likelihood plot data from a csv or txt file.

lognormal.py - includes functions for fitting a lognormal to data.

meta2out.py - this generates the csv file with all the likelihood plot data. It converts the flat file per iteration data to a 2 dimensional matrix. the meta.txt file is written during vector compression.

networks.py - function for generating co-expression networks. Still in beta.


CP:
batch.py - tries to fix batch mode issues.

chosen_dilutions.py - analyzes tab files. Can deal with multiple of the same prefraction at different dilutions, unlike cp.py. Also includes IC50 curve fitting function.
cp.py - analyzes tab files. CANNOT include the same named prefractions. Name portion after underscore is discarded, so concentration information is removed. I did this because our extracts don't include concentration (unless we do a dilution series, which is why there is chosen_dilutions.py)


cluster.py - old module for analyzes cp data. All functionality is now included in cp.py
cp_iterations.py - I honestly don't remember why this file exists. As far as I can tell, all functionality is now in cp.py

difscore.py - analyzes a cdt file. exports a csv of a nxn heatmap with the difscore between prefractions where the prefractions meet in the nxn matrix.
kmeans.py - same thing as above. Uses kmeans algorithm instead though.
pearson.py - same thing as above. Uses pearson correlation instead though.



BioMAP:
biomap_crude.py - analyzes crude data from HiTS. Command line tool.
biomap.py - analyzes dilution series data from HiTS. Command line tool.

pyeq2 - I DID NOT WRITE THIS. This is a curve fitting package I found online.



plotting.py - generates beautiful lineplots from data. command line tool.
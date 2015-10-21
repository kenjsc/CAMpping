
#command line tool for parsing meta.txt files

import sys
import numpy as np

def main(f, g):
	tables = dict()
	iteration = 0
	iterations = []
	for line in f:
		line = line.strip()
		if line == "":
			continue
		elif line.find("Average") != -1:
			iterations.append(np.float(line.split(": ")[1]))
			for l in tables.values():
			    l.append(0.0)
			iteration += 1
		elif len(line) == 5:
			run = line
		elif len(line.split("; ")) == 4:
			data = line.split("; ")
			name = data[0].split(", ")
			tables.setdefault((run, name[0], name[1], name[2]), [0.0])
			tables[(run, name[0], name[1], name[2])][iteration] = np.float(data[-1])
		elif line.find("Run") != -1:
#			tables.setdefault(run, [0.0])
#			tables[run][iteration] = np.float(line.split(": ")[1])
                        continue
			
	for row in sorted(tables.keys(), key=lambda x: repr(x).replace("(", "").replace(")", "")):
		g.write("{}\t{}\n".format(str(row).replace("(", "").replace(")", "").replace("'", "").\
				replace(", ", "; "), "\t".join(map(str, tables[row]))))
			
	
	
	
	
	
if __name__ == "__main__":
	main(open(sys.argv[1]), open(sys.argv[2], 'w'))
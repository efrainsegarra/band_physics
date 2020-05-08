import numpy as np
import sys

if len(sys.argv) != 4:
	print( "Incorrect number of arguments. Please use:\n" )
	print( "\tpython code.py [FADC Offset] [TDC Offset] [Output File]\n" )
	print( "" )
	exit(-1)

Offsets_FADC = {};
Offsets_TDC = {};
with open(sys.argv[1],"r") as f, open(sys.argv[2],"r") as g:
	for line in f:
		parse = line.strip().split(" ")
		sector = parse[0]
		layer = parse[1]
		component = parse[2]
		bkgrd = float(parse[3])
		norm = float(parse[4])
		mean = float(parse[5])
		sig = float(parse[6])
		integral = float(parse[7])
		flag = int(parse[8])
		
		Offsets_FADC[sector+" "+layer+" "+component] = mean
	for line in g:
		parse = line.strip().split(" ")
		sector = parse[0]
		layer = parse[1]
		component = parse[2]
		bkgrd = float(parse[3])
		norm = float(parse[4])
		mean = float(parse[5])
		sig = float(parse[6])
		integral = float(parse[7])
		flag = int(parse[8])
		
		Offsets_TDC[sector+" "+layer+" "+component] = mean

# Make veto bar offsets to 0


with open(sys.argv[3],"w") as out:
	for key in Offsets_FADC:
		out.write(key+" "+str(Offsets_TDC[key])+" "+str(Offsets_FADC[key])+"\n")

	out.write("1 6 1 0 0\n")
	out.write("1 6 2 0 0\n")
	out.write("1 6 3 0 0\n")

	out.write("2 6 1 0 0\n")
	out.write("2 6 2 0 0\n")
	out.write("2 6 3 0 0\n")
	out.write("2 6 4 0 0\n")
	out.write("2 6 5 0 0\n")
	out.write("2 6 6 0 0\n")
	out.write("2 6 7 0 0\n")

	out.write("3 6 1 0 0\n")
	out.write("3 6 2 0 0\n")
	out.write("3 6 3 0 0\n")
	out.write("3 6 4 0 0\n")
	out.write("3 6 5 0 0\n")
	out.write("3 6 6 0 0\n")

	out.write("4 6 1 0 0\n")
	out.write("4 6 2 0 0\n")
	out.write("4 6 3 0 0\n")
	out.write("4 6 4 0 0\n")
	out.write("4 6 5 0 0\n")
	out.write("4 6 6 0 0\n")

	out.write("5 6 1 0 0\n")
	out.write("5 6 2 0 0\n")

	out.close()

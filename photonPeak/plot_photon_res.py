from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


ids 	= []
res 	= []
shifts 	= []
bkg_lvl	= []
ph_peak = []

with open("build/photon_params.txt","rb") as f:
	for line in f:
		parse = line.strip().split(" ")
		
		if( parse[0] == 3 or parse[0] == 4 ): continue
		# line structure: sector, layer, component, bkg_lvl, const, mean, sigma, numEvents, bool		
		ids.append( int(parse[0])*10 + int(parse[1])*100 + int(parse[2]) )
		bkg_lvl.append( float(parse[3]) )
		ph_peak.append( float(parse[4]) )				
		shifts.append( float(parse[5]) )
		res.append( float(parse[6]) )


plt.scatter( ids, res )
plt.ylim([0,0.7])
plt.xlim([50,600])
plt.axvline(x=100,linestyle='--',color='black',linewidth=2)
plt.show()

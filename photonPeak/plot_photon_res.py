from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


ids 	= []
res 	= []
ids_s 	= []
res_s 	= []
shifts 	= []
bkg_lvl	= []
ph_peak = []

with open("../build/photonPeak/global_offset_fadc_2ndIter-10082019.txt","rb") as f:
	for line in f:
		parse = line.strip().split(" ")
		
		#if( parse[0] == 3 or parse[0] == 4 ): continue
		# line structure: sector, layer, component, bkg_lvl, const, mean, sigma, numEvents, bool		
		ID = int(parse[0])*10 + int(parse[1])*100 + int(parse[2])
		if int(parse[0]) == 3 or int(parse[0])== 4:
			ids_s.append( int(parse[0])*10 + int(parse[1])*100 + int(parse[2]) )
			res_s.append( float(parse[6]) )
		else:
			ids.append( int(parse[0])*10 + int(parse[1])*100 + int(parse[2]) )
			bkg_lvl.append( float(parse[3]) )
			ph_peak.append( float(parse[4]) )				
			shifts.append( float(parse[5]) )
			res.append( float(parse[6]) )
			


plt.scatter( ids, res ,label='Long bars',color='blue')
plt.scatter( ids_s, res_s,label='Short bars',color='red')
plt.ylim([0,0.7])
plt.xlim([50,600])
plt.ylabel("ToF Resolution [ns]",fontsize=16)
plt.xlabel('ID [a.u.]',fontsize=16)
plt.xlim([50,600])
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.legend(numpoints=1,loc=2)
plt.grid(True)
plt.tight_layout()
plt.savefig("tof_resolution.pdf")
plt.show()

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


SpB = []; B = []; MeVCut = [];
with open("spb.txt") as f:
	for line in f:
		parse = line.strip().split(" ")
		MeVCut.append( float( parse[0]) )
		SpB.append( float(parse[1]) )
		B.append( float(parse[2] ) )

SpB = np.asarray(SpB,dtype=float)
B   = np.asarray(B  ,dtype=float)
MeVCut = np.asarray(MeVCut, dtype=float)
S = SpB - B
MeVCut += 0.5

'''

plt.figure(1)
plt.errorbar(MeVCut,S,yerr=np.sqrt(SpB+B),fmt='bo',label='Signal')
plt.errorbar(MeVCut,B,yerr=np.sqrt(B),fmt='ro',label='Background')
#plt.title("")
#plt.legend(numpoints=1,loc='best')
plt.title("Signal and Background Counts\n50nA, Full Day, Neutral Candidates",fontsize=18)
plt.xlim([0,16])
plt.xlabel("MeVee Cut on Bar ADC Geometric Mean",fontsize=16)
plt.grid(True)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()
#plt.savefig("signal_background.png",bbox_inches="tight")

# f = S/B = (SpB - B) / B 
# df/dSpB = 1 / B
# df/dB = - SpB / B^2

plt.figure(2)
plt.errorbar(MeVCut,S/B,yerr=np.sqrt( (1./B)*(1./B)*SpB + (SpB)*(SpB)*(1./B)*(1./B)*(1./B)*(1./B)*B  ), fmt='bo')
plt.title("Signal-to-Background Ratio\n50nA, Full Day, Neutral Candidates",fontsize=18)
plt.xlim([0,16])
plt.xlabel("MeVee Cut on Bar ADC Geometric Mean",fontsize=16)
plt.grid(True)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()
#plt.title("")
#plt.savefig("signal2background.png",bbox_inches="tight")

plt.figure(3)
plt.errorbar(MeVCut,S/np.sqrt(SpB),fmt='bo',label='50nA')
plt.xlabel("MeVee Cut on Bar ADC Geometric Mean")
plt.legend(numpoints=1,loc='best')
plt.grid(True)

plt.figure(4)
plt.errorbar(MeVCut, S / np.sqrt( SpB + B ) ,fmt='bo',label="50nA")
plt.xlabel("MeVee Cut on Bar ADC Geometric Mean")
plt.legend(numpoints=1,loc='best')
plt.grid(True)

'''
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis


ax1.set_title("BAND Neutral Candidates, 50nA, Full Day, DIS",fontsize=18)
ax1.set_ylabel('Counts',fontsize=16) 
ax1.tick_params(labelsize=15)
ax1.set_xlabel('MeVee Cut on Bar ADC Geometric Mean',fontsize=16)
ax1.set_xlim([0,17])
ax1.errorbar(MeVCut,S,yerr=np.sqrt(SpB+B),fmt='bo',label='Signal',linestyle='--')
ax1.errorbar(MeVCut,B,yerr=np.sqrt(B),fmt='ro',label='Background',linestyle='--')
ax1.set_yscale('log')
ax1.tick_params(axis='y')


ax2.tick_params(labelsize=15)
ax2.set_ylabel('Signal-to-Background Ratio',fontsize=16,rotation=270)
ax2.errorbar(MeVCut,S/B,yerr=np.sqrt( (1./B)*(1./B)*SpB + (SpB)*(SpB)*(1./B)*(1./B)*(1./B)*(1./B)*B  ), fmt='go',linestyle='--')
ax2.tick_params(axis='y')

ax2.text(0.6,2.75,"Background",fontsize=18,color='red');
ax2.text(5.6,3.7,"Signal",fontsize=18,color='blue');
ax2.text(12,4.5,"S/B",fontsize=18,color='green');
plt.grid(True)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
fig.savefig('band-spb.pdf',bbox_inches="tight");

plt.show()

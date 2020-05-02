import matplotlib.pyplot as plt
import numpy as np


# Read in acceptance-corrected lowx/highx counts
acc_corr_lox = []
acc_corr_err_lox = []
acc_corr_pn_lox = []
acc_corr_hix = []
acc_corr_err_hix = []
acc_corr_pn_hix = []
with open("../build/crossSection/accCorr-lox.txt") as f, open("../build/crossSection/accCorr-hix.txt") as g:
	for line in f:
		parse = line.strip().split(" ")
		acc_corr_pn_lox.append( parse[0] )
		acc_corr_lox.append( parse[1] )
		acc_corr_err_lox.append( parse[2] )
	for line in g:
		parse = line.strip().split(" ")
		acc_corr_pn_hix.append( parse[0] )
		acc_corr_hix.append( parse[1] )
		acc_corr_err_hix.append( parse[2] )

acc_corr_lox = np.asarray(acc_corr_lox,dtype=float)
acc_corr_hix = np.asarray(acc_corr_hix,dtype=float)
acc_corr_pn_lox = np.asarray(acc_corr_pn_lox,dtype=float)
acc_corr_pn_hix = np.asarray(acc_corr_pn_hix,dtype=float)
acc_corr_err_lox = np.asarray(acc_corr_err_lox,dtype=float)
acc_corr_err_hix = np.asarray(acc_corr_err_hix,dtype=float)

# Read in true lowx/highx counts
tru_lox = []
tru_hix = []
tru_pn_lox = []
tru_pn_hix = []
tru_err_lox = []
tru_err_hix = []

with open("../build/crossSection/true-lox.txt") as f, open("../build/crossSection/true-hix.txt") as g:
	for line in f:
		parse = line.strip().split(" ")
		tru_pn_lox.append( parse[0] )
		tru_lox.append( parse[1] )
		tru_err_lox.append( parse[2])
	for line in g:
		parse = line.strip().split(" ")
		tru_pn_hix.append( parse[0] )
		tru_hix.append( parse[1] )
		tru_err_hix.append( parse[2] )

tru_lox = np.asarray(tru_lox,dtype=float)
tru_hix = np.asarray(tru_hix,dtype=float)
tru_pn_lox = np.asarray(tru_pn_lox,dtype=float)
tru_pn_hix = np.asarray(tru_pn_hix,dtype=float)
tru_err_lox = np.asarray(tru_err_lox,dtype=float)
tru_err_hix = np.asarray(tru_err_hix,dtype=float)

# Create the ratio
ratio_lox = acc_corr_lox / tru_lox
ratio_hix = acc_corr_hix / tru_hix
# Create error in ratio
ratio_err_lox = np.sqrt( np.power(1./tru_lox,2)*np.power(acc_corr_err_lox,2) + np.power(acc_corr_lox,2) / np.power(tru_lox,4) * np.power(tru_err_lox,2) )
ratio_err_hix = np.sqrt( np.power(1./tru_hix,2)*np.power(acc_corr_err_hix,2) + np.power(acc_corr_hix,2) / np.power(tru_hix,4) * np.power(tru_err_hix,2) )

plt.figure(3)
plt.title('Ratio of Det. Acc. Corrected to Gen. Counts')
ratio_err_lox = np.sqrt( np.power(np.sqrt(acc_corr_lox),2)/np.power(tru_lox,2) + np.power(np.sqrt(tru_lox),2)*np.power(acc_corr_lox,2)/np.power(tru_lox,4)    )
ratio_err_hix = np.sqrt( np.power(np.sqrt(acc_corr_hix),2)/np.power(tru_hix,2) + np.power(np.sqrt(tru_hix),2)*np.power(acc_corr_hix,2)/np.power(tru_hix,4)    )
plt.errorbar( tru_pn_lox , ratio_lox , yerr = ratio_err_lox , color = 'purple', label="Low x'", fmt ='o')
plt.errorbar( tru_pn_hix , ratio_hix , yerr = ratio_err_hix , color = 'cornflowerblue', label="High x'", fmt ='o')
plt.legend(numpoints=1,loc='best')
plt.ylim([0.5,1.5])
plt.xlim([0.2,0.60])
plt.xlabel(r"Neutron Momentum $p_n$ MeV/c")
plt.ylabel("Det. Acc. Corr. / Gen.")

plt.grid(True)
plt.axhline(y=1,linestyle='--',linewidth=3,color='black')

plt.savefig("acceptance_corrected.pdf",bbox_inches="tight")

plt.show()

'''
plt.figure(1)
plt.title("Low X': Gen. vs Det. Acc. Corrected Cross Section")
plt.errorbar( pn[:-1], tru_lox[:-1], marker='o',linestyle='--',color='green',label="Gen.",s=100)
plt.errorbar( pn[:-1], acc_corr_lox[:-1], marker='o',linestyle='--',color='red',label="Det. Acc. Corr.",s=100)
plt.yscale('log')
plt.legend(numpoints=1,loc='best')
plt.ylim([1E-6,1E-2])
plt.xlim([200,600])
plt.xlabel(r"Neutron Momentum $p_n$ MeV/c")
plt.ylabel(r"Cross Section (x',Q2,$p_n$)")
plt.grid(True)

plt.savefig('lowx.pdf',bbox_inches="tight")

plt.figure(2)
plt.title("High X': Gen. vs Det. Acc. Corrected Cross Section")
plt.scatter( pn, tru_hix, marker='*',linestyle='--',color='green',s=100,label="Gen.")
plt.scatter( pn, acc_corr_hix, marker='*',linestyle='--',color='red',s=100,label="Det. Acc. Corr.")
plt.yscale('log')
plt.legend(numpoints=1,loc='best')
plt.ylim([1E-6,1E-2])
plt.xlim([200,600])
plt.xlabel(r"Neutron Momentum $p_n$ MeV/c")
plt.ylabel(r"Cross Section (x',Q2,$p_n$)")
plt.grid(True)

plt.savefig('highx.pdf',bbox_inches="tight")

plt.figure(3)
plt.title('Ratio of Det. Acc. Corrected to Gen. Counts')
acc_corr_lox=np.asarray(acc_corr_lox[:-1])
acc_corr_hix=np.asarray(acc_corr_hix)
tru_lox=np.asarray(tru_lox[:-1])
tru_hix=np.asarray(tru_hix)
plt.scatter( pn[:-1], acc_corr_lox / tru_lox , marker ='o', s=100, color='blue',label="Low x'")
plt.scatter( pn, acc_corr_hix / tru_hix ,marker = "*", s=100, color='purple',label="High x'")
plt.legend(numpoints=1,loc='best')
plt.ylim([0,1.5])
plt.xlim([200,600])
plt.xlabel(r"Neutron Momentum $p_n$ MeV/c")
plt.ylabel("Det. Acc. Corr. / Gen.")
plt.grid(True)

plt.savefig('ratio.pdf',bbox_inches="tight")

plt.show()'''

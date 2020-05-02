
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


acc_corr_hix = acc_corr_hix*0.039915337
acc_corr_lox = acc_corr_lox*0.039915337

acc_corr_err_hix = np.sqrt(acc_corr_hix)
acc_corr_err_lox = np.sqrt(acc_corr_lox)


plt.figure(3)
plt.errorbar( acc_corr_pn_hix , acc_corr_hix ,  color = 'cornflowerblue', label="|x'-0.505|<0.005", fmt ='o')
plt.errorbar( acc_corr_pn_lox , acc_corr_lox ,  color = 'purple', label="|x'-0.3|<0.005", fmt ='o')
plt.legend(numpoints=1,loc='best')
plt.ylim([0,120])
plt.title("W'>2, Q2>2, CosTheta_qn < -0.9, -30 < Phi_qn < 30 ")
plt.xlim([0.2,0.60])
plt.xlabel(r"Neutron Momentum $p_n$ MeV/c")
plt.ylabel("Counts")

plt.grid(True)
plt.tight_layout()
plt.savefig("deeps-plot.pdf",bbox_inches="tight")

plt.show()

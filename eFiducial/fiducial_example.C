#include "ElectronFiducial.cc"
#include "TRandom3.h"


void fiducial_example() {

	TRandom3* myRand = new TRandom3();
	int seed = time(0) - getpid();
	myRand->SetSeed(seed);

	TH2F* acc_hist = new TH2F("acc_hist", "", 420, -220, 200, 160, 0, 40); 

	int nEvents = 100000;

	double theta, phi, pe, acc;

	ElectronFiducial* myFiducial = new ElectronFiducial("upper_momentum_fit.dat", "lower_momentum_fit.dat");

	for(int i = 0; i < nEvents; i++) {
		theta = myRand->Uniform(0., 40.);
		phi = myRand->Uniform(-180., 180.);
		pe = myRand->Uniform(2.,10.);
		acc = myFiducial->GetElectronAcceptance(theta, phi, pe);

		acc_hist->Fill(phi, theta, acc);
	}

	gStyle->SetOptStat(0);

	acc_hist->Draw("COLZ");

}

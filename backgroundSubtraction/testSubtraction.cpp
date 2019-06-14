#include <cmath>
#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TRandom3.h"

using namespace std;

int main( int argc, char ** argv){

	if( argc != 3 ){
		cerr << "Wrong number of arguments. Please instead use:\n\t./test [inputFile] [outputFile]\n";
		return -1;
	}
	
	TFile * inFile = new TFile(argv[1]);
	TFile * outFile = new TFile(argv[2],"RECREATE");

	TTree * inTree = (TTree *) inFile->Get("skim");

	const int nS = 253; // this is integrated signal in our in-time region of interest
	const int nB = 230; // this is integrated background in some off-time region
	const double lowerBedge = -8.;
	const double upperBedge = 2;
	
	// Let's generate a signal ToF spectrum just drawing from a trianular PDF:
	TF1 * tri = new TF1("tri",	"(x<6)*0 						+\
					(x>=6 && x<16)*(x-6.)*(1./5)*(1./10)			+\
					(x>=16)*0",		0,40				);
	TH1D * ToF_bkgrd = new TH1D("ToF_bkgrd","ToF_bkgrd",2000,-100,100);
	TH1D * ToF_signal = new TH1D("ToF_signal","ToF_signal",2000,-100,100);
	TH1D * ToF_spb = new TH1D("ToF_spb","ToF_spb",2000,-100,100);

	TH1D * mom_spb = new TH1D("mom_spb","mom_spb",100,0,1);
	TH1D * mom_b = new TH1D("mom_b","mom_b",100,0,1);
	TH1D * mom_s = new TH1D("mom_s","mom_s",100,0,1);
	TH1D * mom_smb = new TH1D("mom_smb","mom_smb",100,0,1);

	TRandom3 * myRand = new TRandom3(0);
	for( int i = 0 ; i < nB ; i++ ){
		double dL = 270. + myRand->Rndm()*(310.-280.);
		double weight = (16.-6.)/(upperBedge-lowerBedge);

		// Fill some psuedo background level based on range we have available
		double tof_bkg = (lowerBedge + myRand->Rndm()*(upperBedge-lowerBedge))*(dL/100.);
		ToF_bkgrd->Fill( tof_bkg  );
		// Now in real life, we take the event in our background level we just generated
		// above, and re-drawn a random ToF based on the real range of interest, and
		// calculate momentum, etc..
		tof_bkg = (6. + myRand->Rndm()*(16.-6.))*dL/100.;
		double mom = 0.938/sqrt(1./pow(dL/(30.*tof_bkg),2) - 1.);
		mom_b->Fill( mom , weight);

		// To fill our signal+background spectrum, we want to re-draw a random ToF
		// in the nominal range so our background isn't exactly the same
		tof_bkg = (6. + myRand->Rndm()*(16.-6.))*dL/100.;
		mom = 0.938/sqrt(1./pow(dL/(30.*tof_bkg),2) - 1.);
		ToF_spb->Fill( tof_bkg , weight );
		mom_spb->Fill( mom , weight);
		mom_smb->Fill( mom );
	}

	for( int i = 0 ; i < nS ; i++ ){
		double dL = 270. + myRand->Rndm()*(310.-280.);
		double tof_sig = (tri->GetRandom())*dL/100.;
		double mom = 0.938/sqrt(1./pow(dL/(30.*tof_sig),2) - 1.);
		ToF_signal->Fill( tof_sig );
		ToF_spb->Fill( tof_sig );
		mom_spb->Fill( mom );
		mom_smb->Fill( mom );
		mom_s->Fill( mom );
	}
	

	mom_smb->Add(mom_b,-1);


	inFile->Close();
	outFile->cd();
	ToF_signal->Write();
	ToF_bkgrd->Write();
	ToF_spb->Write();
	mom_spb->Write();
	mom_smb->Write();
	mom_s->Write();
	mom_b->Write();
	outFile->Close();


	return 1;
}

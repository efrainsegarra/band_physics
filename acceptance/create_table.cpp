#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TH2D.h"
using namespace std;

// Beam energy:
const double Ebeam = 10.6;

// Functions
bool ParticleAcceptance( TVector3 particle , std::string type );

int main( int argc, char** argv){
	if( argc !=3 ){
		cerr << "Wrong number of arguments. Instead please use:\n\t./create_table [inputPhysicsRootFile] [outputRootFile]\n";
		return -1;
	}

	// Setup input file with tree
	TFile * inFile = new TFile(argv[1]);
	TTree * inTree = (TTree*)inFile->Get("skim");
	const int nEv = inTree->GetEntries();
		// Set branches
	int nHits;
	int sector, layer, component;
	double adcLcorr, adcRcorr;
	double meantimeFadc, meantimeTdc;
	double difftimeFadc, difftimeTdc;
	double x,y,z;
	double dL, ToF;
	double beta, pN_mag, theta_n, phi_n;
	double En, pN_cosTheta, phi_en;
	double CosTheta_qn, Xp, Wp, As, theta_qn;
	double adcLraw, adcRraw;
	double tTdcLraw, tTdcRraw;
	double tFadcLraw, tFadcRraw;
	double byHand_adcL, byHand_adcR;
	double byHand_meantimeFadc, byHand_meantimeTdc;
	double byHand_difftimeFadc, byHand_difftimeTdc;
	double byHand_x, byHand_y, byHand_z;
	double byHand_dL;
	double theta_e, phi_e;
	double theta_q, phi_q;
	double Q2, nu, xB, W2, q;
	double start_time;
	inTree->SetBranchAddress("nHits",		&nHits		);
	inTree->SetBranchAddress("sector",		&sector		);
	inTree->SetBranchAddress("layer",		&layer		);
	inTree->SetBranchAddress("component",		&component	);
	inTree->SetBranchAddress("adcLcorr",		&adcLcorr	);
	inTree->SetBranchAddress("adcRcorr",		&adcRcorr	);
	inTree->SetBranchAddress("meantimeFadc",	&meantimeFadc	);
	inTree->SetBranchAddress("meantimeTdc",		&meantimeTdc	);
	inTree->SetBranchAddress("difftimeFadc",	&difftimeFadc	);
	inTree->SetBranchAddress("difftimeTdc",		&difftimeTdc	);
	inTree->SetBranchAddress("x",			&x		);
	inTree->SetBranchAddress("y",			&y		);
	inTree->SetBranchAddress("z",			&z		);
	inTree->SetBranchAddress("dL",			&dL		);
	inTree->SetBranchAddress("ToF",			&ToF		);
	inTree->SetBranchAddress("beta",		&beta		);
	inTree->SetBranchAddress("pN_mag",		&pN_mag		);
	inTree->SetBranchAddress("theta_n",		&theta_n	);
	inTree->SetBranchAddress("phi_n",		&phi_n		);
	inTree->SetBranchAddress("En",			&En		);
	inTree->SetBranchAddress("pN_cosTheta",		&pN_cosTheta	);
	inTree->SetBranchAddress("phi_en",		&phi_en		);
	inTree->SetBranchAddress("CosTheta_qn",		&CosTheta_qn	);
	inTree->SetBranchAddress("Xp",			&Xp		);
	inTree->SetBranchAddress("Wp",			&Wp		);
	inTree->SetBranchAddress("As",			&As		);
	inTree->SetBranchAddress("theta_qn",		&theta_qn	);
	inTree->SetBranchAddress("adcLraw",		&adcLraw	);
	inTree->SetBranchAddress("adcRraw",		&adcRraw	);
	inTree->SetBranchAddress("tTdcLraw",		&tTdcLraw	);
	inTree->SetBranchAddress("tTdcRraw",		&tTdcRraw	);
	inTree->SetBranchAddress("tFadcLraw",		&tFadcLraw	);
	inTree->SetBranchAddress("tFadcRraw",		&tFadcRraw	);
	inTree->SetBranchAddress("byHand_adcL",		&byHand_adcL	);
	inTree->SetBranchAddress("byHand_adcR",		&byHand_adcR	);
	inTree->SetBranchAddress("byHand_meantimeFadc",	&byHand_meantimeFadc	);
	inTree->SetBranchAddress("byHand_meantimeTdc",	&byHand_meantimeTdc	);
	inTree->SetBranchAddress("byHand_difftimeFadc",	&byHand_difftimeFadc	);
	inTree->SetBranchAddress("byHand_difftimeTdc",	&byHand_difftimeTdc	);
	inTree->SetBranchAddress("byHand_x",		&byHand_x	);
	inTree->SetBranchAddress("byHand_y",		&byHand_y	);
	inTree->SetBranchAddress("byHand_z",		&byHand_z	);
	inTree->SetBranchAddress("byHand_dL",		&byHand_dL	);
	inTree->SetBranchAddress("theta_e",		&theta_e	);
	inTree->SetBranchAddress("phi_e",		&phi_e		);
	inTree->SetBranchAddress("theta_q",		&theta_q	);
	inTree->SetBranchAddress("phi_q",		&phi_q		);
	inTree->SetBranchAddress("Q2",			&Q2		);
	inTree->SetBranchAddress("nu",			&nu		);
	inTree->SetBranchAddress("xB",			&xB		);
	inTree->SetBranchAddress("W2",			&W2		);
	inTree->SetBranchAddress("q",			&q		);
	inTree->SetBranchAddress("STTime",		&start_time	);

	// Setup output file
	//ofstream outfile;
	TFile * outFile = new TFile(argv[2],"RECREATE");
	TH2D * HiMap = new TH2D("HiMap","HiMap",1000,0,10,1000,-1,-0);
	TH2D * LoMap = new TH2D("LoMap","LoMap",1000,0,10,1000,-1,-0);
	TH2D * HiCnt = new TH2D("HiCnt","HiCnt",1000,0,10,1000,-1,-0);
	TH2D * LoCnt = new TH2D("LoCnt","LoCnt",1000,0,10,1000,-1,-0);
	
	// Setup random generator
	TRandom3 * rand = new TRandom3(0);
	
	// Loop over all events
	for( int event = 0 ; event < nEv ; event++){
		nHits,sector,layer,component,adcLcorr,adcRcorr			= 0.;
		meantimeFadc,meantimeTdc,difftimeFadc,difftimeTdc 		= 0.;
		x,y,z,dL,ToF,beta,pN_mag,theta_n,phi_n,En,pN_cosTheta 		= 0.;
		phi_en,CosTheta_qn,Xp,Wp,As,theta_qn 				= 0.;
		theta_e,phi_e,theta_q,phi_q,Q2,nu,xB,W2,q 			= 0.;
		adcLraw, adcRraw, tTdcLraw, tTdcRraw, tFadcLraw, tFadcRraw 	= 0.;
		byHand_adcL, byHand_adcR, byHand_meantimeFadc			= 0.;
		byHand_meantimeTdc, byHand_difftimeFadc, byHand_difftimeTdc	= 0.;
		byHand_x, byHand_y, byHand_z, byHand_dL				= 0.;

		inTree->GetEntry(event);
		if( event % 1000 == 0 ) cerr << "Working on event " << event << "\n";
	
		// Generate a bunch of random phi_e, phi_qn for this event:
		const int nDraws = 1000;
		double acceptance = 0.;
		for( int i = 0; i < nDraws; i++){
			double newPhi_e = -M_PI + (2.*M_PI)*rand->Rndm();
			double newPhi_qn = -M_PI + (2.*M_PI)*rand->Rndm();

			// For each new phi_e, re-create an electron vector and ask if it is accepted or not
			TVector3 eVec;
			double mag = Ebeam - nu;
			eVec.SetMagThetaPhi( mag, theta_e, newPhi_e );
	
			bool eAccepted = ParticleAcceptance( eVec, "electron" );

			// For each new phi_qn, re-create a neutron vector and ask if it is accepted or not
			TVector3 qVec;
			qVec.SetMagThetaPhi( q, theta_q, phi_q );

			TVector3 nVec;
			nVec.SetMagThetaPhi( pN_mag, theta_n, phi_n );
			
			nVec.Rotate( newPhi_qn , qVec );
			
			bool nAccepted = ParticleAcceptance( nVec , "neutron" );

			if( eAccepted && nAccepted ) acceptance++;
		}
		acceptance /= (double)nDraws;
	
		if( Xp > 0.5 ){
			HiMap->Fill( Q2, CosTheta_qn, acceptance );
			HiCnt->Fill( Q2, CosTheta_qn );
		}
		if( fabs(Xp - 0.3) < 0.05 ){
			LoMap->Fill( Q2, CosTheta_qn, acceptance );
			LoCnt->Fill( Q2, CosTheta_qn );
		}

		// We have some defined grid in Q2, Xp, CosTheta_nq, and Pn, and so count up events
		// that fall in each bin, correcting for acceptance of that bin.

	}

	// Now for each bin in my histogram, I need to divde the acceptance map by the number of times
	// I went to that bin
	for( int binx = 1 ; binx < LoCnt->GetXaxis()->GetNbins() ; binx++){
		for( int biny = 1 ; biny < LoCnt->GetYaxis()->GetNbins(); biny++){
			int bin = LoCnt->GetBin(binx,biny);
			double count = LoCnt->GetBinContent(bin);
			LoMap->SetBinContent(bin, LoMap->GetBinContent(bin)/count );
		}
	}
	for( int binx = 1 ; binx < HiCnt->GetXaxis()->GetNbins() ; binx++){
		for( int biny = 1 ; biny < HiCnt->GetYaxis()->GetNbins(); biny++){
			int bin = HiCnt->GetBin(binx,biny);
			double count = HiCnt->GetBinContent(bin);
			HiMap->SetBinContent(bin, HiMap->GetBinContent(bin)/count );
		}
	}

	
	// Clean up
	delete rand;

	// Close files
	inFile->Close();
	//outfile.close();
	outFile->cd();
	HiMap->Write();
	LoMap->Write();
	outFile->Close();

	return 1;
}

bool ParticleAcceptance( TVector3 particle , std::string type ){
	double theta = particle.Theta();
	double phi = particle.Phi();
	double mom = particle.Mag();
	if( type == "electron"){
		if( mom < 2 ) return false;
		if( mom > 8 ) return false;
		if( phi > M_PI ) return false;
		if( phi < -M_PI ) return false;
		if( theta > 35*M_PI/180. ) return false;
		if( theta < 5*M_PI/180. ) return false;
	}
	if( type == "neutron"){
		if( mom > 0.6 ) return false;
		if( mom < 0.2 ) return false;
		if( phi > M_PI ) return false;
		if( phi < -M_PI ) return false;
		if( theta > 178*M_PI/180. ) return false;
		if( theta < 150*M_PI/180. ) return false;
	}
	return true;
}

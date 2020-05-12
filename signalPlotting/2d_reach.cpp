#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TVectorT.h"

using namespace std;

void FillHists( TTree* inTree , double weight , TH2D* hist );

int main(int argc, char ** argv){
	if (argc != 4){
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\t./code [S+B File] [B File] [Out RootFile]\n";
		return -1;
	}
	TFile * inFile_spb = new TFile(argv[1]);
	TFile * inFile_bac = new TFile(argv[2]);

	TVectorT<double> * norm = (TVectorT<double> *)inFile_spb->Get("norm_0_0");
	TTree * inTree_spb = (TTree*)inFile_spb->Get("calib");
	TTree * inTree_bac = (TTree*)inFile_bac->Get("mixed");

	double bac_weight = (*norm)[1] / inTree_bac->GetEntries() ;

	TH2D * spb_AsXp = new TH2D("spb_AsXp","",6,0.2,0.8,4,1.2,1.6);
	TH2D * bac_AsXp = new TH2D("bac_AsXp","",6,0.2,0.8,4,1.2,1.6);
	TH2D * err_AsXp = new TH2D("err_AsXp","",6,0.2,0.8,4,1.2,1.6);

	FillHists( inTree_spb , 1.0 		, spb_AsXp );
	FillHists( inTree_bac , bac_weight	, bac_AsXp );

	for( int xbin = 1 ; xbin <= spb_AsXp->GetXaxis()->GetNbins() ; xbin++ ){
		for( int ybin = 1 ; ybin <= spb_AsXp->GetYaxis()->GetNbins() ; ybin++ ){
			double cnts_spb = spb_AsXp->GetBinContent(xbin,ybin);
			double cnts_bac = bac_AsXp->GetBinContent(xbin,ybin);
			double err = sqrt( cnts_spb + cnts_bac );
			err_AsXp->SetBinContent(xbin,ybin, 100*err/(cnts_spb - cnts_bac));
			cout << err/(cnts_spb-cnts_bac)*100 << " " 
				<< err_AsXp->GetXaxis()->GetBinCenter(xbin) << " " 
				<< err_AsXp->GetYaxis()->GetBinCenter(ybin) << "\n" ;
		}
	}

	TFile * outFile = new TFile(argv[3],"RECREATE");
	outFile->cd();
	err_AsXp->Write();
	spb_AsXp->Write();
	bac_AsXp->Write();
	outFile->Close();
	

	inFile_spb->Close();
	inFile_bac->Close();

	return 0;
}

void FillHists( TTree* inTree , double weight , TH2D* hist ){
	double p_e       		 = 0;
	double theta_e   		 = 0;
	double phi_e     		 = 0;
	double q         		 = 0;
	double theta_q   		 = 0;
	double phi_q     		 = 0;
	double Q2        		 = 0;
	double nu        		 = 0;
	double xB        		 = 0;
	double W2        		 = 0;
	double theta_n  		 = 0; 
	double phi_n    		 = 0; 
	double ToF      		 = 0; 
	double beta     		 = 0; 
	double p_n      		 = 0; 
	double phi_nq   		 = 0; 
	double theta_nq 		 = 0; 
	double E_n      		 = 0; 
	double CosTheta_nq 		 = 0;
	double Xp       		 = 0; 
	double Wp       		 = 0; 
	double As       		 = 0; 
	inTree->SetBranchAddress("p_e",		&p_e       		 );
	inTree->SetBranchAddress("theta_e",	&theta_e   		 );
	inTree->SetBranchAddress("phi_e",	&phi_e     		 );
	inTree->SetBranchAddress("q",		&q         		 );
	inTree->SetBranchAddress("theta_q",	&theta_q   		 );
	inTree->SetBranchAddress("phi_q",	&phi_q     		 );
	inTree->SetBranchAddress("Q2",		&Q2        		 );
	inTree->SetBranchAddress("nu",		&nu        		 );
	inTree->SetBranchAddress("xB",		&xB        		 );
	inTree->SetBranchAddress("W2",		&W2        		 );
	inTree->SetBranchAddress("theta_n",	&theta_n  		 ); 
	inTree->SetBranchAddress("phi_n",	&phi_n    		 ); 
	inTree->SetBranchAddress("ToF",		&ToF      		 ); 
	inTree->SetBranchAddress("beta",	&beta     		 ); 
	inTree->SetBranchAddress("p_n",		&p_n      		 ); 
	inTree->SetBranchAddress("phi_nq",	&phi_nq   		 ); 
	inTree->SetBranchAddress("theta_nq",	&theta_nq 		 ); 
	inTree->SetBranchAddress("E_n",		&E_n      		 ); 
	inTree->SetBranchAddress("CosTheta_nq",	&CosTheta_nq 	 	 );
	inTree->SetBranchAddress("Xp",		&Xp       		 ); 
	inTree->SetBranchAddress("Wp",		&Wp       		 ); 
	inTree->SetBranchAddress("As",		&As       		 ); 

	for( int ev = 0 ; ev < inTree->GetEntries() ; ev++ ){
		if( ev % 1000000 == 0 ) cout << "\ton event " << ev << "\n";
		p_e       		 = 0;
		theta_e   		 = 0;
		phi_e     		 = 0;
		q         		 = 0;
		theta_q   		 = 0;
		phi_q     		 = 0;
		Q2        		 = 0;
		nu        		 = 0;
		xB        		 = 0;
		W2        		 = 0;
		theta_n  		 = 0; 
		phi_n    		 = 0; 
		ToF      		 = 0; 
		beta     		 = 0; 
		p_n      		 = 0; 
		phi_nq   		 = 0; 
		theta_nq 		 = 0; 
		E_n      		 = 0; 
		CosTheta_nq 	 	 = 0;
		Xp       		 = 0; 
		Wp       		 = 0; 
		As       		 = 0; 

		inTree->GetEntry(ev);

		hist->Fill( Xp , As , weight );

	}

	return;
}

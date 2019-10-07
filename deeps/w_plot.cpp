#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TVectorT.h"
#include "TFitResultPtr.h"
#include "TStyle.h"
#include "../../deuteron_dis/include/constants.h"

using namespace std;

int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};
void LoadGlobalShift();
void LoadRunByRunShift();
void LoadGlobalShift2ndIter();
double FADC_GLOBSHIFT[600] = {0.};
double FADC_RUNBYRUNSHIFT[10000] = {0.};
double FADC_GLOBSHIFT2NDITER[600] = {0.};
char* getRunNumber( char* parse );
double getBeamEnergy( int runNum );


int main(int argc, char ** argv){
	if (argc < 2){
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\t./w_plot [outputRootfile] [inputDatafiles]\n";
		return -1;
	}

	TFile * outFile = new TFile(argv[1],"RECREATE");
	LoadGlobalShift();
	LoadRunByRunShift();
	LoadGlobalShift2ndIter();

	TH1D * hW = new TH1D("hW","hW",120,0.5,2);
	TH1D * hWp = new TH1D("hWp","hWp",120,0.5,2);
	
	// Loop over all the files that are given to me
	for( int i = 2 ; i < argc ; i++ ){
		TFile * inFile = new TFile(argv[i]);
		if (!(inFile->GetListOfKeys()->Contains("skim"))){
			cerr << "File has no entries\n";
			return -2;
		}
		TTree * inTree = (TTree*)inFile->Get("skim");

		int thisRun = atoi(getRunNumber(argv[i]));
		// Get beam energy from this run number:
		double fixed_Ebeam = getBeamEnergy(thisRun);

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
		double EoP       		 = 0;
		double STTime    		 = 0;
		int sector_e	  		 = 0;
		double vrt_x_e   		 = 0;
		double vrt_y_e   		 = 0;
		double vrt_z_e   		 = 0;
		double t_e       		 = 0;
		double beta_e    		 = 0;
		double chi2pid_e 		 = 0;
		int nHits     			 = 0;
		int sector    			 = 0;
		int layer     			 = 0;
		int component 			 = 0;
		double adcLcorr  		 = 0;
		double adcRcorr  		 = 0;
		double meantimeFadc		 = 0;
		double meantimeTdc		 = 0;
		double difftimeFadc		 = 0;
		double difftimeTdc		 = 0;
		double x        		 = 0; 
		double y        		 = 0; 
		double z        		 = 0; 
		double dL       		 = 0; 
		double theta_n  		 = 0; 
		double phi_n    		 = 0; 
		double ToF      		 = 0; 
		double beta     		 = 0; 
		double p_n      		 = 0; 
		double phi_nq   		 = 0; 
		double theta_nq 		 = 0; 
		double E_n      		 = 0; 
		double phi_en   		 = 0; 
		double CosTheta_nq 		 = 0;
		double Xp       		 = 0; 
		double Wp       		 = 0; 
		double As       		 = 0; 
		int nADC     			 = 0; 
		int nTDC     			 = 0; 
		double phaseCorr		 = 0; 
		double adcLraw  		 = 0; 
		double adcRraw  		 = 0; 
		double tTdcLraw 		 = 0; 
		double tTdcRraw 		 = 0; 
		double tFadcLraw		 = 0; 
		double tFadcRraw		 = 0; 
		double byHand_adcL		 = 0; 
		double byHand_adcR		 = 0; 
		double byHand_meantimeFadc 	 = 0;
		double byHand_meantimeTdc  	 = 0;
		double byHand_difftimeFadc 	 = 0;
		double byHand_difftimeTdc 	 = 0;
		double byHand_x 		 = 0; 
		double byHand_y 		 = 0; 
		double byHand_z 		 = 0; 
		double byHand_dL		 = 0; 
		inTree->SetBranchAddress("p_e",	&p_e       		 );
		inTree->SetBranchAddress("theta_e",	&theta_e   		 );
		inTree->SetBranchAddress("phi_e",	&phi_e     		 );
		inTree->SetBranchAddress("q",	&q         		 );
		inTree->SetBranchAddress("theta_q",	&theta_q   		 );
		inTree->SetBranchAddress("phi_q",	&phi_q     		 );
		inTree->SetBranchAddress("Q2",	&Q2        		 );
		inTree->SetBranchAddress("nu",	&nu        		 );
		inTree->SetBranchAddress("xB",	&xB        		 );
		inTree->SetBranchAddress("W2",	&W2        		 );
		inTree->SetBranchAddress("EoP",	&EoP       		 );
		inTree->SetBranchAddress("STTime",	&STTime    		 );
		inTree->SetBranchAddress("sector_e",	&sector_e  		 );
		inTree->SetBranchAddress("vrt_x_e",	&vrt_x_e   		 );
		inTree->SetBranchAddress("vrt_y_e",	&vrt_y_e   		 );
		inTree->SetBranchAddress("vrt_z_e",	&vrt_z_e   		 );
		inTree->SetBranchAddress("t_e",	&t_e       		 );
		inTree->SetBranchAddress("beta_e",	&beta_e    		 );
		inTree->SetBranchAddress("chi2pid_e",	&chi2pid_e 		 );
		inTree->SetBranchAddress("nHits",	&nHits     		 );
		inTree->SetBranchAddress("sector",	&sector    		 );
		inTree->SetBranchAddress("layer",	&layer     		 );
		inTree->SetBranchAddress("component",	&component 		 );
		inTree->SetBranchAddress("adcLcorr",	&adcLcorr  		 );
		inTree->SetBranchAddress("adcRcorr",	&adcRcorr  		 );
		inTree->SetBranchAddress("meantimeFadc",	&meantimeFadc		 );
		inTree->SetBranchAddress("meantimeTdc",	&meantimeTdc		 );
		inTree->SetBranchAddress("difftimeFadc",	&difftimeFadc		 );
		inTree->SetBranchAddress("difftimeTdc",	&difftimeTdc		 );
		inTree->SetBranchAddress("x",	&x        		 ); 
		inTree->SetBranchAddress("y",	&y        		 ); 
		inTree->SetBranchAddress("z",	&z        		 ); 
		inTree->SetBranchAddress("dL",	&dL       		 ); 
		inTree->SetBranchAddress("theta_n",	&theta_n  		 ); 
		inTree->SetBranchAddress("phi_n",	&phi_n    		 ); 
		inTree->SetBranchAddress("ToF",	&ToF      		 ); 
		inTree->SetBranchAddress("beta",	&beta     		 ); 
		inTree->SetBranchAddress("p_n",	&p_n      		 ); 
		inTree->SetBranchAddress("phi_nq",	&phi_nq   		 ); 
		inTree->SetBranchAddress("theta_nq",	&theta_nq 		 ); 
		inTree->SetBranchAddress("E_n",	&E_n      		 ); 
		inTree->SetBranchAddress("phi_en",	&phi_en   		 ); 
		inTree->SetBranchAddress("CosTheta_nq",	&CosTheta_nq 		 );
		inTree->SetBranchAddress("Xp",	&Xp       		 ); 
		inTree->SetBranchAddress("Wp",	&Wp       		 ); 
		inTree->SetBranchAddress("As",	&As       		 ); 
		inTree->SetBranchAddress("nADC",	&nADC     		 ); 
		inTree->SetBranchAddress("nTDC",	&nTDC     		 ); 
		inTree->SetBranchAddress("phaseCorr",	&phaseCorr		 ); 
		inTree->SetBranchAddress("adcLraw",	&adcLraw  		 ); 
		inTree->SetBranchAddress("adcRraw",	&adcRraw  		 ); 
		inTree->SetBranchAddress("tTdcLraw",	&tTdcLraw 		 ); 
		inTree->SetBranchAddress("tTdcRraw",	&tTdcRraw 		 ); 
		inTree->SetBranchAddress("tFadcLraw",	&tFadcLraw		 ); 
		inTree->SetBranchAddress("tFadcRraw",	&tFadcRraw		 ); 
		inTree->SetBranchAddress("byHand_adcL",	&byHand_adcL		 ); 
		inTree->SetBranchAddress("byHand_adcR",	&byHand_adcR		 ); 
		inTree->SetBranchAddress("byHand_meantimeFadc",	&byHand_meantimeFadc 	 );
		inTree->SetBranchAddress("byHand_meantimeTdc",	&byHand_meantimeTdc  	 );
		inTree->SetBranchAddress("byHand_difftimeFadc",	&byHand_difftimeFadc 	 );
		inTree->SetBranchAddress("byHand_difftimeTdc",	&byHand_difftimeTdc 	 );
		inTree->SetBranchAddress("byHand_x",	&byHand_x 		 ); 
		inTree->SetBranchAddress("byHand_y",	&byHand_y 		 ); 
		inTree->SetBranchAddress("byHand_z",	&byHand_z 		 ); 
		inTree->SetBranchAddress("byHand_dL",	&byHand_dL		 ); 

		cout << "Working on file: " << argv[i] << "\n";
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
				EoP       		 = 0;
				STTime    		 = 0;
				sector_e  		 = 0;
				vrt_x_e   		 = 0;
				vrt_y_e   		 = 0;
				vrt_z_e   		 = 0;
				t_e       		 = 0;
				beta_e    		 = 0;
				chi2pid_e 		 = 0;
				nHits     		 = 0;
				sector    		 = 0;
				layer     		 = 0;
				component 		 = 0;
				adcLcorr  		 = 0;
				adcRcorr  		 = 0;
				meantimeFadc		 = 0;
				meantimeTdc		 = 0;
				difftimeFadc		 = 0;
				difftimeTdc		 = 0;
				x        		 = 0; 
				y        		 = 0; 
				z        		 = 0; 
				dL       		 = 0; 
				theta_n  		 = 0; 
				phi_n    		 = 0; 
				ToF      		 = 0; 
				beta     		 = 0; 
				p_n      		 = 0; 
				phi_nq   		 = 0; 
				theta_nq 		 = 0; 
				E_n      		 = 0; 
				phi_en   		 = 0; 
				CosTheta_nq 		 = 0;
				Xp       		 = 0; 
				Wp       		 = 0; 
				As       		 = 0; 
				nADC     		 = 0; 
				nTDC     		 = 0; 
				phaseCorr		 = 0; 
				adcLraw  		 = 0; 
				adcRraw  		 = 0; 
				tTdcLraw 		 = 0; 
				tTdcRraw 		 = 0; 
				tFadcLraw		 = 0; 
				tFadcRraw		 = 0; 
				byHand_adcL		 = 0; 
				byHand_adcR		 = 0; 
				byHand_meantimeFadc 	 = 0;
				byHand_meantimeTdc  	 = 0;
				byHand_difftimeFadc 	 = 0;
				byHand_difftimeTdc 	 = 0;
				byHand_x 		 = 0; 
				byHand_y 		 = 0; 
				byHand_z 		 = 0; 
				byHand_dL		 = 0; 

				inTree->GetEntry(ev);

				if( nHits != 1 || meantimeFadc==0 ) continue;
		
				// Recalculate momentum based on shifts
				int barID = sector*100 + layer*10 + component;
				if( barID == 315 || barID == 336 || barID == 352 || barID == 413 || barID == 445 ) continue;
				double tof = (meantimeFadc - STTime) - FADC_GLOBSHIFT[barID] - FADC_RUNBYRUNSHIFT[thisRun] - FADC_GLOBSHIFT2NDITER[barID];
				double beta = dL / (tof*cAir);
				double p_n = mN / sqrt( 1./pow(beta,2) - 1. );
				double MeVee_cut = 8.;
				if( sqrt(adcLcorr*adcRcorr) < MeVee_cut*3000. ) continue;
			
				// Recalculate electron quantities based on new beam energy
				TVector3 beamVec(0,0,fixed_Ebeam);
				TVector3 eVec;	eVec.SetMagThetaPhi(p_e,theta_e,phi_e);
				TVector3 qVec;	qVec = beamVec - eVec;
				q 	= qVec.Mag();
				theta_q = qVec.Theta();
				phi_q 	= qVec.Phi();
				nu 	= fixed_Ebeam - sqrt( p_e*p_e + mE*mE );
				Q2 	= q*q - nu*nu;
				xB	= Q2 / (2.*mP*nu);
				W2 = mP*mP - Q2 + 2*nu*mP;

				// Recalculate neutron quantities
				TVector3 nVec;	nVec.SetMagThetaPhi(p_n,theta_n,phi_n);
				TVector3 norm_scatter = qVec.Cross( beamVec );
				norm_scatter 	= norm_scatter.Unit();

				TVector3 norm_reaction = qVec.Cross( nVec );
				norm_reaction 	= norm_reaction.Unit();

				phi_nq 	= norm_scatter.Angle( norm_reaction );
				theta_nq = nVec.Angle( qVec );
				TVector3 direction = norm_scatter.Cross(norm_reaction);
				if( direction.Z() > 0 ){ // this means the phi_rq should be between 0 and pi
				}
				else if( direction.Z() < 0 ){ // this means the phi_rq should be between -pi and 0
					phi_nq *= (-1);
				}
				E_n = sqrt( mN*mN + p_n*p_n );
				double W_primeSq = mD*mD - Q2 + mN*mN + 2.*mD*(nu-E_n) - 2.*nu*E_n + 2.*q*p_n*cos(theta_nq);
				Wp = sqrt(W_primeSq);


				// Fill histograms
				//if( xB > 0.3 ){
					hW->Fill( sqrt(W2) );
					hWp->Fill( Wp );
				//}

			

		} // end loop over events

		inFile->Close();
	}// end loop over files

	outFile->cd();
	hW->Write();
	hWp->Write();

		
	outFile->Close();


	return 0;
}

void LoadGlobalShift(){
	ifstream f;
	int sector, layer, component, barId;
	double pol0, height, mean, sig, temp;

	f.open("global_offset_fadc.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		FADC_GLOBSHIFT[barId] = mean;
		f >> temp;
		f >> temp;
	}
	f.close();
}

char* getRunNumber( char* parse ){
	char * parse_copy = (char*) malloc( strlen(parse)+1 );
	char * parsed;

	strcpy( parse_copy, parse );
	char * loop = strtok( parse_copy, ".");
	while(loop){
		char* equals_sign = strchr(loop, '_');
		if (equals_sign){
			*equals_sign = 0;
			equals_sign++;
			parsed = (char*)malloc(strlen(equals_sign) + 1);
			strcpy(parsed, equals_sign);
		}
		loop = strtok(NULL, ".");
	}
	free(parse_copy);
	return parsed;
}

void LoadRunByRunShift(){
	ifstream f;
	int runnum;
	double pol0, height, mean, sig, temp;

	f.open("runByrun_offset_fadc.txt");
	while(!f.eof()){
		f >> runnum;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		FADC_RUNBYRUNSHIFT[runnum] = mean;
		f >> temp;
		f >> temp;
	}
	f.close();
}

void LoadGlobalShift2ndIter(){
	ifstream f;
	int sector, layer, component, barId;
	double pol0, height, mean, sig, temp;

	f.open("global_offset_fadc_2ndIter.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		FADC_GLOBSHIFT2NDITER[barId] = mean;
		f >> temp;
		f >> temp;
	}
	f.close();
}

double getBeamEnergy( int runNum ){
	double thisEn = 0;
	
	if( runNum <= 6399 ) thisEn = 10.6;
	else{ thisEn = 10.2; }

	return thisEn;
}

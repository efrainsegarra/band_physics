// Hacky AF way to align bars due to the residual RF calibration throughout the run period.
// The idea:
// 	- since using photon peak, all bars are globally aligned to something.
// 	- use a single run to put all bars together (IGNORE 315,336,352,413,445)
// 	and do another alignment (to get a run-by-run offset) which should correct for 
// 	the missing RF calibration in some runs.

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
double FADC_GLOBSHIFT[600] = {0.};
double TDC_GLOBSHIFT[600] = {0.};
char* getRunNumber( char* parse );

int main(int argc, char ** argv){

	gStyle->SetOptFit(1);

	if (argc < 3){
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\t./single_photon [outputTxtfile] [outputRootfile] [inputDatafiles]\n";
		return -1;
	}

	TFile * outFile = new TFile(argv[2],"RECREATE");

	LoadGlobalShift();

	// Histograms for every single bar
	const int nRuns = argc-3;
	TH1D ** ToF_spec 	= new TH1D*[nRuns];
	TF1 ** ToF_fits		= new TF1*[nRuns];
	TF1 ** ToF_fits_it	= new TF1*[nRuns];
	TCanvas ** cRun		= new TCanvas*[nRuns];
	for(int run = 0; run < nRuns; run++){
		int thisRun = atoi(getRunNumber(argv[run+3]));
		cout << thisRun << "\n";
		ToF_spec[run] = new TH1D( Form("ToF_spec_%i",thisRun), Form("ToF_spec_%i",thisRun), 640, -15, 25);
	}


	// Loop over all the files that are given to me
	for( int i = 3 ; i < argc ; i++ ){
		TFile * inFile = new TFile(argv[i]);
		if (!(inFile->GetListOfKeys()->Contains("skim"))){
			cerr << "File has no entries\n";
			return -2;
		}
		TTree * inTree = (TTree*)inFile->Get("skim");


		int nHits;
		double adcLcorr, adcRcorr;
		double meantimeFadc, STTime, dL;
		int sector, layer, component;

		inTree->SetBranchAddress("nHits",			&nHits		);
		inTree->SetBranchAddress("adcLcorr",			&adcLcorr	);
		inTree->SetBranchAddress("adcRcorr",			&adcRcorr	);
		//inTree->SetBranchAddress("meantimeFadc",		&meantimeFadc	);
		inTree->SetBranchAddress("meantimeTdc",		&meantimeFadc	);
		inTree->SetBranchAddress("STTime",			&STTime		);
		inTree->SetBranchAddress("dL",				&dL		);
		inTree->SetBranchAddress("sector",			&sector		);
		inTree->SetBranchAddress("layer",			&layer		);
		inTree->SetBranchAddress("component",			&component	);

		cout << "Working on file: " << argv[i] << "\n";
		for( int ev = 0 ; ev < inTree->GetEntries() ; ev++ ){
			if( ev % 1000000 == 0 ) cout << "\ton event " << ev << "\n";

			nHits = 0;
			adcLcorr = 0; adcRcorr = 0; meantimeFadc = 0;
			STTime = 0; dL = 0;
			sector = 0; component = 0; layer = 0;
			inTree->GetEntry(ev);

			if( nHits != 1 || meantimeFadc==0 ) continue;

			double MeVee_cut = 0.;
			if( sqrt(adcLcorr*adcRcorr) < MeVee_cut*2500. ) continue;

			int barID = sector*100 + layer*10 + component;
			//if( barID == 315 || barID == 336 || barID == 352 || barID == 413 || barID == 445 ) continue;
			//if( barID == 314 ) continue;
			double tof = (meantimeFadc - STTime - dL/cAir) - TDC_GLOBSHIFT[barID];
			//double tof = (meantimeFadc - STTime - dL/cAir) - FADC_GLOBSHIFT[barID];
			// correct run-by-run shift:
			//tof = tof - RUN_SHIFT[run];
			// add back in path-length for actual tof:
			//tof = tof + dL/cAir;

			ToF_spec[i-3]->Fill(tof);
		} // end loop over events

		inFile->Close();
	}// end loop over files




	ofstream out_file;
	out_file.open(argv[1]);
	TCanvas * c0 = new TCanvas("c0","c0",900,900);
	c0 -> Print("results_runbyrun.pdf(");
	for(int run = 0; run < nRuns; run++){
		// Create canvas
		cRun[run] = new TCanvas(Form("Run %i",run+1),Form("Run %i",run+1),900,900);
		int thisRun = atoi(getRunNumber(argv[run+3]));


		if( ToF_spec[run]->Integral() == 0 ){
			out_file << thisRun <<  " " 
				<< 0 << " " << 0 << " " 
				<< 0 << " " << 0 << " "
				<< 0 << " " << -1 << "\n";
			continue;
		}


		// Get the min and max of the fit based on assuming 0.3ns resolution and the peak position
		double sig_guess = 0.3;
		double max = ToF_spec[run]->GetMaximum();
		double max_pos = ToF_spec[run]->GetXaxis()->GetBinCenter( ToF_spec[run]->GetMaximumBin() );

		double min_fit = max_pos - 5;
		double max_fit = max_pos + 2.*sig_guess;
		ToF_fits[run] = new TF1(Form("ToF_fits_%i",run),"pol0+gaus(1)",min_fit,max_fit);

		// Set parameters of the fit before fitting:
		double background_lvl = 0.;
		for( int i = 1; i < 25; i++){
			background_lvl += ToF_spec[run]->GetBinContent(i);
		}
		background_lvl /= 24;
		// background level:
		ToF_fits[run]->SetParameter(0,background_lvl);
		// constant of gaus:
		ToF_fits[run]->SetParameter(1,max);
		// mean of gaus:
		ToF_fits[run]->SetParameter(2,max_pos);
		// sigma of gaus:
		ToF_fits[run]->SetParameter(3,sig_guess);

		ToF_spec[run]->Fit(ToF_fits[run],"QESR");

		// Now do another iteration of fits with the current parameters only if the parameters are reasonable
		if( ToF_fits[run]->GetParameter(3) < 0.6 ){
			sig_guess = ToF_fits[run]->GetParameter(3);
			max_fit = ToF_fits[run]->GetParameter(2) - 5;
			min_fit = ToF_fits[run]->GetParameter(2) + 1.5*sig_guess;
			ToF_fits_it[run] = new TF1(Form("ToF_fits_%i_it",run),"pol0+gaus(1)",min_fit,max_fit);
			ToF_fits_it[run]->SetLineColor(4);
			ToF_fits_it[run]->SetParameter(0, ToF_fits[run]->GetParameter(0) );
			ToF_fits_it[run]->SetParameter(1, ToF_fits[run]->GetParameter(1) );
			ToF_fits_it[run]->SetParameter(2, ToF_fits[run]->GetParameter(2) );
			ToF_fits_it[run]->SetParameter(3, ToF_fits[run]->GetParameter(3) );
			ToF_spec[run]->Fit(ToF_fits_it[run],"QESR");

			out_file << (thisRun) << " "  
				<< ToF_fits_it[run]->GetParameter(0) << " " << ToF_fits_it[run]->GetParameter(1) << " " 
				<< ToF_fits_it[run]->GetParameter(2) << " " << ToF_fits_it[run]->GetParameter(3) << " "
				<< ToF_spec[run]->Integral() << " " << 1 << "\n";
		}
		else{
			out_file << (thisRun) << " "  
				<< ToF_fits[run]->GetParameter(0) << " " << ToF_fits[run]->GetParameter(1) << " " 
				<< ToF_fits[run]->GetParameter(2) << " " << ToF_fits[run]->GetParameter(3) << " "
				<< ToF_spec[run]->Integral() << " " << 0 << "\n";

		}

		cRun[run]->cd(1);
		ToF_spec[run]->Draw();
		//ToF_fits[run]->Draw("SAME");
		cRun[run]->Update();
		cRun[run]->Modified();
		cRun[run] -> Print("results_runbyrun.pdf");
		//cRun[run]->Write();
	}
	c0 -> Print("results_runbyrun.pdf)");
	out_file.close();

	outFile->cd();
	for(int run = 0; run < nRuns; run++){
		ToF_spec[run] -> Write();
	}
	outFile->Close();


	return 0;
}

void LoadGlobalShift(){
	ifstream f;
	int sector, layer, component, barId;
	double pol0, height, mean, sig, temp;

	f.open("global_offset_fadc-10082019.txt");
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
	f.open("global_offset_tdc_1stIter-04132020.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		TDC_GLOBSHIFT[barId] = mean;
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

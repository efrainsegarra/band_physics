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
void ReadPreviousShifts();
double PREV_FADC_POL0[600];
double PREV_FADC_HEIGHT[600];
double PREV_FADC_GLOBSHIFT[600];
double PREV_FADC_SIGMA[600];

int main(int argc, char ** argv){

	gStyle->SetOptFit(1);

	if (argc < 3){
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\t./single_photon [outputTxtfile] [outputRootfile] [inputDatafiles]\n";
		return -1;
	}

	TFile * outFile = new TFile(argv[2],"RECREATE");

	// Histograms for every single bar
	
	TH1D **** ToF_spec = new TH1D***[5];
	TF1  **** ToF_fits = new TF1 ***[5];
	TF1  **** ToF_fits_it = new TF1 ***[5];
	TCanvas **** cSLC = new TCanvas***[5];	
	for( int sector = 0; sector < 5; sector++){
		ToF_spec[sector] = new TH1D**[5];
		ToF_fits[sector] = new TF1 **[5];
		ToF_fits_it[sector] = new TF1 **[5];
		cSLC[sector] = new TCanvas**[5];
		for( int layer = 0; layer < 5; layer++){
			ToF_spec[sector][layer] = new TH1D*[7];
			ToF_fits[sector][layer] = new TF1 *[7];
			ToF_fits_it[sector][layer] = new TF1 *[7];
			cSLC[sector][layer] = new TCanvas*[7];
			for(int component = 0; component < slc[layer][sector]; component++){
				//ToF_spec[sector][layer][component] = new TH1D(Form("ToF_spec_%i_%i_%i",(sector+1),(layer+1),(component+1)),Form("ToF_spec_%i_%i_%i",(sector+1),(layer+1),(component+1)),640,30,70);
				ToF_spec[sector][layer][component] = new TH1D(Form("ToF_spec_%i_%i_%i",(sector+1),(layer+1),(component+1)),Form("ToF_spec_%i_%i_%i",(sector+1),(layer+1),(component+1)),640,250,290);
			}
		}
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
		inTree->SetBranchAddress("meantimeTdc",		&meantimeFadc	);
		//inTree->SetBranchAddress("meantimeFadc",		&meantimeFadc	);
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

				double tof = (meantimeFadc - STTime - dL/cAir);

				ToF_spec[sector-1][layer-1][component-1]->Fill(tof);
		} // end loop over events

		inFile->Close();
	}// end loop over files

	outFile->cd();

		
	ofstream out_file;
	out_file.open(argv[1]);
	ReadPreviousShifts();

	TCanvas * c0 = new TCanvas("c0","c0",900,900);
	c0 -> Print("results_single_photon_tdc.pdf(");
	for( int sector = 0; sector < 5; sector++){
		for( int layer = 0; layer < 5; layer++){
			for(int component = 0; component < slc[layer][sector]; component++){

				// Create canvas
				cSLC[sector][layer][component] = new TCanvas(Form("S%iL%iC%i",sector+1,layer+1,component+1),Form("Sector %i, Layer %i, Component %i",sector+1,layer+1,component+1),900,900);

				
				if( ToF_spec[sector][layer][component]->Integral() == 0 ){
					out_file << (sector+1) << " " << (layer+1) << " " << (component+1) << " " 
							<< 0 << " " << 0 << " " 
							<< 0 << " " << 0 << " "
								<< 0 << " " << -1 << "\n";
					continue;
				}

				// Get the min and max of the fit based on assuming 0.3ns resolution and the peak position
				double sig_guess = 0.3;
				
				double max = ToF_spec[sector][layer][component]->GetMaximum();
				if( sector == 3 || sector == 4) max = max*0.75;
				else{ max = max*0.5; }
				double max_pos = ToF_spec[sector][layer][component]->GetXaxis()->GetBinCenter( ToF_spec[sector][layer][component]->GetMaximumBin() );
				max_pos = ToF_spec[sector][layer][component]->GetXaxis()->GetBinCenter( ToF_spec[sector][layer][component]->FindFirstBinAbove(max) );
				
				// Set parameters of the fit before fitting:
				double background_lvl = 0.;
				for( int i = 1; i < 25; i++){
					background_lvl += ToF_spec[sector][layer][component]->GetBinContent(i);
				}
				background_lvl /= 24;

				// Just use previous definitions now
				int barID = (sector+1)*100 + (layer+1)*10 + (component+1);
				//background_lvl = PREV_FADC_POL0[barID];
				//max = PREV_FADC_HEIGHT[barID];
				//max_pos = PREV_FADC_GLOBSHIFT[barID];
				//sig_guess = PREV_FADC_SIGMA[barID];

				double min_fit = max_pos - 10;
				double max_fit = max_pos + 2.*sig_guess;
				ToF_fits[sector][layer][component] = new TF1(Form("ToF_fits_%i_%i_%i",sector,layer,component),"pol0+gaus(1)",min_fit,max_fit);
					// background level:
				ToF_fits[sector][layer][component]->SetParameter(0,background_lvl);
					// constant of gaus:
				ToF_fits[sector][layer][component]->SetParameter(1,max);
					// mean of gaus:
				ToF_fits[sector][layer][component]->SetParameter(2,max_pos);
					// sigma of gaus:
				ToF_fits[sector][layer][component]->SetParameter(3,sig_guess);

				ToF_spec[sector][layer][component]->Fit(ToF_fits[sector][layer][component],"QESR");

				// Now do another iteration of fits with the current parameters only if the parameters are reasonable
				if( ToF_fits[sector][layer][component]->GetParameter(3) < 0.6 ){
					sig_guess = ToF_fits[sector][layer][component]->GetParameter(3);
					max_fit = ToF_fits[sector][layer][component]->GetParameter(2) - 5;
					min_fit = ToF_fits[sector][layer][component]->GetParameter(2) + 1.5*sig_guess;
					ToF_fits_it[sector][layer][component] = new TF1(Form("ToF_fits_%i_%i_%i_it",sector,layer,component),"pol0+gaus(1)",min_fit,max_fit);
					ToF_fits_it[sector][layer][component]->SetLineColor(4);
					ToF_fits_it[sector][layer][component]->SetParameter(0, ToF_fits[sector][layer][component]->GetParameter(0) );
					ToF_fits_it[sector][layer][component]->SetParameter(1, ToF_fits[sector][layer][component]->GetParameter(1) );
					ToF_fits_it[sector][layer][component]->SetParameter(2, ToF_fits[sector][layer][component]->GetParameter(2) );
					ToF_fits_it[sector][layer][component]->SetParameter(3, ToF_fits[sector][layer][component]->GetParameter(3) );
					ToF_spec[sector][layer][component]->Fit(ToF_fits_it[sector][layer][component],"QESR");
		
					out_file << (sector+1) << " " << (layer+1) << " " << (component+1) << " " 
							<< ToF_fits_it[sector][layer][component]->GetParameter(0) << " " << ToF_fits_it[sector][layer][component]->GetParameter(1) << " " 
							<< ToF_fits_it[sector][layer][component]->GetParameter(2) << " " << ToF_fits_it[sector][layer][component]->GetParameter(3) << " "
								<< ToF_spec[sector][layer][component]->Integral() << " " << 1 << "\n";
				}
				else{
					out_file << (sector+1) << " " << (layer+1) << " " << (component+1) << " " 
							<< ToF_fits[sector][layer][component]->GetParameter(0) << " " << ToF_fits[sector][layer][component]->GetParameter(1) << " " 
							<< ToF_fits[sector][layer][component]->GetParameter(2) << " " << ToF_fits[sector][layer][component]->GetParameter(3) << " "
								<< ToF_spec[sector][layer][component]->Integral() << " " << 0 << "\n";

				}

				cSLC[sector][layer][component]->cd(1);
				ToF_spec[sector][layer][component]->Draw();
				//ToF_fits[sector][layer][component]->Draw("SAME");
				
				cSLC[sector][layer][component]->Update();
				cSLC[sector][layer][component]->Modified();
				cSLC[sector][layer][component] -> Print("results_single_photon_tdc.pdf");
				//cSLC[sector][layer][component]->Write();
			}
		}
	}
	c0 -> Print("results_single_photon_tdc.pdf)");

	out_file.close();
	outFile->Close();


	return 0;
}

void ReadPreviousShifts(){
	ifstream f;
	int sector, layer, component, barId;
	double pol0, height, mean, sig, temp;

	f.open("global_offset_fadc-10032019.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		f >> temp;
		f >> temp;
		PREV_FADC_POL0[barId] = pol0;
		PREV_FADC_HEIGHT[barId] = height;
		PREV_FADC_GLOBSHIFT[barId] = mean;
		PREV_FADC_SIGMA[barId] = sig;
		//cout << barId << " " <<  PREV_FADC_POL0[barId] << " " << PREV_FADC_HEIGHT[barId] << " " << PREV_FADC_GLOBSHIFT[barId] << " " << PREV_FADC_SIGMA[barId] << "\n";
	}
	f.close();
}

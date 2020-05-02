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
void LoadLROffsets();
void LoadPaddleOffsets();
void LoadLayerOffsets();

double FADC_GLOBSHIFT[600] = {0.};
double FADC_LROFFSET[600] = {0.};
double FADC_PADDLEOFFSET[600] = {0.};
double FADC_LAYEROFFSET[600] = {0.};
double FADC_RUNBYRUNSHIFT[10000] = {0.};

double TDC_GLOBSHIFT[600] = {0.};
double TDC_LROFFSET[600] = {0.};
double TDC_PADDLEOFFSET[600] = {0.};
double TDC_LAYEROFFSET[600] = {0.};
double TDC_RUNBYRUNSHIFT[10000] = {0.};

char* getRunNumber( char* parse );


int main(int argc, char ** argv){

	gStyle->SetOptFit(1);

	if (argc < 3){
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\t./single_photon [outputTxtfile] [outputRootfile] [inputDatafiles]\n";
		return -1;
	}

	TFile * outFile = new TFile(argv[2],"RECREATE");

	// Load timing shifts
	LoadGlobalShift();
	LoadRunByRunShift();
	LoadLROffsets();
	LoadPaddleOffsets();
	LoadLayerOffsets();

	// Histograms for every single bar
	TH1D * ToF_pre = new TH1D(Form("ToF_pre"),Form("ToF_pre"),1000,-10,90);
	TH1D * ToF_all = new TH1D(Form("ToF_all"),Form("ToF_all"),1000,-10,90);
	TH1D ** ToF_SB = new TH1D*[20];
	TH1D ** ToF_SB_m = new TH1D*[20];
	for( int scale = 1 ; scale < 20 ; scale++){
		ToF_SB[scale-1] = new TH1D( Form("ToF_SB_%i",scale), Form("ToF_SB_%i",scale), 3200,-100,300 );
		ToF_SB_m[scale-1] = new TH1D( Form("ToF_SB_m_%i",scale), Form("ToF_SB_m_%i",scale),4000,-20,100 );
	}

	TH1D * ToF_m = new TH1D(Form("ToF_m"),Form("ToF_m"),1200,-20,100);
	TH1D * ToF_perf = new TH1D(Form("ToF_perf"),Form("ToF_perf"),3200,-100,300);
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
				ToF_spec[sector][layer][component] = new TH1D(Form("ToF_spec_%i_%i_%i",(sector+1),(layer+1),(component+1)),Form("ToF_spec_%i_%i_%i",(sector+1),(layer+1),(component+1)),640,-15,25);
			}
		}
	}
	
	// Loop over all the files that are given to me
	for( int i = 3 ; i < argc ; i++ ){
		int thisRun = atoi(getRunNumber(argv[i]));
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


				int barID = sector*100 + layer*10 + component;
				//if( barID == 315 || barID == 336 || barID == 352 || barID == 413 || barID == 445 ) continue;
				//if( barID == 314 ) continue;
				//double tof = (meantimeFadc - STTime - dL/cAir) + FADC_PADDLEOFFSET[barID] + FADC_LAYEROFFSET[barID] + FADC_LROFFSET[barID]/2.;
				double tof = (meantimeFadc - STTime - dL/cAir) + TDC_PADDLEOFFSET[barID] + TDC_LAYEROFFSET[barID] + fabs(TDC_LROFFSET[barID])/2.;

				// For S:B plot
				for( int scale = 1 ; scale < 20 ; scale++){
					if( sqrt(adcLcorr*adcRcorr) > scale*2300. ){
						double tof_fix = (meantimeFadc - STTime - dL/cAir) 
								- TDC_GLOBSHIFT[barID] - TDC_RUNBYRUNSHIFT[thisRun];
						ToF_SB[scale-1]->Fill( tof_fix + dL/cAir );
						ToF_SB_m[scale-1]->Fill( (tof_fix+dL/cAir)/(dL/100.) );
					}
				}	

				// For plots of before-after given calibrations. Only very small E cut
				double MeVee_cut = 2;
				if( sqrt(adcLcorr*adcRcorr) < MeVee_cut*2300. ) continue;
				
				ToF_pre->Fill(tof-274+50);
				tof = (meantimeFadc - STTime - dL/cAir) - TDC_GLOBSHIFT[barID] - TDC_RUNBYRUNSHIFT[thisRun];
				ToF_all->Fill(tof);
				// For bar resolution in TDC:
				ToF_spec[sector-1][layer-1][component-1]->Fill(tof);
			

				// For performance plots
				MeVee_cut = 5.;
				if( sqrt(adcLcorr*adcRcorr) < MeVee_cut*2300. ) continue;
				ToF_m->Fill( (tof+dL/cAir)/(dL/100.) );
				ToF_perf->Fill( tof + dL/cAir );
		} // end loop over events

		inFile->Close();
	}// end loop over files

	outFile->cd();
	for( int scale = 1 ; scale < 20 ; scale++){
		ofstream tof_sb;
		ofstream tof_m_sb;
		tof_sb.open(Form("tof_sb_%i.txt",scale));
		tof_m_sb.open(Form("tof_m_sb_%i.txt",scale));
			
		for(int binx = 1; binx<ToF_SB[scale-1]->GetXaxis()->GetNbins(); binx++){
			tof_sb << ToF_SB[scale-1]->GetXaxis()->GetBinCenter(binx) << " " <<
				  ToF_SB[scale-1]->GetBinContent(binx) << "\n";
		}
		for(int binx = 1; binx<ToF_SB_m[scale-1]->GetXaxis()->GetNbins(); binx++){
			tof_m_sb << ToF_SB_m[scale-1]->GetXaxis()->GetBinCenter(binx) << " " <<
				    ToF_SB_m[scale-1]->GetBinContent(binx) << "\n";
		}
		tof_sb.close();
		tof_m_sb.close();


		ToF_SB[scale-1]->Write();
		ToF_SB_m[scale-1]->Write();
	}
	ofstream tof_pre, tof_all, tof_m, tof_perf;
	tof_pre.open("tof_pre.txt");
	for( int binx = 1; binx < ToF_pre->GetXaxis()->GetNbins() ; binx++ ){
		tof_pre << ToF_pre->GetXaxis()->GetBinCenter(binx) << " " << 
				ToF_pre->GetBinContent(binx) << "\n";
	}
	tof_pre.close();

	tof_all.open("tof_all.txt");
	for( int binx = 1; binx < ToF_all->GetXaxis()->GetNbins() ; binx++ ){
		tof_all << ToF_all->GetXaxis()->GetBinCenter(binx) << " " << 
				ToF_all->GetBinContent(binx) << "\n";
	}
	tof_all.close();


	tof_m.open("tof_m.txt");
	for( int binx = 1; binx < ToF_m->GetXaxis()->GetNbins() ; binx++ ){
		tof_m << ToF_m->GetXaxis()->GetBinCenter(binx) << " " << 
				ToF_m->GetBinContent(binx) << "\n";
	}
	tof_m.close();

	tof_perf.open("tof_perf.txt");
	for( int binx = 1; binx < ToF_perf->GetXaxis()->GetNbins() ; binx++ ){
		tof_perf << ToF_perf->GetXaxis()->GetBinCenter(binx) << " " << 
				ToF_perf->GetBinContent(binx) << "\n";
	}
	tof_perf.close();


	ToF_pre->Write();
	ToF_all->Write();
	ToF_m->Write();
	ToF_perf->Write();

	ofstream out_file;
	out_file.open(argv[1]);

	TCanvas * c0 = new TCanvas("c0","c0",900,900);
	c0 -> Print("results_single_photon_tdc_2ndIter.pdf(");
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
				double max_pos = ToF_spec[sector][layer][component]->GetXaxis()->GetBinCenter( ToF_spec[sector][layer][component]->GetMaximumBin() );
				
				double min_fit = max_pos - 5;
				double max_fit = max_pos + 2.*sig_guess;
				ToF_fits[sector][layer][component] = new TF1(Form("ToF_fits_%i_%i_%i",sector,layer,component),"pol0+gaus(1)",min_fit,max_fit);
				
				// Set parameters of the fit before fitting:
				double background_lvl = 0.;
				for( int i = 1; i < 25; i++){
					background_lvl += ToF_spec[sector][layer][component]->GetBinContent(i);
				}
				background_lvl /= 24;
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
				cSLC[sector][layer][component] -> Print("results_single_photon_tdc_2ndIter.pdf");
				//cSLC[sector][layer][component]->Write();
			}
		}
	}
	c0 -> Print("results_single_photon_tdc_2ndIter.pdf)");
	out_file.close();
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

void LoadRunByRunShift(){
	ifstream f;
	int runnum;
	double pol0, height, mean, sig, temp;

	f.open("runByrun_offset_fadc-10082019.txt");
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
	f.open("runByrun_offset_tdc_firstIter-04132020.txt");
	while(!f.eof()){
		f >> runnum;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		TDC_RUNBYRUNSHIFT[runnum] = mean;
		f >> temp;
		f >> temp;
	}
	f.close();
}
void LoadLROffsets(){
	ifstream f;
	int sector, layer, component, barId;
	double ftdc, tdc, temp;

	f.open("../../parameters/lr_offsets.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> tdc;
		f >> ftdc;
		FADC_LROFFSET[barId] = ftdc;
		TDC_LROFFSET[barId] = tdc;
		f >> temp;
		f >> temp;
	}
	f.close();
}

void LoadPaddleOffsets(){
	ifstream f;
	int sector, layer, component, barId;
	double off, temp;

	f.open("../../parameters/paddle_offsets.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> off;
		FADC_PADDLEOFFSET[barId] = off;
		f >> temp;
	}
	f.close();
	f.open("../../parameters/paddle_offsets_tdc.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> off;
		TDC_PADDLEOFFSET[barId] = off;
		f >> temp;
	}
	f.close();
}
void LoadLayerOffsets(){
	ifstream f;
	int sector, layer, component, barId;
	double off, temp;

	f.open("../../parameters/layer_offsets.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> off;
		FADC_LAYEROFFSET[barId] = off;
		f >> temp;
	}
	f.close();
	f.open("../../parameters/layer_offsets_tdc.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> off;
		TDC_LAYEROFFSET[barId] = off;
		f >> temp;
	}
	f.close();
}

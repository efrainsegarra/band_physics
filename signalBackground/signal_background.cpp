#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TVector3.h"
#include "TF1.h"
#include "TVectorT.h"
#include "TLatex.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

using namespace std;

double shifts[600] = {0.};
int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};
void LoadToFShifts();

int main(int argc, char ** argv){
	if (argc < 3){
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\t./finalHists [path/to/output/file] [path/to/data/files]\n";
		return -1;
	}

	TFile * outFile = new TFile(argv[1],"RECREATE");

	// ToF histograms:
	cout << "Init histograms\n";
	TH1D * ToF_overall = new TH1D("ToF_overall","ToF_overall",80,0,40);
	TH1D * ToF_overall_pre = new TH1D("ToF_overall_pre","ToF_overall_pre",80,0,40);
	TH1D ** ToF_full = new TH1D*[15];
	TH1D *** ToF_dis = new TH1D**[15];
	for( int i = 0 ; i < 15 ; i++){
		ToF_full[i] = new TH1D(Form("ToF_full_%iMeVee",(i+1)),Form("ToF_full_%iMeVee",(i+1)),320,-40,120);
		ToF_dis[i] = new TH1D*[6];
		for( int j = 0; j < 6; j++){
			ToF_dis[i][j] = new TH1D(Form("ToF_full_%iMeVee_%ixB",(i+1),j),Form("ToF_full_%iMeVee_%ixB",(i+1),j),80,0,40);	
		}
	}
	cout << "...done!\n";

	// Loop over all the files that are given to me
	LoadToFShifts();
	for( int i = 2 ; i < argc ; i++ ){
		TFile * inFile = new TFile(argv[i]);
		if (!(inFile->GetListOfKeys()->Contains("skim"))){
			cerr << "File has no entries\n";
			return -2;
		}
		TTree * inTree = (TTree*)inFile->Get("skim");


		double Q2, W2, xB;
		int nHits, nTDC, nADC;
		double adcLcorr, adcRcorr;
		double meantimeFadc, STTime, dL;
		int sector, layer, component;

		inTree->SetBranchAddress("Q2",		&Q2		);
		inTree->SetBranchAddress("W2",		&W2		);
		inTree->SetBranchAddress("xB",		&xB		);
		inTree->SetBranchAddress("nHits",	&nHits		);
		inTree->SetBranchAddress("nADC",	&nADC		);
		inTree->SetBranchAddress("nTDC",	&nTDC		);
		inTree->SetBranchAddress("sector",	&sector		);
		inTree->SetBranchAddress("byHand_adcL",	&adcLcorr	);
		inTree->SetBranchAddress("byHand_adcR",	&adcRcorr	);
		inTree->SetBranchAddress("byHand_meantimeFadc",&meantimeFadc	);
		inTree->SetBranchAddress("STTime",	&STTime		);
		inTree->SetBranchAddress("byHand_dL",		&dL		);
		inTree->SetBranchAddress("sector",	&sector		);
		inTree->SetBranchAddress("layer",	&layer		);
		inTree->SetBranchAddress("component",	&component	);

		cout << "Working on file: " << argv[i] << "\n";

		for( int ev = 0 ; ev < inTree->GetEntries() ; ev++ ){
			if( ev % 1000000 == 0 ) cout << "\ton event " << ev << "\n";

			Q2 = 0; W2 = 0; xB = 0; nHits = 0; nADC = 0; nTDC = 0;
			adcLcorr = 0; adcRcorr = 0; meantimeFadc = 0;
			STTime = 0; dL = 0;
			sector = 0; component = 0; layer = 0;

			inTree->GetEntry(ev);

			if( nHits != 1) continue;
			if( nADC != 2) continue;
			if( nTDC != 2) continue;
			

			int id = sector*100 + layer*10 + component;
			double tof = ( meantimeFadc - STTime - shifts[id] ) / (dL/100.) ;
			double adc = sqrt(adcLcorr*adcRcorr)/2500.;

			ToF_overall_pre->Fill( (meantimeFadc - STTime)/(dL/100.) );
			ToF_overall->Fill( tof );

			if( sector == 3 || sector == 4 ) continue;

			if( Q2 > 2 && sqrt(W2) > 2.2 ){
				for( int MeVeeCut = 1; MeVeeCut < 16; MeVeeCut++ ){
					if( adc > MeVeeCut ){
						ToF_full[MeVeeCut-1] -> Fill( tof );
						for( int j = 0; j < 6; j+=1 ){ // xB loop cut
							double xLo = 0.1+0.1*j;
							double xHi = 0.1+0.1*j+0.1;
							if( xB > xLo && xB < xHi ){
								ToF_dis[MeVeeCut-1][j] -> Fill( tof );
							} // end xB if
						} // end xB loop
					} // end MeV if
				} // end MeV loop
			} // end Q2/W if

		} // end loop over events

		inFile->Close();
	}// end loop over files


	TF1 * bkg = new TF1("bkg","pol0",-20,0);

	ofstream outtext;
	outtext.open("/home/segarrae/software/band/analysis/signalBackground/spb.txt");
	outFile->cd();

	ToF_overall_pre->Write();
	ToF_overall->Write();
	for( int i = 0 ; i < 15 ; i++){
		ToF_full[i]->Write();

		// Get integral of S+B in ToF spectrum:
		int binLo = ToF_full[i]->FindBin(6.); // roughly 600 MeV/c
		int binHi = ToF_full[i]->FindBin(16.); // roughly 200 MeV/c

		double spb = ToF_full[i]->Integral(binLo,binHi);
		TFitResultPtr ptr = ToF_full[i]->Fit("bkg","QESR","",-20,0);
		double b = ptr->Parameter(0)*(binHi-binLo);
		
		outtext << (i+1) << " " << spb << " " << b << "\n";

		for( int j = 0; j < 6; j+=1 ){ 
			ToF_dis[i][j]->Write();
		}
	}
	outtext.close();
	outFile->Close();

}

void LoadToFShifts(){
	ifstream in_file;
	in_file.open("photon_params_fadc.txt");
	char line[256]; 
	int sector, layer, component;
	double bkg_lvl, cons, mean, sigma;
	int numEv, good;
	if (in_file.is_open() == true){ 
		while (!in_file.eof()) { 
			in_file >> sector >> layer >> component >> bkg_lvl >> cons >> mean >> sigma >> numEv >> good;
			int id = sector*100 + layer*10 + component;
			shifts[id] = mean;
		} 
	} 
	return;
}

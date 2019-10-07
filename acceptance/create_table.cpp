#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "../../deuteron_dis/include/constants.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TH2D.h"
using namespace std;

// Beam energy:
TVector3 e0(0,0,Ebeam);

// Functions
bool ParticleAcceptance( TVector3 particle , std::string type );
double kinLimit( double *x , double *p);
void normalizeMap( TH2D * count_map, TH2D * normMe );

int main( int argc, char** argv){
	if( argc <3 ){
		cerr << "Wrong number of arguments. Instead please use:\n\t./create_table [outputRootFile] [inputPhysicsRootFiles]\n";
		return -1;
	}


	// Setup output file
	//ofstream outfile;
	TFile * outFile = new TFile(argv[1],"RECREATE");
	TTree * outTree = new TTree("AccOut","Output Tree");
	double new_phi_e, new_phi_n, new_phi_q;
	outTree->Branch("phi_e",&new_phi_e);
	outTree->Branch("phi_n",&new_phi_n);
	outTree->Branch("phi_q",&new_phi_q);
	TH2D ** HiMap = new TH2D*[50];
	TH2D ** HiCnt = new TH2D*[50];
	TH2D ** LoMap = new TH2D*[50];
	TH2D ** LoCnt = new TH2D*[50];
	for( int i = 0 ; i < 50; i++){
		HiMap[i] = new TH2D(Form("HiMap_%i",i),"",100,0,10,100,-1,-0.4);
		HiCnt[i] = new TH2D(Form("HiCnt_%i",i),"",100,0,10,100,-1,-0.4);
		LoMap[i] = new TH2D(Form("LoMap_%i",i),"",100,0,10,100,-1,-0.4);
		LoCnt[i] = new TH2D(Form("LoCnt_%i",i),"",100,0,10,100,-1,-0.4);
	}	
	
	// Setup random generator
	TRandom3 * rand = new TRandom3(0);

	for( int file = 2 ; file < argc ; file++ ){
		// Setup input file with tree
		TFile * inFile = new TFile(argv[file]);
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
			if( event % 10000 == 0 ) cerr << "Working on event " << event << " (" << (double)event/nEv*100 << "%" << " done)" << "\n";
			//if( event > 100000 ) break;
		
			// Generate a bunch of random phi_e, phi_qn for this event:
			int momBin = (int)((pN_mag - 0.20)/0.01);
			if( momBin < 0 || momBin > 49) continue;

			const int nDraws = 100;
			double acceptance = 0.;

			// Only for our events in hix/lowx, do the map procedure
			if( fabs(Xp - 0.505) < 0.005 || fabs(Xp - 0.3) < 0.005 ){
		
				for( int i = 0; i < nDraws; i++){
					
					TVector3 eVec, qVec, nVec;
					eVec.SetMagThetaPhi( Ebeam-nu , theta_e, phi_e );
					qVec.SetMagThetaPhi( q, theta_q, phi_q );
					nVec.SetMagThetaPhi( pN_mag, theta_n, phi_n );
					
					// With this new phi_e, I need to rotate ALL the vectors
					double newPhi_e = -M_PI + (2.*M_PI)*rand->Rndm();
					eVec.Rotate( newPhi_e , e0.Unit() );
					qVec.Rotate( newPhi_e , e0.Unit() );
					nVec.Rotate( newPhi_e , e0.Unit() );

					// Now I pull a new random phi_qn, and rotate just the neutron vector
					// around q vector by that number
					double newPhi_qn = -M_PI + (2.*M_PI)*rand->Rndm();
					nVec.Rotate( newPhi_qn , qVec.Unit() );
					
					// Ask if both new draws are accepted or not
					bool eAccepted = ParticleAcceptance( eVec, "electron" );
					bool nAccepted = ParticleAcceptance( nVec , "neutron" );

					if( eAccepted && nAccepted ) acceptance++;
					
					//outTree->Fill();
				}
				acceptance /= (double)nDraws;

				// Calculate phi_nq properly and cut
				// around it to simulate bins in phi_rq
				TVector3 eVec, qVec, nVec;
				TVector3 beamVec(0,0,Ebeam);
				eVec.SetMagThetaPhi( Ebeam-nu , theta_e, phi_e );
				qVec.SetMagThetaPhi( q, theta_q, phi_q );
				nVec.SetMagThetaPhi( pN_mag, theta_n, phi_n );

				TVector3 norm_scatter = qVec.Cross( beamVec );
				norm_scatter 	= norm_scatter.Unit();
				TVector3 norm_reaction = qVec.Cross( nVec );
				norm_reaction 	= norm_reaction.Unit();

				double phi_rq 	= norm_scatter.Angle( norm_reaction );

				TVector3 direction = norm_scatter.Cross(norm_reaction);
				if( direction.Z() > 0 ){ // this means the phi_rq should be between 0 and pi
				}
				else if( direction.Z() < 0 ){ // this means the phi_rq should be between -pi and 0
					phi_rq *= (-1);
				}
			
				// In future we need Xp and Phi bins!
				if( fabs(Xp - 0.505) < 0.005 ){
					HiMap[momBin]->Fill( Q2, CosTheta_qn, acceptance );
					HiCnt[momBin]->Fill( Q2, CosTheta_qn );
				}
				if( fabs(Xp - 0.3) < 0.005 ){
					LoMap[momBin]->Fill( Q2, CosTheta_qn, acceptance );
					LoCnt[momBin]->Fill( Q2, CosTheta_qn );
				}

			} // end if on hix/lox events
		}// end loop over events in file

		// Close file
		inFile->Close();
	} // end loop over files
	
	// Now for each bin in my histogram, I need to divde the acceptance map by the number of times
	// I went to that bin
	for( int i = 0 ; i < 50; i++){
		// first loop through to normalize the acceptance map
		normalizeMap( HiCnt[i], HiMap[i] );
		normalizeMap( LoCnt[i], LoMap[i] );
		// second loop through to fill in the 0s by averaging nearby bin
		//fixZeros( HiCnt[i], HiMap[i] );
		//fixZeros( LoCnt[i], LoMap[i] );
		
	}

	// Print largest theta ranges drawn in random map generator to make sure my generated
	// cross section covers all these possible values
	//cout << "Min/Max Electron Theta: " << smallest_theta_e*180./M_PI << " " << largest_theta_e*180./M_PI << "\n";
	//cout << "Min/Max Neutron Theta: " << smallest_theta_n*180./M_PI << " " << largest_theta_n*180./M_PI << "\n";

	
	// Clean up
	delete rand;

	//outfile.close();
	outFile->cd();
	outTree->Write();
	
	TCanvas ** draws = new TCanvas*[50];
	TF1 ** model_hi = new TF1*[50];
	TF1 ** model_lo = new TF1*[50];
	for( int i = 0 ; i < 50; i++){
		// Create x' hi kinematic edge
		double p_n = 0.2 + (i)*0.01;
		model_hi[i] = new TF1(Form("edge_hi_%i",i),kinLimit,0,10,3);
		model_hi[i]->SetParameters(0.9,Ebeam,p_n);
			// write
		HiMap[i]->Write();
		model_hi[i]->Write();
	
		// Create x' lo kinematic edge
		p_n = 0.2 + (i)*0.01;
		model_lo[i] = new TF1(Form("edge_lo_%i",i),kinLimit,0,10,3);
		model_lo[i]->SetParameters(0.35,Ebeam,p_n);
		
		LoMap[i]->Write();
		model_lo[i]->Write();

	}
	outFile->Close();
		
	return 1;
}

bool ParticleAcceptance( TVector3 particle , std::string type ){
	double theta = particle.Theta();
	double phi = particle.Phi();
	double mom = particle.Mag();
	if( type == "e-"){
		if( mom < 2 ) return false;
		if( mom > 8 ) return false;
		if( phi < -175.*M_PI/180. ) return false;
		if( phi < -115.*M_PI/180. && phi > -125.*M_PI/180. ) return false;
		if( phi < -55.*M_PI/180. && phi > -65.*M_PI/180. ) return false;
		if( phi < 5.*M_PI/180. && phi > -5.*M_PI/180. ) return false;
		if( phi < 65.*M_PI/180. && phi > 55.*M_PI/180. ) return false;
		if( phi < 125.*M_PI/180. && phi > 115.*M_PI/180. ) return false;
		if( phi > 175.*M_PI/180. ) return false;
		if( theta > 35*M_PI/180. ) return false;
		if( theta < 8*M_PI/180. ) return false;
	}
	if( type == "neutron"){
		if( mom > 0.550 ) return false;
		if( mom < 0.250 ) return false;
		if( phi < -30.*M_PI/180. && phi > -150.*M_PI/180. ) return false;
		if( theta > 178*M_PI/180. ) return false;
		if( theta < 150*M_PI/180. ) return false;
	}
	return true;
}

double kinLimit( double *x , double *p){
	double var = *x;
	// p[0] = x'
	// p[1] = nu
	// p[2] = p_n
	double num = - var + 2*mD*p[0]*p[1] - 2.*sqrt(mN*mN + p[2]*p[2])*p[0]*p[1];
	double denom = 2.*p[2]*p[0]*sqrt( var*var + p[1]*p[1]);
	return - num / denom ;
}

void normalizeMap( TH2D * count_map, TH2D * normMe ){
	for( int binx = 1 ; binx < count_map->GetXaxis()->GetNbins() ; binx++){
		for( int biny = 1 ; biny < count_map->GetYaxis()->GetNbins(); biny++){
			int bin = count_map->GetBin(binx,biny);
			double count = count_map->GetBinContent(bin);
			if( count == 0 ) continue;
			normMe->SetBinContent(bin, normMe->GetBinContent(bin)/count );
			normMe->GetZaxis()->SetRangeUser(0,1);
		}
	}
	return;
}

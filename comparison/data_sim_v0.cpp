#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TVectorT.h"
#include "TVector3.h"
#include "../../deuteron_dis/include/constants.h"

using namespace std;

void FillHistsSim( TTree* inTree , double weight , TH1D* h_pn, TH1D** hs_xp , TH1D * h_ashi, TH1D* h_aslo , TH1D* h_xpfull , TH1D** hs_xpnu , TH1D * h_nuhi, TH1D* h_nulo );
void FillHists( TTree* inTree , double weight , TH1D* h_pn, TH1D** hs_xp , TH1D * h_ashi, TH1D* h_aslo , TH1D* h_xpfull , TH1D ** hs_xpnu , TH1D * h_nuhi, TH1D* h_nulo );
void Pretty1D( TH1D * hist , int color=4 , int linewidth=2 , double scale=1 );

int main(int argc, char ** argv){
	if (argc != 5){
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\t./code [S+B File] [B File] [Sim File] [Out RootFile]\n";
		return -1;
	}
	TFile * inFile_spb = new TFile(argv[1]);
	TFile * inFile_bac = new TFile(argv[2]);
	TFile * inFile_sim = new TFile(argv[3]);

	TVectorT<double> * norm = (TVectorT<double> *)inFile_spb->Get("norm_0_0");
	TTree * inTree_spb = (TTree*)inFile_spb->Get("calib");
	TTree * inTree_bac = (TTree*)inFile_bac->Get("mixed");
	TTree * inTree_sim = (TTree*)inFile_sim->Get("T");

	double bac_weight = (*norm)[1] / inTree_bac->GetEntries() ;
	
	TH1D * spb_pn = new TH1D("spb_pn","",15,0.2,0.5);
	TH1D * bac_pn = new TH1D("bac_pn","",15,0.2,0.5);
	TH1D * sim_pn = new TH1D("sim_pn","",15,0.2,0.5);


	TH1D * spb_xpfull = new TH1D("spb_xpfull","",50,0,1);
	TH1D * bac_xpfull = new TH1D("bac_xpfull","",50,0,1);
	TH1D * sim_xpfull = new TH1D("sim_xpfull","",50,0,1);
	TH1D ** spb_xpas = new TH1D*[5];
	TH1D ** bac_xpas = new TH1D*[5];
	TH1D ** sim_xpas = new TH1D*[5];
	for( int i = 0; i < 5 ; i++){
		spb_xpas[i] = new TH1D(Form("spb_xpas_%i",i),"",50,0,1);
		bac_xpas[i] = new TH1D(Form("bac_xpas_%i",i),"",50,0,1);
		sim_xpas[i] = new TH1D(Form("sim_xpas_%i",i),"",50,0,1);
	}
	TH1D * spb_asHi = new TH1D("spb_asHi","",5,1.2,1.7);
	TH1D * spb_asLo = new TH1D("spb_asLo","",5,1.2,1.7);
	TH1D * bac_asHi = new TH1D("bac_asHi","",5,1.2,1.7);
	TH1D * bac_asLo = new TH1D("bac_asLo","",5,1.2,1.7);
	TH1D * sim_asHi = new TH1D("sim_asHi","",5,1.2,1.7);
	TH1D * sim_asLo = new TH1D("sim_asLo","",5,1.2,1.7);
	TH1D ** spb_xpnu = new TH1D*[5];
	TH1D ** bac_xpnu = new TH1D*[5];
	TH1D ** sim_xpnu = new TH1D*[5];
	for( int i = 0; i < 5 ; i++){
		spb_xpnu[i] = new TH1D(Form("spb_xpnu_%i",i),"",50,0,1);
		bac_xpnu[i] = new TH1D(Form("bac_xpnu_%i",i),"",50,0,1);
		sim_xpnu[i] = new TH1D(Form("sim_xpnu_%i",i),"",50,0,1);
	}
	TH1D * spb_nuHi = new TH1D("spb_nuHi","",4,-0.9,-0.1);
	TH1D * spb_nuLo = new TH1D("spb_nuLo","",4,-0.9,-0.1);
	TH1D * bac_nuHi = new TH1D("bac_nuHi","",4,-0.9,-0.1);
	TH1D * bac_nuLo = new TH1D("bac_nuLo","",4,-0.9,-0.1);
	TH1D * sim_nuHi = new TH1D("sim_nuHi","",4,-0.9,-0.1);
	TH1D * sim_nuLo = new TH1D("sim_nuLo","",4,-0.9,-0.1);

	FillHists( inTree_spb , 1.0 		, spb_pn , spb_xpas , spb_asHi , spb_asLo , spb_xpfull , spb_xpnu , spb_nuHi , spb_nuLo);
	FillHists( inTree_bac , bac_weight	, bac_pn , bac_xpas , bac_asHi , bac_asLo , bac_xpfull , bac_xpnu , bac_nuHi , bac_nuLo);
	FillHistsSim( inTree_sim , 1.0		, sim_pn , sim_xpas , sim_asHi , sim_asLo , sim_xpfull , sim_xpnu , sim_nuHi , sim_nuLo);

	TFile * outFile = new TFile(argv[4],"RECREATE");
	outFile->cd();
	spb_pn->Write();
	bac_pn->Write();
	sim_pn->Write();
	
	Pretty1D( spb_xpfull );
	Pretty1D( bac_xpfull , 2 );
	Pretty1D( sim_xpfull , 3 );
	spb_xpfull->Write();
	bac_xpfull->Write();
	sim_xpfull->Write();


	for( int i = 0; i < 5 ; i++){
		Pretty1D( spb_xpas[i] );
		Pretty1D( bac_xpas[i] , 2 );
		Pretty1D( sim_xpas[i] , 3 );
		Pretty1D( spb_xpnu[i] );
		Pretty1D( bac_xpnu[i] , 2 );
		Pretty1D( sim_xpnu[i] , 3 );

		spb_xpas[i]->Write();
		bac_xpas[i]->Write();
		sim_xpas[i]->Write();
		spb_xpnu[i]->Write();
		bac_xpnu[i]->Write();
		sim_xpnu[i]->Write();
	}
	Pretty1D(spb_asHi);
	Pretty1D(bac_asHi, 2 );
	Pretty1D(sim_asHi, 3 );
	Pretty1D(spb_asLo);
	Pretty1D(bac_asLo, 2 );
	Pretty1D(sim_asLo, 3 );
	Pretty1D(spb_nuHi);
	Pretty1D(bac_nuHi, 2 );
	Pretty1D(sim_nuHi, 3 );
	Pretty1D(spb_nuLo);
	Pretty1D(bac_nuLo, 2 );
	Pretty1D(sim_nuLo, 3 );
	
	spb_asHi->Write();
	spb_asLo->Write();
	bac_asHi->Write();
	bac_asLo->Write();
	sim_asHi->Write();
	sim_asLo->Write();
	spb_nuHi->Write();
	spb_nuLo->Write();
	bac_nuHi->Write();
	bac_nuLo->Write();
	sim_nuHi->Write();
	sim_nuLo->Write();
	outFile->Close();
	

	inFile_spb->Close();
	inFile_bac->Close();
	inFile_sim->Close();

	return 0;
}

void Pretty1D( TH1D * hist , int color , int linewidth , double scale ){
	hist->SetLineColor(color);
	hist->SetLineWidth(linewidth);
	hist->SetStats(0);
	hist->Scale(scale);
	return;
}

void FillHists( TTree* inTree , double weight , TH1D* h_pn, TH1D** hs_xp , TH1D * h_ashi, TH1D* h_aslo , TH1D* h_xpfull , TH1D ** hs_xpnu , TH1D * h_nuhi, TH1D* h_nulo ){
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
	double virt			 = 0;
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
		virt			 = 0;

		inTree->GetEntry(ev);

		virt = ( pow(mD-E_n,2) - p_n*p_n - mP*mP )/(mP*mP);
		if( Wp > 1.8 ) {
			int Nu_bin = (int) ((virt+0.9)/0.2);
			if( Nu_bin >= 0 || Nu_bin < 4 ){
				hs_xpnu[Nu_bin]->Fill( Xp , weight );
				if( Xp > 0.6 && Xp < 1 ){
					h_nuhi->Fill( virt , weight );
				}
				if( fabs(Xp-0.3)<0.05) {
					h_nulo->Fill( virt , weight );
				}
			}
			h_xpfull->Fill( Xp , weight );


			int As_bin = (int) ((As-1.2)/0.1);
			if( As_bin < 0 || As_bin > 4 ) continue;
			hs_xp[As_bin]->Fill( Xp , weight );


			if( Xp > 0.6 && Xp < 1 ){
				h_ashi->Fill( As , weight );
			}


			if( fabs(Xp-0.3)<0.05) {
				h_pn->Fill( p_n , weight );
				h_aslo->Fill( As , weight );
			}
		}

	}

	return;
}
void FillHistsSim( TTree* inTree , double weight , TH1D* h_pn, TH1D** hs_xp , TH1D * h_ashi, TH1D* h_aslo , TH1D* h_xpfull , TH1D** hs_xpnu , TH1D * h_nuhi, TH1D* h_nulo ){
	double p_e       		 = 0;
	double theta_e   		 = 0;
	double phi_e     		 = 0;
	double p_n      		 = 0; 
	double theta_n  		 = 0; 
	double phi_n    		 = 0; 
	double Q2        		 = 0;
	double xB        		 = 0;
	double W2        		 = 0;
	double theta_q   		 = 0;
	double phi_q     		 = 0;
	inTree->SetBranchAddress("p_e",		&p_e       		 );
	inTree->SetBranchAddress("theta_e",	&theta_e   		 );
	inTree->SetBranchAddress("phi_e",	&phi_e     		 );
	inTree->SetBranchAddress("theta_q",	&theta_q   		 );
	inTree->SetBranchAddress("phi_q",	&phi_q     		 );
	inTree->SetBranchAddress("Q2",		&Q2        		 );
	inTree->SetBranchAddress("xB",		&xB        		 );
	inTree->SetBranchAddress("W2",		&W2        		 );
	inTree->SetBranchAddress("theta_n",	&theta_n  		 ); 
	inTree->SetBranchAddress("phi_n",	&phi_n    		 ); 
	inTree->SetBranchAddress("p_n",		&p_n      		 ); 

	double q         		 = 0;
	double nu        		 = 0;
	double ToF      		 = 0; 
	double beta     		 = 0; 
	double phi_nq   		 = 0; 
	double theta_nq 		 = 0; 
	double E_n      		 = 0; 
	double CosTheta_nq 		 = 0;
	double Xp       		 = 0; 
	double Wp       		 = 0; 
	double As       		 = 0; 
	double virt			 = 0;

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
		virt 			 = 0; 

		inTree->GetEntry(ev);
		theta_e *= M_PI/180;
		theta_n *= M_PI/180;
		phi_e *= M_PI/180;
		phi_n *= M_PI/180;

		double Ebeam = 10.6;
		TVector3 beamVec(0,0,Ebeam);
		TVector3 eVec;	eVec.SetMagThetaPhi(p_e,theta_e,phi_e);
		TVector3 qVec;	qVec = beamVec - eVec;
		q 	= qVec.Mag();
		theta_q = qVec.Theta();
		phi_q 	= qVec.Phi();
		nu 	= Ebeam - sqrt( p_e*p_e + mE*mE );
		Q2 	= q*q - nu*nu;
		xB	= Q2 / (2.*mP*nu);
		W2	= mP*mP - Q2 + 2*nu*mP;

		TVector3 nVec;	nVec.SetMagThetaPhi(1.,theta_n,phi_n);

		// Calculate the nq angles
		TVector3 norm_scatter = qVec.Cross( beamVec );
		norm_scatter 	= norm_scatter.Unit();

		TVector3 norm_reaction = qVec.Cross( nVec );
		norm_reaction 	= norm_reaction.Unit();

		phi_nq 	= norm_scatter.Angle( norm_reaction );
		theta_nq = nVec.Angle( qVec );
		CosTheta_nq = cos(theta_nq);
		TVector3 direction = norm_scatter.Cross(norm_reaction);
		if( direction.Z() > 0 ){ // this means the phi_rq should be between 0 and pi
		}
		else if( direction.Z() < 0 ){ // this means the phi_rq should be between -pi and 0
			phi_nq *= (-1);
		}
		nVec.Clear();

		nVec.SetMagThetaPhi(p_n,theta_n,phi_n);
		E_n 	= sqrt( mN*mN + p_n*p_n );
		double W_primeSq = mD*mD - Q2 + mN*mN + 2.*mD*(nu-E_n) - 2.*nu*E_n + 2.*q*p_n*cos(theta_nq);
		Wp = sqrt(W_primeSq);
		Xp = Q2/(2.*( nu*(mD-E_n) + p_n*q*CosTheta_nq));
		As = (E_n - p_n*CosTheta_nq)/mN;



		virt = ( pow(mD-E_n,2) - p_n*p_n - mP*mP )/(mP*mP);
		if( Wp > 1.8 && Q2 < 10 && Q2 > 2 && sqrt(W2)>2 && CosTheta_nq < -0.8 ) {
			int Nu_bin = (int) ((virt+0.9)/0.2);
			if( Nu_bin >= 0 || Nu_bin < 4 ){
				hs_xpnu[Nu_bin]->Fill( Xp , weight );
				if( Xp > 0.6 && Xp < 1 ){
					h_nuhi->Fill( virt , weight );
				}
				if( fabs(Xp-0.3)<0.05) {
					h_nulo->Fill( virt , weight );
				}
			}
			h_xpfull->Fill( Xp , weight );

			int As_bin = (int) ((As-1.2)/0.1);
			if( As_bin < 0 || As_bin > 4 ) continue;
			hs_xp[As_bin]->Fill( Xp , weight );


			if( Xp > 0.6 && Xp < 1 ){
				h_ashi->Fill( As , weight );
			}


			if( fabs(Xp-0.3)<0.05) {
				h_pn->Fill( p_n , weight );
				h_aslo->Fill( As , weight );
			}
		}

	}

	return;
}

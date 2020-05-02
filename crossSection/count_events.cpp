#include <cmath>
#include <iostream>
#include <fstream>

#include "../../deuteron_dis/include/constants.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TVectorT.h"
#include "TVector3.h"
#include "TF1.h"

using namespace std;

double GetAcceptance( TH2D * map, double xVal, double yVal );
double kinLimit( double *x , double *p);

int main( int argc, char ** argv){
	if( argc != 5 ){
		cerr << "Wrong number of arguments. Instead please use:\n\t./count_events [inputPhysicsRootFile] [acceptanceMap] [nEv gen in map per pt] [1 for Acc Corr]\n";
		return -1;
	}
	
	int n_gen = atoi(argv[3]);
	int accCorr = atoi(argv[4]);

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

	// Load input acceptance map
	TFile * mapFile = new TFile(argv[2]);
	TH2D ** HiMap = new TH2D*[50];
	TH2D ** HiCnt = new TH2D*[50];
	TH2D ** LoMap = new TH2D*[50];
	TH2D ** LoCnt = new TH2D*[50];
	for( int i = 0 ; i < 50; i++){
		HiMap[i] = (TH2D*)mapFile->Get(Form("HiMap_%i",i));
		HiCnt[i] = (TH2D*)mapFile->Get(Form("HiCnt_%i",i));
		LoMap[i] = (TH2D*)mapFile->Get(Form("LoMap_%i",i));
		LoCnt[i] = (TH2D*)mapFile->Get(Form("LoCnt_%i",i));
	}	

	// Storage for my counts:
	double counts[2][1][1][50] = {0.};
	double errs[2][1][1][50] = {0.};

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
		if( event % 100000 == 0 ) cerr << "Working on event " << event << " (" << int((double)event/nEv*100) << "%" << " done)" << "\n";


		// Grab corresponding momentum bin for acceptance map
		int momBin = (int)((pN_mag - 0.20)/0.01);

		// Kinematic cut on Wp and Q2 min
		if( Wp < 2. ) continue;
		if( Q2 < 2. ) continue;
		if( momBin < 0 || momBin > 49) continue; // if it's outside the map momentum (i.e. < 0.2 or > 0.7)

		// Basic fiducial cut on CosTheta_qn
		if( CosTheta_qn > -0.9 ) continue; // common thetq_qn region to both hi and lo Xp

		int bin_xp, bin_q2, bin_cthqn, bin_pn;
		double acceptance = 0.;

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

		// Low Xp Region:
		if( fabs(Xp - 0.3) < 0.005 ){
			bin_xp = 0;

			// Create the kinematic line and make sure we are less than that
			TF1 * model = new TF1("edge",kinLimit,0,10,3);
			model->SetParameters(0.3,Ebeam-2.2,(momBin+1)*0.01+0.2);	//nu max 8.6 because 10.6 beam energy and generate |p_e|>2
			
			if( Q2 < model->GetX( CosTheta_qn ) ){
				acceptance = GetAcceptance( LoMap[momBin], Q2, CosTheta_qn );
				//if( accCorr ) acceptance = GetAcceptance( LoMap[momBin], Q2, CosTheta_qn );
				//else acceptance = 1;
			}
			delete model;
		}
		// High Xp Region:
		else if( fabs(Xp - 0.505) < 0.005 ){
			bin_xp = 1;

			//// Create the kinematic line and make sure we are less than that
			TF1 * model = new TF1("edge",kinLimit,0,10,3);
			model->SetParameters(0.5,Ebeam-2.2,(momBin+1)*0.01+0.2);	//nu max 8.6 because 10.6 beam energy and generate |p_e|>2

			if( Q2 < model->GetX( CosTheta_qn ) ){
				acceptance = GetAcceptance( HiMap[momBin], Q2, CosTheta_qn );
				//if( accCorr) acceptance = GetAcceptance( HiMap[momBin], Q2, CosTheta_qn );
				//else acceptance = 1;
			}
			delete model;
		}
			

		// TODO: need to fix up map so I don't have a bunch of nan's if they aren't filled
		//if( acceptance != acceptance || acceptance == 0){
		//	cerr << "Acceptance error [A, X', Q2, CosThQN, Pn]: " << acceptance << " " << Xp << " " <<  Q2 << " " << CosTheta_qn << " " << momBin << "\n";
		//	continue;
		//}

		//acceptance = 1;	
		if( acceptance > 0 && acceptance==acceptance ){
			bin_q2		= 0;
			bin_cthqn 	= 0;
			bin_pn		= momBin;
			if( !accCorr ) acceptance = 1;
			//cout << theta_e << " " << Ebeam - nu << " " << theta_n << " " << pN_mag << " " << phi_en << "\n";
			//cout << Xp << " " << Q2 << " " << CosTheta_qn << "\n\n";
			counts[bin_xp][bin_q2][bin_cthqn][bin_pn] = counts[bin_xp][bin_q2][bin_cthqn][bin_pn] + 1./acceptance;
			errs[bin_xp][bin_q2][bin_cthqn][bin_pn] = errs[bin_xp][bin_q2][bin_cthqn][bin_pn] + pow(1./acceptance,4)*(acceptance/n_gen + pow(acceptance,2)/n_gen );
		}

	}
	
	ofstream output_lox;
	ofstream output_hix;
	if( accCorr){
		output_lox.open("accCorr-lox.txt");
		output_hix.open("accCorr-hix.txt");
	}	
	else{
		output_lox.open("true-lox.txt");
		output_hix.open("true-hix.txt");
	}
			

	for( int pBin = 0; pBin < 50; pBin++) output_lox << pBin*0.01+0.200 << " " << counts[0][0][0][pBin] << " " << errs[0][0][0][pBin] << "\n";
	for( int pBin = 0; pBin < 50; pBin++) output_hix << pBin*0.01+0.200 << " " << counts[1][0][0][pBin] << " " << errs[1][0][0][pBin] << "\n";

	output_lox.close();
	output_hix.close();
	
	return 1;
}

double GetAcceptance( TH2D * map, double xVal, double yVal ){
	return map->GetBinContent( map->GetXaxis()->FindBin(xVal), map->GetYaxis()->FindBin(yVal) );
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

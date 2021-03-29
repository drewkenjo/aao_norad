#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TF1.h>
#include <time.h>
#include <TSystem.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TStyle.h>

using namespace std;

//Double_t EBeam = 10.6;
Double_t EBeam = 7.535;
const Double_t ProtonMass = 0.93827;
const Double_t Pi0mass = 0.134976;
const Double_t Etamass = 0.5473;
const Double_t me2 = TMath::Power(0.000511,2.);
const Double_t alpha_em = 1./137.036;

TFile *outfile;
TTree *outtree;

TRandom3 *ThisMCrandom;
TLorentzVector v4Beam, v4Target, v4Elec, v4Prot, v4Phot1, v4Phot2, v4Pi0, v4Q;
Double_t gen_Q2, gen_xB, gen_t, gen_phi, gen_phie;
Double_t gen_elec_p, gen_elec_t, gen_elec_f;
Double_t gen_prot_p, gen_prot_t, gen_prot_f;
Double_t gen_phot1_E, gen_phot1_theta, gen_phot1_phi, gen_phot1_el_ang;
Double_t gen_phot2_E, gen_phot2_theta, gen_phot2_phi, gen_phot2_el_ang;
Double_t gen_pi0_p, gen_pi0_t, gen_pi0_f, gen_gg_angle;
Double_t gen_pi0_CMt, gen_pi0_CMf;
Double_t pi0_xsec;
Double_t gen_vx, gen_vy, gen_vz;
Double_t rad_dE1, rad_dE2;

void Kinematics(Double_t x, Double_t q, Double_t t, Double_t phi);
Bool_t GenAcc(Double_t x, Double_t q, Double_t t);
void MakeTree(void);
void FillTree(void);

Double_t GetMinQ2(Double_t Xbj);
Double_t GetMaxQ2(Double_t Xbj);
Double_t GetQ2_thetaX(Double_t Xbj, Double_t theta);
Double_t GetQ2_WX(Double_t Xbj, Double_t W);
Double_t GetQ2_EpX(Double_t Xbj, Double_t Ep);

Bool_t AboveTMin(Double_t DVCS_t_Pr, Double_t DVCS_Xbj, Double_t DVCS_Q2);
double dvmpx(double t, double xb, double Q2, double PHI_G , double E, double heli);

int main(Int_t argc, Char_t *argv[]){
cout << "---------------------------------------------------" << endl;
cout << "|               pi0 CLAS12 generator               |" << endl;
cout << "---------------------------------------------------" << endl;

gROOT->SetStyle("Plain");
gStyle->SetOptTitle(0);
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
gStyle->SetOptFit(0);
gStyle->SetOptDate(0);

int expo = 3;
if(argc>1)expo=atoi(argv[1]);
if(argc>2)EBeam=atof(argv[2]);

cout << "!!!!!!!!!!!!!! EBeam = ";cout<<Form("%1.2f",EBeam);cout<<" GeV !!!!!!!!!!!!!" << endl;

ThisMCrandom = new TRandom3(0);
for(int MCinit=0;MCinit<100;MCinit++)ThisMCrandom->Uniform();
cout << "My seed " << ThisMCrandom->Uniform() << endl;

Double_t minxB, minQ2, mint, minphi;
Double_t maxxB, maxQ2, maxt, maxphi;

minxB = 0.05;
maxxB = 0.8;
minQ2 = 0.5;
maxQ2 = 14.;
mint = 0.01;
maxt = 5.;

minxB = 0.3;
maxxB = 0.4;
minQ2 = 3;
maxQ2 = 4.;
mint = 0.3;
maxt = 0.8;

minphi = 0.5*TMath::DegToRad();
maxphi = 359.5*TMath::DegToRad();

Long_t NeventsTot = (Long_t)TMath::FloorNint(TMath::Power(10.,expo));

Double_t thisX, thisQ, thisT, thisPhi;

MakeTree();

ofstream outlund;
outlund.open("rename_pi0_lund.txt");   
outlund.precision(4);
outlund << right;
for(int ie=0;ie<NeventsTot;){
      thisX = minxB+(maxxB-minxB)*ThisMCrandom->Uniform();
      thisQ = minQ2+(maxQ2-minQ2)*ThisMCrandom->Uniform();
      thisT = mint+(maxt-mint)*ThisMCrandom->Uniform();
      thisPhi = minphi+(maxphi-minphi)*ThisMCrandom->Uniform();
      if(GenAcc(thisX,thisQ,thisT)){
	pi0_xsec = dvmpx(-thisT,thisX,thisQ,thisPhi,EBeam,1);
	if( pi0_xsec > 5 * ThisMCrandom->Uniform() ){
		Kinematics(thisX,thisQ,thisT,thisPhi);
		FillTree();
		ie++;
                 outlund << "4 1 1 0 -1 11 " << EBeam << " 2212 3 " << pi0_xsec << endl;
                 outlund << "1 -1 1   11 0 0 " 
                         << std::setw(8) << v4Elec.Px() << " " << setw(8) << v4Elec.Py() << " " << setw(8) << v4Elec.Pz() << " " << setw(8) << v4Elec.E() << " 0.0005 " << gen_vx << " " << gen_vy << " " << gen_vz << endl;
                 outlund << "2  1 1 2212 0 0 " 
                         << setw(8) << v4Prot.Px() << " " << setw(8) << v4Prot.Py() << " " << setw(8) << v4Prot.Pz() << " " << setw(8) << v4Prot.E() << " 0.9383 " << gen_vx << " " << gen_vy << " " << gen_vz << endl;
                 outlund << "3  0 1   22 0 0 " 
                         << setw(8) << v4Phot1.Px() << " " << setw(8) << v4Phot1.Py() << " " << setw(8) << v4Phot1.Pz() << " " << setw(8) << v4Phot1.E() << " 0.0000 " << gen_vx << " " << gen_vy << " " << gen_vz << endl;
                 outlund << "4  0 1   22 0 0 " 
                         << setw(8) << v4Phot2.Px() << " " << setw(8) << v4Phot2.Py() << " " << setw(8) << v4Phot2.Pz() << " " << setw(8) << v4Phot2.E() << " 0.0000 " << gen_vx << " " << gen_vy << " " << gen_vz << endl;

	}
      }
 }
 outtree->Write();
 outfile->Close();
 outlund.close();
 return 0;
}

void MakeTree(void){
 outfile = new TFile("pi0.root","recreate");
 outtree = new TTree("pi0","pi0");
 outtree->Branch("gen_Q2",&gen_Q2,"gen_Q2/D");
 outtree->Branch("gen_xB",&gen_xB,"gen_xB/D");
 outtree->Branch("gen_t",&gen_t,"gen_t/D");
 outtree->Branch("gen_phi",&gen_phi,"gen_phi/D");
 outtree->Branch("gen_phie",&gen_phie,"gen_phie/D");
 outtree->Branch("gen_elec_p",&gen_elec_p,"gen_elec_p/D");
 outtree->Branch("gen_elec_t",&gen_elec_t,"gen_elec_t/D");
 outtree->Branch("gen_elec_f",&gen_elec_f,"gen_elec_f/D");
 outtree->Branch("gen_prot_p",&gen_prot_p,"gen_prot_p/D");
 outtree->Branch("gen_prot_t",&gen_prot_t,"gen_prot_t/D");
 outtree->Branch("gen_prot_f",&gen_prot_f,"gen_prot_f/D");
 outtree->Branch("gen_phot1_E",&gen_phot1_E,"gen_phot1_E/D");
 outtree->Branch("gen_phot1_theta",&gen_phot1_theta,"gen_phot1_theta/D");
 outtree->Branch("gen_phot1_phi",&gen_phot1_phi,"gen_phot1_phi/D");
 outtree->Branch("gen_phot2_E",&gen_phot2_E,"gen_phot2_E/D");
 outtree->Branch("gen_phot2_theta",&gen_phot2_theta,"gen_phot2_theta/D");
 outtree->Branch("gen_phot2_phi",&gen_phot2_phi,"gen_phot2_phi/D");
 outtree->Branch("gen_pi0_p",&gen_pi0_p,"gen_pi0_p/D");
 outtree->Branch("gen_pi0_t",&gen_pi0_t,"gen_pi0_t/D");
 outtree->Branch("gen_pi0_f",&gen_pi0_f,"gen_pi0_f/D");
 outtree->Branch("gen_gg_angle",&gen_gg_angle,"gen_gg_angle/D");
 outtree->Branch("gen_pi0_CMt",&gen_pi0_CMt,"gen_pi0_CMt/D");
 outtree->Branch("gen_pi0_CMf",&gen_pi0_CMf,"gen_pi0_CMf/D");
 outtree->Branch("gen_vx",&gen_vx,"gen_vx/D");
 outtree->Branch("gen_vy",&gen_vy,"gen_vy/D");
 outtree->Branch("gen_vz",&gen_vz,"gen_vz/D");
 outtree->Branch("pi0_xsec",&pi0_xsec,"pi0_xsec/D");
}

void FillTree(void){
 outtree->Fill();
}

void Kinematics(Double_t x, Double_t q, Double_t t, Double_t phi){
 //kinematic variables in the lab
 //Double_t lnq2me2 = TMath::Log(q/me2);
 //Double_t delta_s = alpha_em/TMath::Pi()*(lnq2me2-1.);
 rad_dE1 = 0;//0.00001 + EBeam * TMath::Power(ThisMCrandom->Uniform(),1./delta_s);	
 rad_dE2 = 0;//0.00001 + ( EBeam - q/(2*ProtonMass*x) -rad_dE1 ) * TMath::Power(ThisMCrandom->Uniform(),1./delta_s);

 v4Beam.SetXYZT(0,0,EBeam-rad_dE1,TMath::Sqrt( (EBeam-rad_dE1)*(EBeam-rad_dE1)+me2));
 v4Target.SetXYZT(0.,0.,0.,ProtonMass);
 Double_t Px, Py, Pz, P, E, Theta, Phi;
 E = EBeam-rad_dE1 - q/(2*ProtonMass*x)-rad_dE2;
 P = TMath::Sqrt(E*E-me2);
 Theta = 2.*TMath::ASin( TMath::Sqrt(q/(4*EBeam*E)) );
 Phi = ThisMCrandom->Uniform()*2.*TMath::Pi();
 //cout << "Setting electron phi=0" << endl;
 //Phi = 0;
 Px = P*TMath::Sin(Theta)*TMath::Cos(Phi);
 Py = P*TMath::Sin(Theta)*TMath::Sin(Phi);
 Pz = P*TMath::Cos(Theta);
 v4Elec.SetXYZT(Px,Py,Pz,E);
 v4Q = v4Beam-v4Elec;

 gen_Q2 = q;
 gen_xB = x;
 gen_t = t;
 gen_phi = TMath::RadToDeg()*phi;
 gen_phie = TMath::RadToDeg()*Phi;
 float W2 = ProtonMass*ProtonMass+gen_Q2*(1./gen_xB-1.);
 float E1CM = (W2+gen_Q2+ProtonMass*ProtonMass)/(2*TMath::Sqrt(W2));
 float P1CM = TMath::Sqrt(E1CM*E1CM-ProtonMass*ProtonMass);
 float E3CM = (W2+ProtonMass*ProtonMass-Pi0mass*Pi0mass)/(2*TMath::Sqrt(W2));
 float P3CM = TMath::Sqrt(E3CM*E3CM-ProtonMass*ProtonMass);
 float TMIN = TMath::Power(gen_Q2+Pi0mass*Pi0mass,2)/(4.*W2)-TMath::Power(P1CM-P3CM,2);
 float thetacm = 2.*TMath::ASin(TMath::Sqrt( (gen_t-TMath::Abs(TMIN))/(4*P1CM*P3CM) ));
 TLorentzVector v4ProtCM, v4Pi0CM, v4Phot1CM, v4Phot2CM, v4QCM, v4TargetCM ;
 v4QCM = v4Q;
 v4TargetCM = v4Target;
 TVector3 BETA, CMboost;
 BETA = (v4Q+v4Target).BoostVector();
 TVector3 LABboost = BETA;
 CMboost = - BETA;
 v4QCM.Boost(CMboost);
 Double_t QThetaCM = (v4QCM.Vect()).Theta();
 Double_t QPhiCM = (v4QCM.Vect()).Phi();
 v4TargetCM.Boost(CMboost);
 v4TargetCM.RotateZ(-QPhiCM);
 v4QCM.RotateZ(-QPhiCM);
 v4TargetCM.RotateY(-QThetaCM);
 v4QCM.RotateY(-QThetaCM);

 E = TMath::Sqrt(Pi0mass*Pi0mass+P3CM*P3CM);
 P = P3CM;
 Px = P*TMath::Sin(thetacm)*TMath::Cos(phi+TMath::Pi());
 Py = P*TMath::Sin(thetacm)*TMath::Sin(phi+TMath::Pi());
 Pz = P*TMath::Cos(thetacm);
 v4Pi0CM.SetXYZT(Px,Py,Pz,E);

 thetacm = TMath::Pi()-thetacm;
 E = TMath::Sqrt(P*P+ProtonMass*ProtonMass);
 Px = P*TMath::Sin(thetacm)*TMath::Cos(phi);
 Py = P*TMath::Sin(thetacm)*TMath::Sin(phi);
 Pz = P*TMath::Cos(thetacm);
 v4ProtCM.SetXYZT(Px,Py,Pz,E);

 v4ProtCM.RotateY(QThetaCM);
 v4Pi0CM.RotateY(QThetaCM);
 v4ProtCM.RotateZ(QPhiCM);
 v4Pi0CM.RotateZ(QPhiCM);
 v4ProtCM.Boost(LABboost);
 v4Pi0CM.Boost(LABboost);
 v4Prot = v4ProtCM;
 v4Pi0 = v4Pi0CM;

 //generate photons here
 BETA = v4Pi0.BoostVector();
 Theta = 2*ThisMCrandom->Uniform()-1.;
 Phi = 2*TMath::Pi()*ThisMCrandom->Uniform();
 E = Pi0mass/2.;
 P = E;
 Px = P*TMath::Sqrt(1-Theta*Theta)*TMath::Cos(Phi);
 Py = P*TMath::Sqrt(1-Theta*Theta)*TMath::Sin(Phi);
 Pz = P*Theta;
 v4Phot1.SetXYZT(Px,Py,Pz,E);
 v4Phot2.SetXYZT(-Px,-Py,-Pz,E);
 v4Phot1.Boost(BETA);
 v4Phot2.Boost(BETA);

 gen_elec_p = v4Elec.P();
 gen_elec_t = v4Elec.Theta()*TMath::RadToDeg();
 gen_elec_f = v4Elec.Phi()*TMath::RadToDeg();
 gen_prot_p = v4Prot.P();
 gen_prot_t = v4Prot.Theta()*TMath::RadToDeg();
 gen_prot_f = v4Prot.Phi()*TMath::RadToDeg();
 gen_phot1_E = v4Phot1.E();
 gen_phot1_theta = TMath::RadToDeg()*v4Phot1.Theta();
 gen_phot1_phi = TMath::RadToDeg()*v4Phot1.Phi();
 gen_phot1_el_ang = TMath::RadToDeg()*(v4Phot1.Vect()).Angle(v4Elec.Vect());
 gen_phot2_E = v4Phot2.E();
 gen_phot2_theta = TMath::RadToDeg()*v4Phot2.Theta();
 gen_phot2_phi = TMath::RadToDeg()*v4Phot2.Phi();
 gen_phot2_el_ang = TMath::RadToDeg()*(v4Phot2.Vect()).Angle(v4Elec.Vect());
 gen_pi0_p = v4Pi0.P();
 gen_pi0_t = TMath::RadToDeg()*v4Pi0.Theta();
 gen_pi0_f = TMath::RadToDeg()*v4Pi0.Phi();
 gen_gg_angle = TMath::RadToDeg() * (v4Phot1.Vect()).Angle(v4Phot2.Vect()) ;
 gen_pi0_CMt = TMath::RadToDeg()*TMath::ACos(Theta);
 gen_pi0_CMf = TMath::RadToDeg()*Phi;

 gen_vx =  0.0;//0.04*ThisMCrandom->Gaus();
 gen_vy =  0.0;//0.04*ThisMCrandom->Gaus();
 gen_vz = -2.5 + 5.0*ThisMCrandom->Uniform();

}

Bool_t GenAcc(Double_t x, Double_t q, Double_t t){
 Double_t Eprime = EBeam-q/(2*ProtonMass*x);
 Double_t El_Theta = 2*TMath::RadToDeg()*TMath::ASin(TMath::Sqrt(q/(4*EBeam*Eprime) ));
 Double_t W = TMath::Sqrt(ProtonMass*ProtonMass+q*(1./x-1.));
 if( q>0.5 && x>0.05 && El_Theta>4. && W>1.7 && Eprime>0.5 && AboveTMin(t,x,q) )return kTRUE;
 return kFALSE;
}

Bool_t AboveTMin(Double_t DVCS_t_Pr, Double_t DVCS_Xbj, Double_t DVCS_Q2){
 Double_t DVCS_W2 = ProtonMass*ProtonMass+DVCS_Q2*(1./DVCS_Xbj-1.); 
 Double_t E1CM = (DVCS_W2+DVCS_Q2+ProtonMass*ProtonMass)/(2*TMath::Sqrt(DVCS_W2));
 Double_t P1CM = TMath::Sqrt(E1CM*E1CM-ProtonMass*ProtonMass);
 Double_t E3CM = (DVCS_W2-Pi0mass*Pi0mass+ProtonMass*ProtonMass)/(2*TMath::Sqrt(DVCS_W2));
 Double_t P3CM = TMath::Sqrt(E3CM*E3CM-ProtonMass*ProtonMass);
 Double_t TMIN = TMath::Power(DVCS_Q2+Pi0mass*Pi0mass,2)/(4.*DVCS_W2)-TMath::Power(P1CM-P3CM,2);
 Double_t TMAX = TMath::Power(DVCS_Q2+Pi0mass*Pi0mass,2)/(4.*DVCS_W2)-TMath::Power(P1CM+P3CM,2);
 return (TMath::Abs(DVCS_t_Pr)>TMath::Abs(TMIN)&&TMath::Abs(DVCS_t_Pr)<TMath::Abs(TMAX));
}

Double_t GetTmin(Double_t DVCS_Xbj, Double_t DVCS_Q2){
 Double_t DVCS_W2 = ProtonMass*ProtonMass+DVCS_Q2*(1./DVCS_Xbj-1.); 
 Double_t E1CM = (DVCS_W2+DVCS_Q2+ProtonMass*ProtonMass)/(2*TMath::Sqrt(DVCS_W2));
 Double_t P1CM = TMath::Sqrt(E1CM*E1CM-ProtonMass*ProtonMass);
 Double_t E3CM = (DVCS_W2-Pi0mass*Pi0mass+ProtonMass*ProtonMass)/(2*TMath::Sqrt(DVCS_W2));
 Double_t P3CM = TMath::Sqrt(E3CM*E3CM-ProtonMass*ProtonMass);
 Double_t TMIN = TMath::Power(DVCS_Q2+Pi0mass*Pi0mass,2)/(4.*DVCS_W2)-TMath::Power(P1CM-P3CM,2);
 return TMIN;
}

Double_t GetTmax(Double_t DVCS_Xbj, Double_t DVCS_Q2){
 Double_t DVCS_W2 = ProtonMass*ProtonMass+DVCS_Q2*(1./DVCS_Xbj-1.);
 Double_t E1CM = (DVCS_W2+DVCS_Q2+ProtonMass*ProtonMass)/(2*TMath::Sqrt(DVCS_W2));
 Double_t P1CM = TMath::Sqrt(E1CM*E1CM-ProtonMass*ProtonMass);
 Double_t E3CM = (DVCS_W2-Pi0mass*Pi0mass+ProtonMass*ProtonMass)/(2*TMath::Sqrt(DVCS_W2));
 Double_t P3CM = TMath::Sqrt(E3CM*E3CM-ProtonMass*ProtonMass);
 Double_t TMAX = TMath::Power(DVCS_Q2+Pi0mass*Pi0mass,2)/(4.*DVCS_W2)-TMath::Power(P1CM+P3CM,2);
 return TMAX;
}

Double_t GetT(Double_t DVCS_Xbj, Double_t DVCS_Q2, Double_t theta_CM){
 Double_t DVCS_W2 = ProtonMass*ProtonMass+DVCS_Q2*(1./DVCS_Xbj-1.);
 Double_t E1CM = (DVCS_W2+DVCS_Q2+ProtonMass*ProtonMass)/(2*TMath::Sqrt(DVCS_W2));
 Double_t P1CM = TMath::Sqrt(E1CM*E1CM-ProtonMass*ProtonMass);
 Double_t E3CM = (DVCS_W2-Pi0mass*Pi0mass+ProtonMass*ProtonMass)/(2*TMath::Sqrt(DVCS_W2));
 Double_t P3CM = TMath::Sqrt(E3CM*E3CM-ProtonMass*ProtonMass);
 return GetTmin(DVCS_Xbj,DVCS_Q2)-4.*P1CM*P3CM*TMath::Power(TMath::Sin(TMath::DegToRad()*theta_CM/2.),2.);
}

Double_t Fun_T(Double_t *X, Double_t *PAR){
 return GetT(PAR[0],PAR[1],X[0]);
}

Double_t GetMinQ2(Double_t Xbj){
 Double_t result = 1.;
 Double_t thetalim = GetQ2_thetaX(Xbj,21.);
 Double_t Wlim = GetQ2_WX(Xbj,2.);
 //cout << "thetalim = " << thetalim << endl;
 //cout << "Wlim = " << Wlim << endl;
 if(TMath::Max(thetalim,Wlim)>1.)result = TMath::Max(thetalim,Wlim);
 return result;
}

Double_t GetMaxQ2(Double_t Xbj){
 Double_t result = 1.;
 Double_t eplim = GetQ2_EpX(Xbj,0.8);
 Double_t thetalim = GetQ2_thetaX(Xbj,45.);
 //cout << "thetalim = " << thetalim << endl;
 //cout << "elim = " << eplim << endl;
 if(TMath::Min(eplim,thetalim)>1.)result = TMath::Min(eplim,thetalim);
 return result;
}

Double_t GetQ2_thetaX(Double_t Xbj, Double_t theta){
 Double_t THETA  = 2*EBeam*TMath::Power(TMath::Sin(TMath::DegToRad()*theta/2.),2.);
 return 2.*EBeam*THETA*ProtonMass*Xbj/(THETA+ProtonMass*Xbj);
}

Double_t GetQ2_WX(Double_t Xbj, Double_t W){
 return (W*W-ProtonMass*ProtonMass)/(1./Xbj-1.);
}

Double_t GetQ2_EpX(Double_t Xbj, Double_t Ep){
 return 2.*ProtonMass*(EBeam-Ep)*Xbj;
}

/**********************************************************************************/
/*                             pi0/eta 2019 Model begins here                              */
/**********************************************************************************/
/* https://github.com/vkubarovsky/Model-for-the-pi0-eta-exclusive-cross-section   */
/*

     REAL FUNCTION DVMPX(del2,xb,Q2, Phi_g,E,heli,MESONMASS)
C
c  dsigma/dQ2 dX dt dphi for ep-->ep pi0/eta
C
C exc pi0/eta x-section
c
c input:
c del2=t (negative GeV^2)           NEGATIVE !!!!
c xb,Q2 x and Q^2 (GeV^2)
c Phi_g angle in the photon frame (radians)
c E energy of the electron in GeV
c heli electron helicity -1.0 or +1.0
c MESONMASS is the mass of the pi0 or eta

PARAMETERS pi0/eta

       DATA P/
     *   6.2624,  0.8032, -0.1152,  1.7431,  0.0,
     *  92.9727,  3.4071, -1.9099, -0.1199,  0.0, 16.8603,  2.1722,
     *   6.9439,  1.7523, -1.2875,  0.6822,  0.0,
     *  17.0423,  1.1264,  0.0491,  1.6594,  0.0, 21.6379,  3.9873/
*/

double dvmpx_phi(double *X, double *P){
 return dvmpx(P[2],P[0],P[1],X[0]*TMath::DegToRad(),P[3],1);
}

double dvmpx(double t, double xb, double Q2, double PHI_G , double E, double heli){
 //  dsigma/dQ2 dX dt dphi for ep-->ep pi0

  /* pi0 */  	
double p[12]={6.2059, 0.8019, -0.1064, 1.8361, 0.00, 94.6402,  3.4433, -1.9773, -0.1315,  0.0000, 16.9977,  2.1240};

/* eta 
double p[12]={6.8921,  1.7848, -1.3187,  0.7766,  0.00,17.6413,  1.1400,  0.0135,  1.6014,  0.0000, 21.9200,  4.0043} ;
*/
 
 if(xb<Q2/(2.0*ProtonMass*E)||xb>1.0)return 0.;
 double nu  = Q2/(2.*ProtonMass*xb);
 double y = nu/E;
 double e1 = TMath::Power(y*xb*ProtonMass,2)/Q2;
 double EPS = (1.0-y-e1)/(1-y+y*y/2+e1);
 if(EPS<0.||EPS>1.)return 0.;
 double W2  = ProtonMass*ProtonMass + 2.0*ProtonMass*nu - Q2;
 if(W2<TMath::Power(ProtonMass+Pi0mass,2.))return 0.;
 double W = TMath::Sqrt(W2);
 double E1cm = ProtonMass*(ProtonMass + nu)/W;
 double P1cm = ProtonMass*TMath::Sqrt(nu*nu+Q2)/W;
 double E2cm = (W2 + ProtonMass*ProtonMass-Pi0mass*Pi0mass)/(2.*W);
 if(E2cm<ProtonMass)return 0.;
 double P2cm = TMath::Sqrt(E2cm*E2cm-ProtonMass*ProtonMass);

 double E3cm = (W2-Pi0mass*Pi0mass+ProtonMass*ProtonMass)/(2*W);
 double P3cm = TMath::Sqrt(E3cm*E3cm-ProtonMass*ProtonMass);

 double tmax = 2.0*(ProtonMass*ProtonMass - E1cm*E2cm + P1cm*P2cm);
 double tmin = 2.0*(ProtonMass*ProtonMass - E1cm*E2cm - P1cm*P2cm);
 if(t<tmin||t>tmax)return 0.;
 double FLUXW = alpha_em/(2*TMath::Pi()) * y*y/(1-EPS)*(1-xb)/xb/Q2;

 tmin = -(TMath::Power(Q2+Pi0mass*Pi0mass,2)/4./W2-TMath::Power(P1cm-P3cm,2));
 double SLOPE = 2.*1.1*TMath::Log(xb);

 double T  = TMath::Abs(t);
 double T0 = TMath::Abs(tmin);
// cout << "T-T0=" << T-T0 << " , T=" << T << " , T0=" << T0 << endl;

 double HT = p[0]*TMath::Exp(-(p[1]+p[2]*(TMath::Log(xb)-TMath::Log(0.15)))*T)*TMath::Power(Q2,p[3]/2.);
 double ET = p[5]*TMath::Exp(-(p[6]+p[7]*(TMath::Log(xb)-TMath::Log(0.15)))*T)*TMath::Power(Q2,p[8]/2.);
 double HTEBAR = p[10]*TMath::Exp(-p[11]*T);
 
 double pi = TMath::Pi();  
 double hc2= 389379.36;
 double ProtonMass2 = ProtonMass*ProtonMass;
 double ksi = xb/(2-xb)*(1.+ProtonMass2/Q2);
 double phase = 16.*pi*(W2-ProtonMass2)*TMath::Sqrt(W2*W2+Q2*Q2+ProtonMass2*ProtonMass2+2.*W2*Q2-2.*W2*ProtonMass2+2.*Q2*ProtonMass2);

 double S_T  = hc2*4.*pi*alpha_em/(2.*phase*Q2*Q2)*((1.-ksi*ksi)*HT*HT+(T-T0)/(8.*ProtonMass2) * ET*ET);
 double S_L  = 0.0;
 double S_LT = hc2*4.*pi*alpha_em/(TMath::Sqrt(2.)*phase*TMath::Power(Q2,1.5)) * ksi*TMath::Sqrt(1.-ksi*ksi)*TMath::Sqrt(T-T0)/(2.*ProtonMass)*HTEBAR*HTEBAR;
 double S_TT = -hc2*4.*pi*alpha_em/(2.*phase*Q2*Q2)*(T-T0)/(8.*ProtonMass2) * ET*ET;
 double S_LTP = 0.;

 /*
 double S_T =  (814.59900+ 0.0*TMath::Sqrt(T-T0))*TMath::Exp(SLOPE*T*0.44703799)*    1./TMath::Power(Q2+1.0876700,1.6934299);
 double S_L = Q2*(1404.0500 + 0.0*(T-T0))*TMath::Exp(SLOPE*T*0.69298601)*     1./TMath::Power(Q2+1.0876700,1.6934299);
 double S_LT =         608.48901*TMath::Sqrt(T-T0)*TMath::Exp(SLOPE*T*1.0290900)*    1./TMath::Power(Q2+1.0876700,1.6934299);
 double S_TT =              -5205.8999*(T-T0)*TMath::Exp(SLOPE*T*0.98690498)* 1./TMath::Power(Q2+1.0876700,1.6934299);
 double S_LTP = 0.;
 */
// cout << FLUXW << " S_T=" << S_T << " " << EPS << " " << S_L << " " << S_TT << " S_LT=" << S_LT << endl;


 double DVMPX = FLUXW/(2.*TMath::Pi())*( S_T + EPS*S_L + EPS * S_TT  * TMath::Cos(2*PHI_G) + TMath::Sqrt(2.*EPS*(1.+EPS))*S_LT * TMath::Cos(PHI_G) + heli*TMath::Sqrt(2.*EPS*(1.-EPS))*S_LTP * TMath::Sin(2*PHI_G) ) ;
      if(DVMPX<0.) DVMPX=0.;
 return DVMPX;
}


/**********************************************************************************/


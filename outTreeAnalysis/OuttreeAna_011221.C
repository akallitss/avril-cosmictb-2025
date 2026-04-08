#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

#ifdef __MAKECINT__
#pragma link C++ class vector<double>+;
#endif
//extrpolation of RD3 tracker to which det
const int DETEXTRA = 0;

const bool detinfit[] = {false,true,false,true};

//Analysis constant
const double ZDET[] = {686,//TPOT
		       209,//RD3 haut
		       252,//RD3 mid
		       303};//RD3 low

const double XDET[] = {1.70886e+02-1.07499e+00,//TPOT
		       3.98319299999999998e+01,//RD3 haut
		       3.98319299999999998e+01+6.15302e-02,//+0.140,//RD3 mid
		       3.98319299999999998e+01};//RD3 low

const double THETA[] = {0.8,//TPOT
		       1.1,//RD3 haut
		       1.1,//RD3 mid
		       1.1};//RD3 low

const double SCALE[] = {1.005,//TPOT
		       1.,//RD3 haut
		       1.,//RD3 mid
		       1.};//RD3 low

const double nsigma = 4.5;
const double pitch = 2.;
const double pitchrd3 = 0.5;
const int nrd3 = 4;


//for efficency
const double xdet = 2.78681e+02;
const double xres = 4.5;//mm

//Global for ray tree
double rayxup;
double rayyup;
double rayzup;
double rayxdn;
double rayydn;
double rayzdn;
double chi2x;
double chi2y;
vector <double> * fovrayx, * fovrayy;
vector <double> * fovclustmaxamp;
vector <double> * fovclustmaxstrip;
vector <double> * fovclustx;
vector <double> * fovextra;
vector <double> * fovextradet;
vector <int> * fovclustn;
vector <int> * fovclustsize;

vector <double> * clustamp,* clustmaxamp, *clustmaxstrip,  * clusttime,  * clustx;
vector <int>  * clustsize, *clustn;

double getXray(double z){
  //solve X=Az+B
  double Z_Up= rayzup;  double Z_Down= rayzdn;
  double X_Up= rayxup;  double X_Down= rayxdn;
  if(Z_Up-Z_Down == 0) {cout << "ERROR: Z_up == Z_Down" << endl; return -1;}  
  double A=(X_Up-X_Down)/(Z_Up-Z_Down);
  double B=X_Up-A*Z_Up;
  double x=A*z+B;//tadaaaa
  return x;
}
double getYray(double z){
  //solve Y=Az+B
  double Z_Up= rayzup;  double Z_Down= rayzdn;
  double Y_Up= rayyup;  double Y_Down= rayydn;
  if(Z_Up-Z_Down == 0) {cout << "ERROR: Z_up == Z_Down" << endl; return -1;}  
  double A=(Y_Up-Y_Down)/(Z_Up-Z_Down);
  double B=Y_Up-A*Z_Up;
  double y=A*z+B;//tadaaaa
  return y;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void OuttreeAna(TString filename, int maxev = -1)
{
  cout << "Opening " << filename << endl;
  TTree * treAna; TFile * fana;
  fana = new TFile(filename.Data());
  treAna =  (TTree*)fana->Get("Tout");
  int nev = treAna->GetEntries();
  if (maxev>0) nev=maxev;
  cout << "Tree = " << treAna->GetName() << "\t"<< treAna->GetEntries() << " events" << endl;

  rayxup=0; treAna->SetBranchAddress("RayXup", &rayxup);
  rayyup=0; treAna->SetBranchAddress("RayYup", &rayyup);
  rayxdn=0; treAna->SetBranchAddress("RayXdn", &rayxdn);
  rayydn=0; treAna->SetBranchAddress("RayYdn", &rayydn);
  rayzdn=0; treAna->SetBranchAddress("RayZdn", &rayzdn);
  rayzup=0; treAna->SetBranchAddress("RayZup", &rayzup);
  chi2x=0; treAna->SetBranchAddress("Chi2X", &chi2x);
  chi2y=0; treAna->SetBranchAddress("Chi2Y", &chi2y);
  
  clustamp=new vector <double>;treAna->SetBranchAddress("ClustAmp", &clustamp);
  clustmaxamp=new vector <double>;treAna->SetBranchAddress("ClustMaxAmp", &clustmaxamp);
  clustmaxstrip=new vector <double>;treAna->SetBranchAddress("ClustMaxStrip", &clustmaxstrip);
  clusttime=new vector <double>;treAna->SetBranchAddress("ClustTime", &clusttime);
  clustx=new vector <double>;treAna->SetBranchAddress("ClustX", &clustx);
  clustn=new vector <int>;treAna->SetBranchAddress("ClustN", &clustn);
  clustsize=new vector <int>;treAna->SetBranchAddress("ClustSize", &clustsize);

  //outtree
  TFile * ffinaltree = new TFile("finaltree.root","recreate");
  TTree * finaltree = new TTree("TF","Event");
  fovrayx=new vector <double>; finaltree->Branch("RayX", &fovrayx);
  fovrayy=new vector <double>; finaltree->Branch("RayY", &fovrayy);
  fovclustmaxamp=new vector <double>; finaltree->Branch("ClustMaxAmp", &fovclustmaxamp);
  fovclustmaxstrip=new vector <double>; finaltree->Branch("ClustMaxStrip", &fovclustmaxstrip);
  fovclustx=new vector <double>; finaltree->Branch("ClustX", &fovclustx);
  fovclustn=new vector <int>; finaltree->Branch("ClustN", &fovclustn);
  fovclustsize=new vector <int>; finaltree->Branch("ClustSize", &fovclustsize);
  fovextra=new vector <double>; finaltree->Branch("Extra", &fovextra);
  fovextradet=new vector <double>; finaltree->Branch("ExtraX", &fovextradet);

  /////////////////////////////////////////////////////////////////////////////////////////////////
  for(int iev=0; iev<nev;iev++){
    if(iev%1000==0) cout << iev << endl; 
    treAna->GetEntry(iev);
    fovrayx->clear();
    fovrayy->clear();
    fovclustmaxstrip->clear();
    fovclustx->clear();
    fovclustmaxamp->clear();
    fovclustn->clear();
    fovclustsize->clear();
    fovextra->clear();
    fovextradet->clear();

    double xrd3[nrd3], ard3[nrd3]; for(int id=0;id<nrd3;id++) {xrd3[id]=-1; ard3[id]=-1;} 
    if(clustn->size()>50) {cout << iev << " " << clustn->size() <<  endl;continue;}
    for(unsigned int icl=0;icl<clustn->size();icl++){//cluster loop
      int detn = clustn->at(icl); 
      if(detn<=0 || detn>4) {cout << detn << endl; continue;}
      
      double zdet = ZDET[detn-1];
      double rayx = getXray(zdet);
      double rayy = getYray(zdet);
      //rotation of tracker
      double theta = THETA[detn-1]*TMath::Pi()/180.;//in radian
      double xray = rayx*cos(theta)-rayy*sin(theta);
      double yray = rayx*sin(theta)+rayy*cos(theta);
      
      fovrayx->push_back(xray+XDET[detn-1]);//35000
      fovrayy->push_back(yray);//35000
      //cout << iev<< "] " << detn << "\t" << zdet << "\t" << xray << endl;
      fovclustx->push_back(clustx->at(icl)*SCALE[detn-1]);
      fovclustmaxamp->push_back(clustmaxamp->at(icl));
      fovclustmaxstrip->push_back(clustmaxstrip->at(icl));
      fovclustn->push_back(clustn->at(icl));
      fovclustsize->push_back(clustsize->at(icl));
      if(clustmaxamp->at(icl)>ard3[detn-1]) 
	{ard3[detn-1]=clustmaxamp->at(icl); xrd3[detn-1]=clustx->at(icl)*SCALE[detn-1]+XDET[detn-1];}
    }//cluster loop
    
    double Ex=0, Ez=0, Ex2=0, Exz=0, npt=0; bool trok=true;
    for(int id=0;id<nrd3;id++) 
      {
	if(!detinfit[id]) continue;
	//least sq method //Z = M.x + P
	if(ard3[id]<0) {trok=false;break;}
	//cout << "["<<iev<<"]["<<id<< "] " << xrd3[id]<< ", " <<ZDET[id] <<endl;
	npt++;
	Ez += xrd3[id];
	Ex += ZDET[id];
	Exz += xrd3[id]*ZDET[id];
	Ex2 += ZDET[id]*ZDET[id];
      } 
    if(trok)
      {
	
	double m=(npt*Exz-Ex*Ez)/(npt*Ex2-Ex*Ex);
	double p=(Ez-m*Ex)/npt;
	double xextra = m*ZDET[DETEXTRA]+p;
	fovextra->push_back(xextra);
	fovextradet->push_back(xrd3[DETEXTRA]);
	//cout << "["<<iev<<"] m=" << m << " ,p=" << p << " at z=" << ZDET[DETEXTRA] << "=" << xextra << " ~ " << xrd3[DETEXTRA] <<endl; 

      }
    finaltree->Fill();
  }//end for events
  finaltree->Write();

  //return;
  gStyle->SetOptFit(111);
  TCanvas * cres = new TCanvas("cres","cres",1400,600);
  cres->Divide(4,1);
  cres->cd(1);
  TH1F * hres1 = new TH1F("hres1", "hres1", 200,-20,20);
  finaltree->Draw("RayX-ClustX>>hres1","ClustN==1","");
  hres1->Fit("gaus");
  cres->cd(2);//TCanvas * c2 = new TCanvas();
  TH1F * hres2 = new TH1F("hres2", "hres2", 200,-20,20);
  finaltree->Draw("RayX-ClustX>>hres2","ClustN==2","");
  hres2->Fit("gaus");
  cres->cd(3);//TCanvas * c3 = new TCanvas();
  TH1F * hres3 = new TH1F("hres3", "hres3", 200,-20,20);
  finaltree->Draw("RayX-ClustX>>hres3","ClustN==3","");
  hres3->Fit("gaus");
  // TCanvas * c4 = new TCanvas();
  cres->cd(4);
  TH1F * hres4 = new TH1F("hres4", "hres4", 200,-20,20);
  finaltree->Draw("RayX-ClustX>>hres4","ClustN==4 && ClustMaxAmp>100","");
  hres4->Fit("gaus");
}

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
;// tpot, rd3_1, rd3_2, rd3_3

const bool detinfit[] = {false,true,true,true};// tpot, rd3_1, rd3_2, rd3_3

//Analysis constant
const double ZDET[] = {336,//TPOT upper deck is 686
		       297,//RD3 haut
		       252,//RD3 mid
		       207};//RD3 low

const double XDET[] = {73.3,//TPOT
		       294.55,//RD3 haut
		       294.6,//+0.140,//RD3 mid
		       294.764};//RD3 low

const double THETA[] = {1.,//TPOT
		       1.,//RD3 haut
		       1,//RD3 mid
		       1.};//RD3 low

const double SCALE[] = {1.,//TPOT
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
vector <double> * fovclustx1, * fovclustx2, * fovclustx3, * fovclustx4;
vector <double> * fovrdraym;
vector <int>  * clustsize, *clustn;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void OuttreeAna_NoRays(TString filename, int maxev = -1)
{
  cout << "Opening " << filename << endl;
  TTree * treAna; TFile * fana;
  fana = new TFile(filename.Data());
  treAna =  (TTree*)fana->Get("Tout");
  int nev = treAna->GetEntries();
  if (maxev>0) nev=maxev;
  cout << "Tree = " << treAna->GetName() << "\t"<< treAna->GetEntries() << " events" << endl;
  
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
  fovclustx1=new vector <double>; finaltree->Branch("ClustX1", &fovclustx1);
  fovclustx2=new vector <double>; finaltree->Branch("ClustX2", &fovclustx2);
  fovclustx3=new vector <double>; finaltree->Branch("ClustX3", &fovclustx3);
  fovclustx4=new vector <double>; finaltree->Branch("ClustX4", &fovclustx4);
  fovclustn=new vector <int>; finaltree->Branch("ClustN", &fovclustn);
  fovclustsize=new vector <int>; finaltree->Branch("ClustSize", &fovclustsize);
  fovextra=new vector <double>; finaltree->Branch("Extra", &fovextra);
  fovextradet=new vector <double>; finaltree->Branch("ExtraX", &fovextradet);
  fovrdraym=new vector <double>; finaltree->Branch("RDrayM", &fovrdraym);
  /////////////////////////////////////////////////////////////////////////////////////////////////
  for(int iev=0; iev<nev;iev++){
    if(iev%10000==0) cout << iev << endl; 
    treAna->GetEntry(iev);
    fovclustmaxstrip->clear();
    fovclustx->clear();
    fovclustx1->clear();    fovclustx2->clear();    fovclustx3->clear();    fovclustx4->clear();
    fovclustmaxamp->clear();
    fovclustn->clear();
    fovclustsize->clear();
    fovextra->clear();
    fovrdraym->clear();
    fovextradet->clear();

    double xrd3[nrd3],xrd3ms[nrd3], ard3[nrd3]; for(int id=0;id<nrd3;id++) {xrd3[id]=-1;xrd3ms[id]=-1; ard3[id]=-1;} 
    if(clustn->size()>50) {cout << iev << " " << clustn->size() <<  endl;continue;}
    for(unsigned int icl=0;icl<clustn->size();icl++){//cluster loop
      int detn = clustn->at(icl); 
      if(detn<=0 || detn>4) {cout << detn << endl; continue;}
      
      fovclustx->push_back(clustx->at(icl)*SCALE[detn-1]);
      switch(clustn->at(icl))
	{
	case 1 :       fovclustx1->push_back(clustx->at(icl)*SCALE[detn-1]); break;
	case 2 :       fovclustx2->push_back(clustx->at(icl)*SCALE[detn-1]); break;
	case 3 :       fovclustx3->push_back(clustx->at(icl)*SCALE[detn-1]); break;
	case 4 :       fovclustx4->push_back(clustx->at(icl)*SCALE[detn-1]); break;
	}
      fovclustmaxamp->push_back(clustmaxamp->at(icl));
      fovclustmaxstrip->push_back(clustmaxstrip->at(icl));
      fovclustn->push_back(clustn->at(icl));
      fovclustsize->push_back(clustsize->at(icl));
      if(clustmaxamp->at(icl)>ard3[detn-1]) 
	{ard3[detn-1]=clustmaxamp->at(icl); xrd3[detn-1]=clustx->at(icl)*SCALE[detn-1]+XDET[detn-1];
	  xrd3ms[detn-1]=clustmaxstrip->at(icl);}
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
	fovextradet->push_back(xrd3[DETEXTRA]);//xrd3[DETEXTRA]);
	fovrdraym->push_back(m);
	//cout << "["<<iev<<"] m=" << m << " ,p=" << p << " at z=" << ZDET[DETEXTRA] << "=" << xextra << " ~ " << xrd3[DETEXTRA] <<endl; 

      }
    finaltree->Fill();
  }//end for events
  finaltree->Write();

  return;
 
}
//TF->Draw("ExtraX-Extra:ExtraX>>(300,200,500,300,-5,5)","","colz");
// TF->Draw("ExtraX-Extra>>(300,-5,5)","","");
// TF->Draw("ExtraX-Extra:RDrayM>>(300,-1,1,300,-5,5)","","colz");
//

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLine.h"
#include "TLatex.h"
#include "TProfile.h"
#include "TMath.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

#ifdef __MAKECINT__
#pragma link C++ class vector<double>+;
#endif
//cosmetics
void SetStyle(Bool_t graypalette=false);

//efficiency stuff
int DETEFF = 0;//NPRHI or NZ
const bool PLOTEFF=true;
const bool SAVEEFF=false;
const double xmin = 20, xmax = 480, ymin= 20, ymax =230;//fiducial to measure efficiency
const double ROADEFF = 10;//mm
const double CUTNOISE = 20;//adc
const double CUTCHI2MAX = 10;//track quality
const double CUTCHI2MIN = 0.001;//track quality

//Detector numbers and projection
const int NRPHI = 2;//reads Y coordinate
const int NZ = 1;//reads X coordinate

//globals
TString DETNAME = "";//Rphi or Z

//old rd3 stuff
const int DETEXTRA = 1;// 0 tpot, 1 rd3_1,2 rd3_2, 3 rd3_3
const bool detinfit[] = {false,false,true,true};// tpot, rd3_1, rd3_2, rd3_3
TTree * finaltree;
//Analysis constant
//to do scan todo


double ZDET[] = {748,//TPOT Z
		 755.5,  //TPOT RPHI 
		 285,//RD3 mid
		 235.5};//RD3 low

const double XDET[] = {282.0,//TPOT Z
		       158.8,//TPOT RPHI
		       79.3,//+0.140,//RD3 mid
		       78.4};//RD3 low

double THETA[] = {0.86, //TPOT Z
		  0.82, //TPOT RPHI
		  0.55,//RD3 mid
		  0.55};//RD3 low

double SCALE[] = {1.,//TPOT Z
		  1.,//TPOT RPHI
		  1.,//RD3 mid
		  1.};//RD3 low

const double NSIGMA = 4.5;
const double PITCHRPHI = 1.;//tpot
const double PITCHRZ = 2.;//tpot
const double pitchrd3 = 0.5;
const int nrd3 = 4;

//Global for ray tree
double rayxup;
double rayyup;
double rayzup;
double rayxdn;
double rayydn;
double rayzdn;
double chi2x;
double chi2y;
vector <double> * fovrayx, * fovrayy, * fovchix, * fovchiy,* fovangx, * fovangy;
vector <double> * fovclustmaxamp;
vector <double> * fovclustamp;
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
double getAXray(double z){
  //solve X=Az+B
  double Z_Up= rayzup;  double Z_Down= rayzdn;
  double X_Up= rayxup;  double X_Down= rayxdn;
  if(Z_Up-Z_Down == 0) {cout << "ERROR: Z_up == Z_Down" << endl; return -1;}  
  double A=(X_Up-X_Down)/(Z_Up-Z_Down);
  return A;
}
double getAYray(double z){
  //solve Y=Az+B
  double Z_Up= rayzup;  double Z_Down= rayzdn;
  double Y_Up= rayyup;  double Y_Down= rayydn;
  if(Z_Up-Z_Down == 0) {cout << "ERROR: Z_up == Z_Down" << endl; return -1;}  
  double A=(Y_Up-Y_Down)/(Z_Up-Z_Down);
  return A;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void OuttreeAnaSingle(TString filename, int maxev = -1)
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
  finaltree = new TTree("TF","Event");
  fovrayx=new vector <double>; finaltree->Branch("RayX", &fovrayx);
  fovrayy=new vector <double>; finaltree->Branch("RayY", &fovrayy);
  fovclustmaxamp=new vector <double>; finaltree->Branch("ClustMaxAmp", &fovclustmaxamp);
  fovclustamp=new vector <double>; finaltree->Branch("ClustAmp", &fovclustamp);
  fovclustmaxstrip=new vector <double>; finaltree->Branch("ClustMaxStrip", &fovclustmaxstrip);
  fovclustx=new vector <double>; finaltree->Branch("ClustX", &fovclustx);
  fovclustn=new vector <int>; finaltree->Branch("ClustN", &fovclustn);
  fovclustsize=new vector <int>; finaltree->Branch("ClustSize", &fovclustsize);
  fovextra=new vector <double>; finaltree->Branch("Extra", &fovextra);
  fovextradet=new vector <double>; finaltree->Branch("ExtraX", &fovextradet);
  fovchix=new vector <double>; finaltree->Branch("Chi2X", &fovchix);
  fovchiy=new vector <double>; finaltree->Branch("Chi2Y", &fovchiy);
  fovangx=new vector <double>; finaltree->Branch("AngX", &fovangx);
  fovangy=new vector <double>; finaltree->Branch("AngY", &fovangy);

  /////////////////////////////////////////////////////////////////////////////////////////////////
  //efficiency part
  int nbineff=300;
  TH2F * hefftrphi = new TH2F("hefftrphi","hefftrphi",nbineff,-30,500,nbineff,-30,500);
  TH2F * hefftrz = new TH2F("hefftrz","hefftrz",nbineff,-30,500,nbineff,-30,500);
  TH2F * heffdetphi = new TH2F("heffdetphi","Efficiency RPhi;x(mm);y(mm)",nbineff,-30,500,nbineff,-30,500);
  TH2F * heffdetz = new TH2F("heffdetz","Efficiency Z;x(mm);y(mm)",nbineff,-30,500,nbineff,-30,500);
  TH2F * hamp = new TH2F("hamp","Mean Max Amplitude per Cluster;x(mm);y(mm)",nbineff,-30,500,nbineff,-30,500);
  TH2F * hcl = new TH2F("hcl","Mean Cluster Size;x(mm);y(mm)",nbineff,-30,500,nbineff,-30,500);
  TH2F * hclall = new TH2F("hclall","Amplitude;x(mm);y(mm)",nbineff,-30,500,nbineff,-30,500);
  bool trokphi = false,trokz = false, trrec=false, trdone=false, trdonerphi=false,trdonez=false; 
  int trokcntz=0;int trokcntboth=0;int trokcntphi=0,trseenboth=0 ; 
  int trseencntrphi=0;int trseencntz=0; bool trfilledphi=false;bool trfilledz=false;
  /////////////////////////////////////////////////////////////////////////////////////////////////
  
  for(int iev=0; iev<nev;iev++){
    if(iev%10000==0) cout << iev << endl; 
    treAna->GetEntry(iev);
    trokz=false;trokphi=false; trrec=false; trdonerphi=false;
    trdonez=false; trdone=false;trfilledz = false;trfilledphi = false;
    fovrayx->clear();
    fovrayy->clear();
    fovchix->clear();
    fovchiy->clear();
    fovangx->clear();
    fovangy->clear();
    fovclustmaxstrip->clear();
    fovclustx->clear();
    fovclustmaxamp->clear();
    fovclustamp->clear();
    fovclustn->clear();
    fovclustsize->clear();
    fovextra->clear();
    fovextradet->clear();

    double xrd3[nrd3], ard3[nrd3]; for(int id=0;id<nrd3;id++) {xrd3[id]=-1; ard3[id]=-1;} 
    if(clustn->size()>500) {cout << iev << " " << clustn->size() <<  endl;
      continue;}
    for(unsigned int icl=0;icl<clustn->size();icl++){//cluster loop
      int detn = clustn->at(icl); 
      if(detn<=0 || detn>4) {cout << "det does not exist : " << detn << endl; continue;}
      
      double zdet = ZDET[detn-1];
      double rayx = getXray(zdet);
      double rayy = getYray(zdet);
      //rotation of tracker
      double theta = THETA[detn-1]*TMath::Pi()/180.;//in radian
      double xray = rayx*cos(theta)-rayy*sin(theta);
      double yray = rayx*sin(theta)+rayy*cos(theta);
      double clustxhere = clustx->at(icl)*SCALE[detn-1];
      yray = yray+XDET[NRPHI-1];
      xray = xray+XDET[NZ-1];
      fovrayx->push_back(xray);
      fovrayy->push_back(yray);
      
      //effcalc
      if(chi2x>CUTCHI2MIN && chi2y>CUTCHI2MIN && chi2x<CUTCHI2MAX && chi2y<CUTCHI2MAX)//ok track
	{
	  double thetaz = THETA[NZ-1]*TMath::Pi()/180.;//in radian
	  double xrayz = rayx*cos(thetaz)-rayy*sin(thetaz);
	  double yrayz = rayx*sin(thetaz)+rayy*cos(thetaz);
	  yrayz = yrayz+XDET[NRPHI-1];
	  xrayz = xrayz+XDET[NZ-1];
	  double thetaphi = THETA[NRPHI-1]*TMath::Pi()/180.;//in radian
	  double xrayphi = rayx*cos(thetaphi)-rayy*sin(thetaphi);
	  double yrayphi = rayx*sin(thetaphi)+rayy*cos(thetaphi);
	  yrayphi = yrayphi+XDET[NRPHI-1];
	  xrayphi = xrayphi+XDET[NZ-1];
	  
	  if(!trrec)//might create a pb where the wrong theta is used 
	    {
	      trrec=true;
	      //hefftrphi->Fill(xrayphi, yrayphi);
	      //hefftrz->Fill(xrayz, yrayz);
	      if(xrayz<xmax && xrayz>xmin && yrayz<ymax && yrayz>ymin) //fiducial cut
		{ trokcntz++;trokz=true;}
	      if(xrayphi<xmax && xrayphi>xmin && yrayphi<ymax && yrayphi>ymin) //fiducial cut
		{ trokcntphi++;trokphi=true;}
	      if(xrayz<xmax && xrayz>xmin && yrayz<ymax && yrayz>ymin &&
		 xrayphi<xmax && xrayphi>xmin && yrayphi<ymax && yrayphi>ymin) //fiducial cut both
		{ trokcntboth++;}
	    }
	  
	  if(!trdonerphi  && detn==NRPHI && clustmaxamp->at(icl)>CUTNOISE && fabs(yrayphi-clustxhere)<ROADEFF)
	    {
	      if(!trfilledphi) { /*heffdetphi->Fill(xrayphi, yrayphi);*/ trfilledphi=true;}
	      if(trokphi){trseencntrphi++; trdonerphi=true;} 
	    }
	  
	  if(!trdonez && detn==NZ && clustmaxamp->at(icl)>CUTNOISE && fabs(xrayz-clustxhere)<ROADEFF)//detok
	    { if(!trfilledz) {/*heffdetz->Fill(xray, yray);*/ trfilledz=true;}
	      if(trokz){trseencntz++; trdonez=true;} }
	  
	  if(trdonerphi && trdonez && !trdone) {trseenboth++; trdone=true;}
	  if(icl==clustn->size()-1)//last cluster
	    {
	      heffdetphi->Fill(xrayphi, yrayphi,trfilledphi);
	      heffdetz->Fill(xrayz, yrayz,trfilledz);
	      if(trrec) hefftrphi->Fill(xrayphi, yrayphi,trrec);else hefftrphi->Fill(xrayphi, yrayphi,0.);
	      hefftrz->Fill(xrayz, yrayz,trrec);
	    } 
	}//chi2 trck cut

      if(clustmaxamp->at(icl)>CUTNOISE)
	{
	  hclall->Fill(xray, yray);
	  hamp->Fill(xray, yray, clustmaxamp->at(icl));
	  hcl->Fill(xray, yray, clustsize->at(icl));
	}
      //
      fovchix->push_back(chi2x);
      fovchiy->push_back(chi2y);
      fovangx->push_back(getAXray(zdet));
      fovangy->push_back(getAYray(zdet));

      //cout << iev<< "] " << detn << "\t" << zdet << "\t" << xray << endl;
      fovclustx->push_back(clustxhere);
      fovclustmaxamp->push_back(clustmaxamp->at(icl));
      fovclustamp->push_back(clustamp->at(icl));
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
  //finaltree->Write();

  if(PLOTEFF)
    {
      //gStyle->SetOptFit(000);
      //gStyle->SetOptStat(000);
      //gROOT->LoadMacro("SetStyle.C");
      SetStyle();
      TCanvas * ceff = new TCanvas("ceff","ceff",1400,600);
      ceff->Divide(2,1);
      ceff->cd(1);
      heffdetz->SetStats(0);
      heffdetz->Divide(hefftrz);
      heffdetz->Draw("colz");
      //heffdet->GetZaxis()->SetRangeUser(0,1.);
      TLine * l1 = new TLine(xmin,ymin,xmin,ymax); l1->SetLineStyle(7); l1->SetLineWidth(2); l1->Draw();
      TLine * l2 = new TLine(xmin,ymax,xmax,ymax); l2->SetLineStyle(7); l2->SetLineWidth(2);l2->Draw();
      TLine * l3 = new TLine(xmax,ymax,xmax,ymin); l3->SetLineStyle(7); l3->SetLineWidth(2);l3->Draw();
      TLine * l4 = new TLine(xmax,ymin,xmin,ymin); l4->SetLineStyle(7); l4->SetLineWidth(2);l4->Draw();
      ceff->cd(2);
      heffdetphi->SetStats(0);
      heffdetphi->Divide(hefftrphi);
      heffdetphi->Draw("colz");
      
      l1->SetLineWidth(2); l1->Draw();
      l2->SetLineWidth(2);l2->Draw();
      l3->SetLineWidth(2);l3->Draw();
      l4->SetLineWidth(2);l4->Draw();

      cout << "RPHI : " << trseencntrphi << " / " << trokcntphi << " = " << 1.0*trseencntrphi/trokcntphi << endl;
      cout << "RPHI wZ :" << trseenboth << " / " << trseencntz << " = " << 1.0*trseenboth/trseencntz << endl;

      cout << "Z :" << trseencntz << " / " << trokcntz << " = " << 1.0*trseencntz/trokcntz << endl;
      cout << "Z wRPHI :" << trseenboth << " / " << trseencntrphi << " = " << 1.0*trseenboth/trseencntrphi << endl;

      cout << "MOD :" << trseenboth << " / " << trokcntboth << " = " << 1.0*trseenboth/trokcntboth << endl;
      //      cout << trseencnt << " / " << trokcnt << " = " << 1.0*trseencnt/trokcnt << endl;

      TCanvas * camp = new TCanvas("camp","camp",1000,600);
      hamp->Divide(hclall);//todo list of impart plots
      hamp->Draw("colz");
      hamp->GetZaxis()->SetRangeUser(0,500);

      //TCanvas * ccl2d = new TCanvas("ccl2d","ccl2d",1000,600);
      //hcl->Divide(hclall);//todo list of impart plots
      //hcl->Draw("colz");

      //residuals :
      // gStyle->SetOptFit(111);
      // gStyle->SetOptStat(1111);
      TCanvas * cres = new TCanvas("cres","cres",1400,400);
      cres->Divide(4,1);
      cres->cd(1);
      TString cutstring =  TString::Format("ClustN==%d && ClustMaxAmp>%f && Chi2Y<%f && Chi2X<%f && Chi2X>%f && Chi2Y>%f" 
					   ,NRPHI,       CUTNOISE,   CUTCHI2MAX,CUTCHI2MAX,CUTCHI2MIN,CUTCHI2MIN );
      cout << "CUTS : " << cutstring << endl;
      TH1F * hres1 = new TH1F("hres1", "Residuals Phi;X_{RAY}-X_{DET}(mm);#;", 200,-10,10);
      finaltree->Draw("RayY-ClustX>>hres1",cutstring,"");      
      hres1->Fit("gaus","","",-2,2);
      hres1->SetStats(0);
      hres1->Draw();
      TF1 *st = (TF1*)hres1->FindObject("gaus");
      st->GetParameter(2);
      TString sigma = TString::Format("#sigma_{X} = %.2fmm",st->GetParameter(2));
      double xtp = (hres1->GetXaxis()->GetXmax()-hres1->GetXaxis()->GetXmin())*0.6+hres1->GetXaxis()->GetXmin();
      double ytp = hres1->GetMaximum()*0.8;
      TLatex * tp = new TLatex(xtp,ytp,sigma);
      tp->Draw();
      cres->cd(2);
      TH1F * hclsize = new TH1F("hclsize", "Cluster Size Phi;Cl size (strips)""", 20,-0.5,19.5);
      finaltree->Draw("ClustSize>>hclsize",cutstring,"");//,TString::Format("ClustN==%d && ClustMaxAmp>200 ",DETEFF),"");
      hclsize->SetStats(0);
      cres->GetPad(2)->SetLogy();
      cres->cd(3);
      TH1F * hmaxamp = new TH1F("hmaxamp", "Max Hit Amplitude Phi;Max. Amp (adc)", 300,0,5000);
      finaltree->Draw("ClustMaxAmp>>hmaxamp",cutstring,"");//,TString::Format("ClustN==%d && ClustMaxAmp>200",DETEFF),"");
      hmaxamp->SetStats(0);
      hmaxamp->GetXaxis()->SetNdivisions(508,true);
      cres->GetPad(3)->SetLogy();
      cres->cd(4);
      TH1F * hclamp = new TH1F("hclamp", "Cluster Amplitude Phi;Amp (adc);#", 300,0,20000);
      finaltree->Draw("ClustAmp>>hclamp",cutstring,"");//,TString::Format("ClustN==%d && ClustMaxAmp>200",DETEFF),"");
      hclamp->SetStats(0);
      hclamp->GetXaxis()->SetNdivisions(508,true);
      cres->GetPad(4)->SetLogy();

      TCanvas * cresz = new TCanvas("cresz","cresz",1400,400);
      cresz->Divide(4,1);
      cresz->cd(1);
      TString cutstringz =  TString::Format("ClustN==%d && ClustMaxAmp>%f && Chi2Y<%f && Chi2X<%f && Chi2X>%f && Chi2Y>%f" 
					   ,NZ,       CUTNOISE,   CUTCHI2MAX,CUTCHI2MAX,CUTCHI2MIN,CUTCHI2MIN );
      cout << "CUTS Z: " << cutstringz << endl;
      TH1F * hres1z = new TH1F("hres1z", "Residuals Z;X_{RAY}-X_{DET}(mm);#;", 200,-10,10);
      finaltree->Draw("RayX-ClustX>>hres1z",cutstringz,"");
      hres1z->Fit("gaus","","",-2,2);
      TF1 *stz = (TF1*)hres1z->FindObject("gaus");
      stz->GetParameter(2);
      TString sigmaz = TString::Format("#sigma_{X} = %.2fmm",stz->GetParameter(2));
      double xtpz = (hres1z->GetXaxis()->GetXmax()-hres1z->GetXaxis()->GetXmin())*0.6+hres1z->GetXaxis()->GetXmin();
      double ytpz = hres1z->GetMaximum()*0.8;
      TLatex * tpz = new TLatex(xtpz,ytpz,sigmaz);
      tpz->Draw();



      cresz->cd(2);
      TH1F * hclsizez = new TH1F("hclsizez", "Cluster Size Z;Cl size (strips)""", 20,-0.5,19.5);
      finaltree->Draw("ClustSize>>hclsizez",cutstringz,"");//,TString::Format("ClustN==%d && ClustMaxAmp>200 ",DETEFF),"");
      hclsizez->SetStats(0);
      cresz->GetPad(2)->SetLogy();
      cresz->cd(3);
      TH1F * hmaxampz = new TH1F("hmaxampz", "Max Hit Amplitude Z;Max. Amp (adc)", 300,0,5000);
      finaltree->Draw("ClustMaxAmp>>hmaxampz",cutstringz,"");//,TString::Format("ClustN==%d && ClustMaxAmp>200",DETEFF),"");
      cresz->GetPad(3)->SetLogy();
      hmaxampz->SetStats(0);
      hmaxampz->GetXaxis()->SetNdivisions(508,true);
      cresz->cd(4);
      TH1F * hclampz = new TH1F("hclampz", "Cluster Amplitude Z;Amp (adc);#", 300,0,20000);
      finaltree->Draw("ClustAmp>>hclampz",cutstringz,"");//,TString::Format("ClustN==%d && ClustMaxAmp>200",DETEFF),"");
      hclampz->SetStats(0);
      hclampz->GetXaxis()->SetNdivisions(508,true);
      cresz->GetPad(4)->SetLogy();

      if(SAVEEFF)
	{
	  TString pictfolder = "pict";
	  
	  ceff->SaveAs(TString::Format("%s/eff2d_%s.pdf", pictfolder.Data(), DETNAME.Data()));
	  cres->SaveAs(TString::Format("%s/plot1d_Phi_%s.pdf", pictfolder.Data(), DETNAME.Data()));
	  cresz->SaveAs(TString::Format("%s/plot1d_Z_%s.pdf", pictfolder.Data(), DETNAME.Data()));
	  camp->SaveAs(TString::Format("%s/amp2d_%s.pdf", pictfolder.Data(), DETNAME.Data()));
	  ofstream myfile;
	  myfile.open (TString::Format("%s/eff.txt", pictfolder.Data()));

	  myfile << "RPHI : " << trseencntrphi << " / " << trokcntphi << " = " << 1.0*trseencntrphi/trokcntphi << endl;
	  myfile << "RPHI wZ :" << trseenboth << " / " << trseencntz << " = " << 1.0*trseenboth/trseencntz << endl;
	  
	  myfile << "Z :" << trseencntz << " / " << trokcntz << " = " << 1.0*trseencntz/trokcntz << endl;
	  myfile << "Z wRPHI :" << trseenboth << " / " << trseencntrphi << " = " << 1.0*trseenboth/trseencntrphi << endl;
	  
	  myfile << "MOD :" << trseenboth << " / " << trokcntboth << " = " << 1.0*trseenboth/trokcntboth << endl;

	  myfile.close();
	}
      
    }

  return;
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
  //return;

  for(int idet=1;idet<=1/*nrd3*/;idet++)
    {
      TCanvas * cali = new TCanvas(TString::Format("cali%d",idet),TString::Format("cali%d",idet),
				   1200,600);
      cali->Divide(2,1);
      cali->cd(1);
      TH2F * hali1 = new TH2F(TString::Format("hali%d",idet),TString::Format("hali%d",idet),
			      200,0,500,200,-30,30);
      finaltree->Draw(TString::Format("RayX-ClustX:RayX>>hali%d",idet),TString::Format("ClustN==%d && ClustMaxAmp>100",idet),"colz");
      hali1->Draw("colz");
      cali->cd(2);
      TH2F * hali2 = new TH2F(TString::Format("hali2%d",idet),TString::Format("hali2%d",idet),
			      200,-250,250,200,-30,30);
      finaltree->Draw(TString::Format("RayX-ClustX:RayY>>hali2%d",idet),TString::Format("ClustN==%d && ClustMaxAmp>100",idet),"colz");
      
    }
  

}

void OuttreeAna(TString filename, int layer = 1, int maxev = -1, bool doscan=false)
{  
  int detscan =0;
  if(layer==1) //Z
    {
      detscan = NZ-1;//reference for centering the scan
      DETNAME = "Z";
      DETEFF = NZ; 
    }
  else {
    detscan = NRPHI-1;//reference for centering the scan
    DETNAME = "Rphi";
    DETEFF = NRPHI; 
  }
  if(!doscan)
    {
      OuttreeAnaSingle(filename,maxev);
      return;
    }
  //SCAN PARAM //todo
  double step = 0.5;
  double nstep = 10;
  double from = ZDET[detscan]-step*nstep/2; 
  double to = ZDET[detscan]+step*(nstep+1)/2;
  //double from = THETA[detscan]-step*nstep/2; 
  //double to = THETA[detscan]+step*(nstep+1)/2;
  //double from = SCALE[detscan]-step*nstep/2; 
  //double to = SCALE[detscan]+step*(nstep+1)/2;

  TProfile * scanhistophi = new TProfile("scanhistophi","scanhistophi",(to-from)/step+1,from-step,to+step);
  TProfile * scanhistoz = new TProfile("scanhistoz","scanhistoz",(to-from)/step+1,from-step,to+step);
  TCanvas * c1 = new TCanvas("c1","c1");
  for(double param = from; param<=to; param+=step)
    {
      ZDET[NZ-1]=param; ZDET[NRPHI-1]=param;
      //THETA[NZ-1]=param; THETA[NRPHI-1]=param;
      //SCALE[NZ-1]=param; SCALE[NRPHI-1]=param;
      cout << "PARAM[" << detscan+1 << "] = " << param <<endl; 
      OuttreeAnaSingle(filename,maxev);
      //z
      TH1F * balz = new TH1F("balz","balz",300,-50,50);//to be tuned
      balz->Draw();
      finaltree->Draw("RayX-ClustX>>balz",TString::Format("ClustN==%d && ClustMaxAmp>100 && Chi2X<10 && Chi2X>0.01",NZ).Data(),"");
      balz->Fit("gaus","","");//to be tuned
      balz->Draw();
      TF1 * fitz = balz->GetFunction("gaus");
      scanhistoz->Fill(param,fitz->GetParameter(2));
      delete fitz;
      delete balz;
      //phi
      TH1F * balphi = new TH1F("balphi","balphi",300,-30,30);//to be tuned
      balphi->Draw();
      finaltree->Draw("RayY-ClustX>>balphi",TString::Format("ClustN==%d && ClustMaxAmp>200 && Chi2Y>0.01 && Chi2Y<10",NRPHI).Data(),"");
      //else finaltree->Draw("RayX-ClustX>>bal",TString::Format("ClustN==%d && ClustMaxAmp>100 && Chi2X<0.1 && ClustX<230",detscan+1).Data(),"");
      balphi->Fit("gaus","","");//,-10,10);//to be tuned
      balphi->Draw();
      TF1 * fitphi = balphi->GetFunction("gaus");
      scanhistophi->Fill(param,fitphi->GetParameter(2));
      delete fitphi;
      delete balphi;

      finaltree->Clear();
    }
  TCanvas * cphi = new TCanvas("cphi","cphi");
  scanhistophi->Draw();
  TCanvas * cz = new TCanvas("cz","cz");
  scanhistoz->Draw();
  return;
}


const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7}; // for syst bands
const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};

void SetStyle(Bool_t graypalette) {
  cout << "Setting style!" << endl;
  
  gStyle->Reset("Plain");
  //gStyle->SetOptTitle(0);
  //gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kBlue+2);
  //gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kRed);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y"); 

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);


}


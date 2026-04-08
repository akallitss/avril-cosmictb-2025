#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLine.h"
#include "TProfile.h"
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
  int nbineff=400;
  TH2F * hefftr = new TH2F("hefftr","hefftr",nbineff,-30,500,nbineff,-30,500);
  TH2F * heffdet = new TH2F("heffdet","Efficiency;x(mm);y(mm)",nbineff,-30,500,nbineff,-30,500);
  TH2F * hamp = new TH2F("hamp","Mean Max Amplitude per Cluster;x(mm);y(mm)",nbineff,-30,500,nbineff,-30,500);
  TH2F * hcl = new TH2F("hcl","Mean Cluster Size;x(mm);y(mm)",nbineff,-30,500,nbineff,-30,500);
  TH2F * hclall = new TH2F("hclall","Amplitude;x(mm);y(mm)",nbineff,-30,500,nbineff,-30,500);
  bool trok = false, trrec=false, trdone=false, trdonerphi=false,trdonez=false; int trokcnt=0,trseenboth=0 ; 
  int trseencntrphi=0;int trseencntz=0; bool trfilled=false;
  /////////////////////////////////////////////////////////////////////////////////////////////////
  
  for(int iev=0; iev<nev;iev++){
    if(iev%10000==0) cout << iev << endl; 
    treAna->GetEntry(iev);
    trok=false; trrec=false; trdonerphi=false;trdonez=false; trdone=false;trfilled = false;
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
      //
      //      if(!trdone)
      //	{
      if(chi2x>CUTCHI2MIN && chi2y>CUTCHI2MIN && chi2x<CUTCHI2MAX && chi2y<CUTCHI2MAX)//ok track
	{
	  if(!trrec)//might create a pb where the wrong theta is used 
	    {
	      trrec=true;
	      double thetaeff = THETA[DETEFF-1]*TMath::Pi()/180.;//in radian
	      double xrayeff = rayx*cos(thetaeff)-rayy*sin(thetaeff);
	      double yrayeff = rayx*sin(thetaeff)+rayy*cos(thetaeff);
	      //double clustxhere = clustx->at(icl)*SCALE[detn-1];
	      yrayeff = yrayeff+XDET[NRPHI-1];
	      xrayeff = xrayeff+XDET[NZ-1];
	      hefftr->Fill(xrayeff, yrayeff);
	      if(xrayeff<xmax && xrayeff>xmin && yrayeff<ymax && yrayeff>ymin) //fiducial cut
		{ trokcnt++;trok=true;}
	    }
	}
	  
      if(trrec==true)
	{
	  if(!trdonerphi  && detn==NRPHI && clustmaxamp->at(icl)>CUTNOISE && fabs(yray-clustxhere)<ROADEFF)
	    {
	      if(!trfilled && DETEFF==NRPHI) {heffdet->Fill(xray, yray); trfilled=true;}
	      if(trok){trseencntrphi++; trdonerphi=true;} }
	  
	  if(!trdonez && detn==NZ && clustmaxamp->at(icl)>CUTNOISE && fabs(xray-clustxhere)<ROADEFF)//detok
	    { if(!trfilled && DETEFF==NZ) {heffdet->Fill(xray, yray); trfilled=true;}
	      if(trok){trseencntz++; trdonez=true;} }
	  
	  if(trdonerphi && trdonez && !trdone) {trseenboth++; trdone=true;}
	}
      //}//trdone
      
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
      gStyle->SetOptFit(000);
      gStyle->SetOptStat(000);
      TCanvas * ceff = new TCanvas("ceff","ceff",1000,600);
      heffdet->Divide(hefftr);
      heffdet->Draw("colz");
      //heffdet->GetZaxis()->SetRangeUser(0,1.);
      TLine * l1 = new TLine(xmin,ymin,xmin,ymax); l1->SetLineStyle(7); l1->SetLineWidth(2); l1->Draw();
      TLine * l2 = new TLine(xmin,ymax,xmax,ymax); l2->SetLineStyle(7); l2->SetLineWidth(2);l2->Draw();
      TLine * l3 = new TLine(xmax,ymax,xmax,ymin); l3->SetLineStyle(7); l3->SetLineWidth(2);l3->Draw();
      TLine * l4 = new TLine(xmax,ymin,xmin,ymin); l4->SetLineStyle(7); l4->SetLineWidth(2);l4->Draw();
      cout << "RPHI : " << trseencntrphi << " / " << trokcnt << " = " << 1.0*trseencntrphi/trokcnt << endl;
      cout << "RPHI wZ :" << trseenboth << " / " << trseencntz << " = " << 1.0*trseenboth/trseencntz << endl;

      cout << "Z :" << trseencntz << " / " << trokcnt << " = " << 1.0*trseencntz/trokcnt << endl;
      cout << "Z wRPHI :" << trseenboth << " / " << trseencntrphi << " = " << 1.0*trseenboth/trseencntrphi << endl;

      //      cout << trseencnt << " / " << trokcnt << " = " << 1.0*trseencnt/trokcnt << endl;

      TCanvas * camp = new TCanvas("camp","camp",1000,600);
      hamp->Divide(hclall);//todo list of impart plots
      hamp->Draw("colz");
      hamp->GetZaxis()->SetRangeUser(0,500);

      //TCanvas * ccl2d = new TCanvas("ccl2d","ccl2d",1000,600);
      //hcl->Divide(hclall);//todo list of impart plots
      //hcl->Draw("colz");

      //residuals :
      gStyle->SetOptFit(111);
      gStyle->SetOptStat(1111);
      TCanvas * cres = new TCanvas("cres","cres",1400,400);
      cres->Divide(4,1);
      cres->cd(1);
      TString cutstring =  TString::Format("ClustN==%d && ClustMaxAmp>%f && Chi2Y<%f && Chi2X<%f && Chi2X>%f && Chi2Y>%f" 
					   ,DETEFF,       CUTNOISE,   CUTCHI2MAX,CUTCHI2MAX,CUTCHI2MIN,CUTCHI2MIN );
      cout << "CUTS : " << cutstring << endl;
      TH1F * hres1 = new TH1F("hres1", "Residuals;X_{RAY}-X_{DET}(mm);#;", 200,-10,10);
      if(DETEFF==NRPHI) finaltree->Draw("RayY-ClustX>>hres1",cutstring,"");      
      else  if(DETEFF==NZ)finaltree->Draw("RayX-ClustX>>hres1",cutstring,"");
      hres1->Fit("gaus","","",-2,2);
      cres->cd(2);
      TH1F * hclsize = new TH1F("hclsize", "Cluster Size;Cl size (strips)""", 20,-0.5,19.5);
      finaltree->Draw("ClustSize>>hclsize",cutstring,"");//,TString::Format("ClustN==%d && ClustMaxAmp>200 ",DETEFF),"");
      cres->GetPad(2)->SetLogy();
      cres->cd(3);
      TH1F * hmaxamp = new TH1F("hmaxamp", "Max Hit Amplitude;Max. Amp (adc)", 300,0,5000);
      finaltree->Draw("ClustMaxAmp>>hmaxamp",cutstring,"");//,TString::Format("ClustN==%d && ClustMaxAmp>200",DETEFF),"");
      cres->GetPad(3)->SetLogy();
      cres->cd(4);
      TH1F * hclamp = new TH1F("hclamp", "Cluster Amplitude;Amp (adc);#", 300,0,20000);
      finaltree->Draw("ClustAmp>>hclamp",cutstring,"");//,TString::Format("ClustN==%d && ClustMaxAmp>200",DETEFF),"");
      cres->GetPad(4)->SetLogy();

      if(SAVEEFF)
	{
	  TString pictfolder = "pict";
	  
	  ceff->SaveAs(TString::Format("%s/eff2d_%s.pdf", pictfolder.Data(), DETNAME.Data()));
	  cres->SaveAs(TString::Format("%s/plot1d_%s.pdf", pictfolder.Data(), DETNAME.Data()));
	  camp->SaveAs(TString::Format("%s/amp2d_%s.pdf", pictfolder.Data(), DETNAME.Data()));
	  ofstream myfile;
	  myfile.open (TString::Format("%s/eff_%s.txt", pictfolder.Data(), DETNAME.Data()));
	  if(DETEFF==NZ)
	    myfile << "EFF = " << trseenboth << " / " << trseencntrphi << " = " << 1.0*trseenboth/trseencntrphi << endl;
	  else if(DETEFF==NRPHI)
	    myfile << "EFF = " << trseenboth << " / " << trseencntz << " = " << 1.0*trseenboth/trseencntz << endl;
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

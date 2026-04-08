#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
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

//Analysis constant
const double zdet = 666+25;
const double nsigma = 4.5;
const double pitch = 2.;
const double THETA = 0.9;//degrees

//for efficency
const double xdet = 2.78681e+02;
const double xres = 4.5;


//Global for ray tree
TTree * treRay; TFile * fray;
int evn;
double evttime;
int rayN;
double Z_Up;
double Z_Down;
vector <double>* X_Up;
vector <double>* X_Down;
vector <double>* Y_Up;
vector <double>* Y_Down;
vector <double>* Chi2X;
vector <double>* Chi2Y; 

//Globals for data tree
TTree * treData; TFile * fdata;
int fevnData;
int fNdet, fNstrip, fNsample;
float *data;
//float data[16][64][32];

//RMS:
float *RMS;

//Hitos
TH2F *hRayXY, *hRayXYseen;
TH2F *hRayXDetX, *hRayYDetX;
TH1F *hDetX;

//Out Tree
TFile * foutfile;
TTree * fouttree;
vector <double> * fovdetx;
vector <int> * fovstrip;
vector <int> * fovdet;
vector <int> * fovsample;
vector <double> * fovdetamp;
//vector <double> * fovrayx;
//vector <double> * fovrayy;
double fovrayx;
double fovrayy;
double fochi2x;
double fochi2y;

int foDataEvId, foRayEvId;
int foIdData;


void ReadRMS(string fileRMS, bool saveintree=true){
  // ifstream in;
  // in.open("Pedestal.dat");
  // int nlines=0;
  // int det, strip;
  // float pedestal_value;
  // float Pedestal[Ndet][Nstrip];
  // while (1) { // read the text file
  //   in >> det >> strip >> pedestal_value;
  //   Pedestal[det][strip] = pedestal_value;
  //   if (!in.good()) break;
  //   nlines++;
  // }
  
    
  double fovrms[fNdet][fNstrip];
  char strtmp[50]; sprintf(strtmp,"RMS[%d][%d]/F",fNdet,fNstrip);
  TBranch *newBranch ;  
  if(saveintree) newBranch = fouttree->Branch("RMS", fovrms, strtmp);

  ifstream inrms;
  int nlines=0; int det, strip;
  float rms_value;
  RMS = new float[fNdet*fNstrip]();
  inrms.open(fileRMS.c_str());
  while (1) { // read the text file
    inrms >> det >> strip >> rms_value;
    RMS[fNstrip*det+strip] = rms_value;
    if(saveintree) fovrms[det][strip]=rms_value;
    //cout << det << " " <<  strip << " " << RMS[fNstrip*det+strip]  << endl;    
    if (!inrms.good()) break;
    nlines++;
  }
  if(saveintree)    {newBranch->Fill();}
}


void InitHisto(){
  double xmin = -500, xmax = 500; int nbinx =200;
  double xdmin = -50, xdmax = 600; int nbindx =200;
  double ymin = -500, ymax = 500; int nbiny =200;
  hRayXY = new TH2F("hRayXY","hRayXY;x[mm];y[mm]",nbinx,xmin,xmax,nbiny,ymin,ymax);
  hRayXYseen = new TH2F("hRayXYseen","hRayXYseen;x[mm];y[mm]",nbinx,xmin,xmax,nbiny,ymin,ymax);
  hRayXDetX= new TH2F("hRayXDetX","hRayXDetX;xray[mm];xdet[mm]",nbinx,xmin,xmax,nbindx,xdmin,xdmax);
  hRayYDetX= new TH2F("hRayYDetX","hRayYDetX;yray[mm];xdet[mm]",nbiny,ymin,ymax,nbindx,xdmin,xdmax);
  hDetX=new TH1F("hDetX","hDetX;xdet[mm]",nbindx,xdmin,xdmax);

  //ttree
  TFile * foutfile = new TFile("outtree.root","recreate");
  fouttree = new TTree("Tout","Event");
  fouttree->Branch("EvIdData", &foDataEvId); // event number
  fouttree->Branch("EvIdRay", &foRayEvId); // event number
  fovdetx=new vector <double>;
  fouttree->Branch("DetX", &fovdetx);
  fovstrip=new vector <int>; fouttree->Branch("StripNb", &fovstrip);
  fovdet=new vector <int>; fouttree->Branch("Det", &fovdet);
  fovsample=new vector <int>; fouttree->Branch("Sample", &fovsample);
  fovdetamp=new vector <double>;
  fouttree->Branch("DetAmp", &fovdetamp);
  //fovrayx=new vector <double>;
  fovrayx=0; fouttree->Branch("RayX", &fovrayx);
  //fovrayy=new vector <double>;
  fovrayy=0; fouttree->Branch("RayY", &fovrayy);
  fochi2x=0; fouttree->Branch("Chi2X", &fochi2x);
  fochi2y=0; fouttree->Branch("Chi2Y", &fochi2y);
}

void Draw(){//to display stuf at the end
  hRayXY->Draw("colz");
}

int InitRayTreeRead(TString filename){
  cout << "Read Ray tree : " << filename << endl; 
  fray = new TFile(filename.Data());
  treRay =  (TTree*)fray->Get("T");  
  cout << "Tree = " << treRay->GetName() << "\t"<< treRay->GetEntries() << " events" << endl;
  if(!treRay->GetEntries()) {cout << "FATAL : empty tree " << endl;return -1;}
  evn = 0; treRay->SetBranchAddress("evn",&evn);
  evttime=0.; treRay->SetBranchAddress("evttime",&evttime);
  rayN=0; treRay->SetBranchAddress("rayN",&rayN);
  Z_Up=0; treRay->SetBranchAddress("Z_Up",&Z_Up);
  Z_Down=0; treRay->SetBranchAddress("Z_Down",&Z_Down);
  X_Up=0; treRay->SetBranchAddress("X_Up",&X_Up);
  X_Down=0; treRay->SetBranchAddress("X_Down",&X_Down);
  Y_Up=0; treRay->SetBranchAddress("Y_Up",&Y_Up);
  Y_Down=0; treRay->SetBranchAddress("Y_Down",&Y_Down);
  Chi2X=0; treRay->SetBranchAddress("Chi2X",&Chi2X);
  Chi2Y=0; treRay->SetBranchAddress("Chi2Y",&Chi2Y);

  int iiinit = 5;//display
  for(int iev=0;iev<iiinit;iev++)
    {
      treRay->GetEntry(iev);
      cout << "EVENT "<<iev<<" (" << evn << "): t= " << evttime 
	   << "\tNray=" << rayN << "("<< Z_Up<<","<< Z_Down<< ") => "
	   << X_Up->size() << " "<< X_Down->size() << " " 
	   << Y_Up->size() << " "<< X_Down->size() << endl;
    }
  //f->Close();
  return 0;
}

int InitDataTree(TString filename){
  cout << "Read Data tree : " << filename << endl; 
  fdata = new TFile(filename.Data());
  treData = 0; treData =  (TTree*)fdata->Get("T"); 
  if(!treData) {cout << "empty data file "<< filename<< endl;return -1;} 
  cout << "Data Tree = " << treData->GetName() << "\t"<< treData->GetEntries() << " events" << endl;
  if(!treData->GetEntries()) {cout << "FATAL : empty data tree " << endl;return -1;}
  fevnData = 0; treData->SetBranchAddress("IDEvent",&fevnData);
  const char * branchname =  treData->GetBranch("StripAmpl_corrped")->GetTitle();  
  sscanf(branchname,"StripAmpl_corrped[%d][%d][%d]/F",&fNdet,&fNstrip,&fNsample);
  cout << "Data have " << fNdet << " DREAM, " << fNstrip << " strips, " 
       << fNsample << " samples." << endl;
  //data = new float[16][64][32]; 
  //float ** aaa[16][64] = new float aaa[16][64][32];
  //float aaa[fNdet][fNstrip][fNsample];
  data = new float[fNdet*fNstrip*fNsample]();
 
  treData->SetBranchAddress("StripAmpl_corrped",data); 
  
  //int iiinit = 5;//display
  //for(int iev=0;iev<iiinit;iev++)
  // {
  int iev=1;
  treData->GetEntry(iev);//[3][2][1]
  cout << "DATA EVENT "<<iev<<" (" << fevnData << "): data[0][0][0]= " << data[1*fNsample+2*(fNstrip+3*fNdet)] <<endl;
  for(int i=0;i<10;i++)
    cout << data[i+fNsample*(22+4*fNstrip)] << " ";
  cout << endl;
      // }
  //f->Close();
  return 0;
}

double getXray(double z, int nray){
  //solve X=Az+B
  if(Z_Up-Z_Down == 0) {cout << "ERROR: Z_up == Z_Down" << endl; return -1;}  
  double A=(X_Up->at(nray)-X_Down->at(nray))/(Z_Up-Z_Down);
  double B=X_Up->at(nray)-A*Z_Up;
  double x=A*z+B;//tadaaaa
  return x;
}

double getXdet(int idet,int istr){
  return (idet*fNstrip+istr)*pitch;
}
double getYray(double z, int nray){
  //solve Y=Az+B
  if(Z_Up-Z_Down == 0) {cout << "ERROR: Z_up == Z_Down" << endl; return -1;}  
  double A=(Y_Up->at(nray)-Y_Down->at(nray))/(Z_Up-Z_Down);
  double B=Y_Up->at(nray)-A*Z_Up;
  double y=A*z+B;//tadaaaa
  return y;
}

//-------------------------------------------------
//------------------ACTUAL ANALYSIS PART-----------
void DoAnalysis(int maxevn){
  int ntr=0, ndettr=0; bool seen = false;
  int ncycle=0, ievdata=0; int dist2look = 10;
  bool lostsync = false;
  int nerrcount = 0;
  int nev = treRay->GetEntries();
  int nevd = treData->GetEntries();
  nev = (maxevn<nev && maxevn>0) ? maxevn:nev;
  nev = (nevd<nev) ? nevd:nev;
  int nevdatamax=nevd;
  cout << "building index table for tree synchronisation... "  << endl;
  int turn = 0; int turniev = 0;
  if(0)
  for(int iev=0; iev<nev;iev++){
    treRay->GetEntry(iev);//nev
    treData->GetEntry(iev);
    if(turniev>=4095) {turn++; turniev=0;} else turniev++;
    int realevndata = fevnData + turn*4096 ;//
    cout << "[" << iev << "] = " << evn << " == " <<  realevndata << "\t" <<  fevnData << endl ;
  }
  turn = 0;turniev = 0; int skipped = 0;
  cout << "Analysis loop : " << nev << endl;
  for(int iev=0; iev<nev;iev++){
    treRay->GetEntry(iev);
    treData->GetEntry(iev);
    if(iev%1000 == 0) cout << "." << flush;
    if(turniev>=4095) {turn++; turniev=0;} else turniev++;
    int realevndata = fevnData + turn*4096 ;
    if(evn != realevndata) {skipped++; continue;}
      
    for(int iray=0;iray<rayN;iray++)
      {	
	fovdetx->clear(); fovdetamp->clear();
	fovdet->clear(); fovstrip->clear();fovsample->clear();
	foDataEvId=fevnData;
	foRayEvId=evn;

	double rayx = getXray(zdet,iray);
	double rayy = getYray(zdet,iray);
	
	//rotation of tracker
	double theta = THETA*TMath::Pi()/180.;//in radian
	double rayxprime = rayx*cos(theta)-rayy*sin(theta);
	double rayyprime = rayx*sin(theta)+rayy*cos(theta);
	
	rayx=rayxprime; rayy=rayyprime;

	hRayXY->Fill(rayx, rayy);
	// fovrayx->push_back(getXray(zdet,iray));//find a way to get a bijection; only one track??
	// fovrayy->push_back(getYray(zdet,iray));
	fovrayx = rayx;
	fovrayy = rayy;
	fochi2x = Chi2X->at(iray);
	fochi2y = Chi2Y->at(iray);
	if(fovrayy>-250 && fovrayy<250 && fovrayx>-250 && fovrayx<250  && fochi2x<20 && fochi2y<20) // effiency test ok
	  {seen = false; ntr++;}
	else seen =true;
	for(int idet=0;idet<fNdet;idet++)//detloop
	  for(int istr=0;istr<fNstrip;istr++)
	    {
	      double val=0; double maxsample = 0; double time = -1;       
	      for(int is=0;is<fNsample;is++)
		{
		  if(data[is+fNsample*(istr+fNstrip*idet)]>val) {
		    val = data[is+fNsample*(istr+fNstrip*idet)];//compute max
		    maxsample = is;
		     //compute time at maximum
		    if(is>0 && is<fNsample-1)
		      {
			//Fit
			double atemp = 0.5*(data[is+1+fNsample*(istr+fNstrip*idet)]-2.*data[is+fNsample*(istr+fNstrip*idet)]
			 		    + data[is-1+fNsample*(istr+fNstrip*idet)]);
			double btemp = (data[is+fNsample*(istr+fNstrip*idet)]-data[is-1+fNsample*(istr+fNstrip*idet)])-atemp*(2.*is-1.);
			time = -0.5*btemp/atemp;
			//time += fFineTimeStamp[idet/8]*1./6.;//Coorection avec le Fine Time Stamp
		      }
		    else time=-1;
		  }
		}
	      if(val>=nsigma*RMS[istr+fNstrip*idet])//ZS
		{
		  if(val>40 && getXdet(idet,istr)-fovrayy-xdet>-4.*xres && getXdet(idet,istr)-fovrayy-xdet<4.*xres && !seen) 
		    {seen=true; 
		      ndettr++;		      
		      hRayXYseen->Fill(fovrayx,fovrayy);
		    }//efficiency cut
		  else hRayXYseen->Fill(fovrayx,fovrayy, 0.);
		  fovstrip->push_back(istr);
		  fovdet->push_back(idet);
		  fovsample->push_back(maxsample);
		  fovdetx->push_back(getXdet(idet,istr));
		  fovdetamp->push_back(val);
		  hDetX->Fill(getXdet(idet,istr));
		  if(seen)
		    {
		      hRayXDetX->Fill(fovrayx,getXdet(idet,istr));
		      hRayYDetX->Fill(fovrayy,getXdet(idet,istr));
		    }
		}
	    }
	fouttree->Fill();
      }//nray
  }//iev loop
  fouttree->Write();
  cout << endl <<nev << "events, " << skipped << "skipped." << endl;
  cout << " EFF = " << ndettr << "/" << ntr << " = " << 1.*ndettr/ntr << endl;
  TCanvas * ceff = new TCanvas("ceff","ceff");
  ceff->Divide(2,1);
  ceff->cd(1);
  hRayYDetX->Draw("colz");
  ceff->cd(2);
  hRayXYseen->Divide(hRayXY);
  hRayXYseen->Draw("colz");

  

}//DoAnalysis


//---------------------------------------------------
//-----------MAIN------------------------------------
void TBanalysis(int maxevn = -1) {//todo args...
  InitRayTreeRead("run_rays.root");
  InitDataTree("CodeDamien/output.root");
  InitHisto();
  ReadRMS("CodeDamien/RMSPed.dat",false);
  DoAnalysis(maxevn);
  // Draw();
}

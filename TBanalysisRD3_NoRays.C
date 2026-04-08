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

//DEBUG
bool DEBUGCL = false; //printout clustering

//Analysis constant
const double zdet = 666+25;//TPOT
//const double zdet = 210;
const double nsigma = 4.5;
const double pitch = 2.;
const double pitchrd3 = 0.5;
//const double THETA = 0.9;//degrees //TPOT
const double THETA = 0;//degrees

//for efficency
const double xdet = 2.78681e+02;
const double xres = 4.5;//mm


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
vector <int> * fovdetn;
vector <int> * fovsample;
vector <double> * fovdetamp;
vector <double> * fovdettime;
vector <double> * fovclustamp,* fovclustmaxamp, *fovclustmaxstrip,  * fovclusttime,  * fovclustx;
vector <int>  * fovclustsize, *fovclustn;

double fovrayx;
double fovrayy;
double fovrayxup;
double fovrayyup;
double fovrayzup;
double fovrayxdn;
double fovrayydn;
double fovrayzdn;
double fochi2x;
double fochi2y;

int foDataEvId, foRayEvId;
int foIdData;


void ReadRMS(string fileRMS, bool saveintree=true){
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
  
  //detector hits 
  fovdetx=new vector <double>;fouttree->Branch("DetX", &fovdetx);
  fovdetn=new vector <int>;fouttree->Branch("DetN", &fovdetn);
  fovstrip=new vector <int>; fouttree->Branch("StripNb", &fovstrip);
  fovdet=new vector <int>; fouttree->Branch("Det", &fovdet);
  fovsample=new vector <int>; fouttree->Branch("Sample", &fovsample);
  fovdetamp=new vector <double>;fouttree->Branch("DetAmp", &fovdetamp);
  fovdettime=new vector <double>;fouttree->Branch("DetTime", &fovdettime);
  //detector clusters
  fovclustamp=new vector <double>;fouttree->Branch("ClustAmp", &fovclustamp);
  fovclustmaxamp=new vector <double>;fouttree->Branch("ClustMaxAmp", &fovclustmaxamp);
  fovclustmaxstrip=new vector <double>;fouttree->Branch("ClustMaxStrip", &fovclustmaxstrip);
  fovclusttime=new vector <double>;fouttree->Branch("ClustTime", &fovclusttime);
  fovclustx=new vector <double>;fouttree->Branch("ClustX", &fovclustx);
  fovclustn=new vector <int>;fouttree->Branch("ClustN", &fovclustn);
  fovclustsize=new vector <int>;fouttree->Branch("ClustSize", &fovclustsize);


  //track related
  fovrayx=0; fouttree->Branch("RayX", &fovrayx);
  fovrayy=0; fouttree->Branch("RayY", &fovrayy);
  fovrayxup=0; fouttree->Branch("RayXup", &fovrayxup);
  fovrayyup=0; fouttree->Branch("RayYup", &fovrayyup);
  fovrayxdn=0; fouttree->Branch("RayXdn", &fovrayxdn);
  fovrayydn=0; fouttree->Branch("RayYdn", &fovrayydn);
  fovrayzdn=0; fouttree->Branch("RayZdn", &fovrayzdn);
  fovrayzup=0; fouttree->Branch("RayZup", &fovrayzup);
  fochi2x=0; fouttree->Branch("Chi2X", &fochi2x);
  fochi2y=0; fouttree->Branch("Chi2Y", &fochi2y);
}

void Draw(){//to display stuf at the end
  hRayXY->Draw("colz");
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


int getNdet(int idet){
  switch(idet)
    {
    case 0:
    case 1:
    case 2:
    case 3: return 1; break;//TPOT 01/2022
    case 4:
    case 5:
    case 6:
    case 7: return -1; break;
    case 8:
    case 9:
    case 10: 
    case 11: return 3; break; //RD3 D2 01/2022
    case 12:
    case 13:
    case 14:
    case 15: return 2; break; //RD3 D1 01/2022
    case 16:
    case 17:
    case 18:
    case 19: return -1; break;
    case 20:
    case 21:
    case 22:
    case 23: return 4; break; //RD3 D3 01/2022
    default : return -1;
    };
  return -1;
}

double getXdet(int idet,int istr){//this is so beautiful
  int det = getNdet(idet);
  int newdet =idet;
  switch(idet){//this part is for cable inversion
  case 0: newdet=1; break;
  case 1: newdet=2; break;
  case 2: newdet=3; break;
  case 3: newdet=4; break;
  case 22: newdet=22; break;
  case 23: newdet=23; break;
  }
  idet = newdet;
    
  //test
  if(1)//un inverted for 211130 // 2022 I AM NOW CERTAIN THAT TPOT IS ODD EVEN INVERTED
    if(det==1)
      {
	if(istr%2)//odd
	  istr--;
	else istr++;//even
      }
  switch(det)
    {
    case -1: return -1;
    case 1: //return (idet*fNstrip+(fNstrip-istr))*pitch; break;//inverted
            //return (idet*fNstrip+(istr))*pitch; break;//non inverted
      return 4*fNstrip-(idet*fNstrip+istr)*pitch; break;//full inverted
    case 2: return ((idet-12)*fNstrip+istr)*pitchrd3;break;
    case 3: return ((idet-8)*fNstrip+istr)*pitchrd3;break;
    case 4: return ((idet-20)*fNstrip+istr)*pitchrd3;break;
    default : return -1;
    };
  return -1;
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
  int ntr=0, ndettr=0; bool seen = false, trindet = false;;
  int ncycle=0, ievdata=0; int dist2look = 10;
  bool lostsync = false;
  int nerrcount = 0;
  int nev =  treData->GetEntries();
  int nevd =  treData->GetEntries();
  nev = (maxevn<nev && maxevn>0) ? maxevn:nev;
  nev = (nevd<nev) ? nevd:nev;
  int nevdatamax=nevd;
  
  cout << "Analysis loop : " << nev << endl;
  for(int iev=0; iev<nev;iev++){
    
    treData->GetEntry(iev);
    if(iev%1000 == 0) cout << "." << flush;
    if(DEBUGCL) cout << "EVENT " << iev << endl;	      	      
    //reset
    fovdetx->clear(); fovdetn->clear(); fovdetamp->clear(); fovdettime->clear();
    fovclustx->clear(); fovclustamp->clear();fovclustmaxamp->clear(); fovclusttime->clear(); fovclustmaxstrip->clear(); 
    fovclustsize->clear(); fovclustn->clear();
    double clustx = 0., clustamp=0., clustmaxamp=0., clustmaxstrip=0., clusttime=0.; int clustsize = 0, clustlaststrip = 0; int clustN=-1;
    fovdet->clear(); fovstrip->clear();fovsample->clear();
    foDataEvId=fevnData;

    for(int idet=0;idet<fNdet;idet++)//detloop
      for(int istr=0;istr<fNstrip;istr++)//max = fNdet
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
	  if(val>=nsigma*RMS[istr+fNstrip*idet] && val>30)//ZS
	    {
	      //ndettr++;	
	      fovstrip->push_back(istr);
	      fovdet->push_back(idet);
	      fovdetn->push_back(getNdet(idet));
	      fovsample->push_back(maxsample);
	      fovdetx->push_back(getXdet(idet,istr));
	      fovdetamp->push_back(val);
	      fovdettime->push_back(time);
	      //Clustering
	      if(getNdet(idet)!=clustN && clustsize)//not the same detector
		{
		  if(DEBUGCL) cout << "s[" << getNdet(idet) << "] [" << istr+idet*fNstrip << "] " << val
				   << "\t" << clustx/clustamp 
				   << "\t" << clustamp
				   << "\t" << clustmaxamp
				   << "\t" << clustsize << endl;
		  
		  fovclustx->push_back(clustx/clustamp); clustx=0.;
		  fovclusttime->push_back(clusttime/clustsize); clusttime=0.;//cluster time is average for now ...
		  fovclustamp->push_back(clustamp); clustamp=0.;
		  fovclustsize->push_back(clustsize); clustsize=0;
		  fovclustmaxamp->push_back(clustmaxamp); clustmaxamp=0.;
		  fovclustmaxstrip->push_back(clustmaxstrip); clustmaxstrip=-1;
		  fovclustn->push_back(clustN); clustN=-1;
		}
	      
	      //clustering
	      if(clustsize==0)//creating cluster
		{
		  clustamp=val;
		  clustmaxamp=val;
		  clustmaxstrip = getXdet(idet,istr);
		  clusttime=time;
		  clustsize=1;
		  clustx=getXdet(idet,istr)*val;
		  clustlaststrip = getXdet(idet,istr);
		  clustN=getNdet(idet);
		  if(DEBUGCL) cout << "\nc[" << getNdet(idet) << "] [" << istr+fNstrip*idet 
				   << "]\tgx=" << getXdet(idet,istr)
				   << "\ts=" << istr << " d=" << idet 
				   << "\t" << clustx/clustamp
				   << "\t" << val
				   << "\t" << clustamp
				   << "\t" << clustmaxamp
				   << "\t" << clustsize << endl;
		  
		}
	      else
		{
		  if(fabs(getXdet(idet,istr)-clustlaststrip)<4.5)//we are still in cluster, we tolerate one missing strip
		    {
		      clustamp+=val;
		      if(val>clustmaxamp) {clustmaxamp=val;clustmaxstrip=getXdet(idet,istr);}
		      clusttime+=time;
		      clustsize++;
		      clustx+=getXdet(idet,istr)*val;
		      clustlaststrip = getXdet(idet,istr);
		      if(DEBUGCL) cout << "+[" << getNdet(idet) << "] [" << istr+fNstrip*idet 
				       << "]\tgx=" << getXdet(idet, istr)
				       << "\ts=" << istr << " d=" << idet 
				       << "\t" << val 
				       << "\t" << clustx/clustamp
				       << "\t" << clustamp
				       << "\t" << clustmaxamp
				       << "\t" << clustsize << endl;
		      
		    }
		  else//closing cluster + create new
		    {
		      if(DEBUGCL) cout << "scl[" << getNdet(idet) << "] [" << istr+fNstrip*idet 
				       << "]\tgx=" << getXdet(idet, istr)
				       << "\ts=" << istr << " d=" << idet 
				       << "\t" << val 
				       << "\tspos = " << clustx/clustamp 
				       << "\t" << clustamp
				       << "\t" << clustmaxamp
				       << "\t" << clustsize << endl;
		      
		      fovclustx->push_back(clustx/clustamp); clustx=0.;
		      fovclusttime->push_back(clusttime/clustsize); clusttime=0.;//cluster time is average for now ...
		      fovclustamp->push_back(clustamp); clustamp=0.;
		      fovclustsize->push_back(clustsize); clustsize=0;
		      fovclustmaxamp->push_back(clustmaxamp); clustmaxamp=0.;
		      fovclustmaxstrip->push_back(clustmaxstrip); clustmaxstrip=-1.;
		      fovclustn->push_back(clustN); clustN=-1;
		      
			  //create new
		      clustamp=val;
		      clustmaxamp=val;
		      clustmaxstrip=getXdet(idet,istr);
		      clusttime=time;
		      clustsize=1;
		      clustx=getXdet(idet,istr)*val;
		      clustlaststrip = getXdet(idet,istr);
		      clustN=getNdet(idet);
		      if(DEBUGCL) cout << "c[" << getNdet(idet) << "] [" << istr+fNstrip*idet << "]\t" << val 
				       << "]\tgx=" << getXdet(idet, istr)
				       << "\ts=" << istr << " d=" << idet 
				       << "\t" << clustx/clustamp
				       << "\t" << clustamp
				       << "\t" << clustmaxamp
				       << "\t" << clustsize << endl;
		      
		    }
		}
	      //laststrip
	      if(istr==fNdet*fNstrip-1 && clustsize)//save the last cluster in
		{
		  if(DEBUGCL) cout << "ls[" << getNdet(idet) << "] [" << istr+fNstrip*idet << "]\t" << val 
				   << "\t" << clustx/clustamp
				   << "\t" << clustamp
				   << "\t" << clustmaxamp
				   << "\t" << clustsize << endl;
		  
		  fovclustx->push_back(clustx/clustamp); clustx=0.;
		  fovclusttime->push_back(clusttime/clustsize); clusttime=0.;//cluster time is average for now ...
		  fovclustamp->push_back(clustamp); clustamp=0.;
		  fovclustsize->push_back(clustsize); clustsize=0;
		  fovclustmaxamp->push_back(clustmaxamp); clustmaxamp=0.;
		  fovclustmaxstrip->push_back(clustmaxstrip); clustmaxstrip=-1.;
		  fovclustn->push_back(clustN); clustN=-1;		      
		}
	    }//Zero suppresion on val
	}//strip loop
    
    //we save the last cluster in memory
    if(clustsize)//not the same detector
      {
	if(DEBUGCL) cout << "e[" << -1 <<  "] [" << clustmaxstrip  
			 << "]\t" << clustx/clustamp 
			 << "\t" << clustamp
			 << "\t" << clustmaxamp
			 << "\t" << clustsize << endl;
	
	fovclustx->push_back(clustx/clustamp); clustx=0.;
	fovclusttime->push_back(clusttime/clustsize); clusttime=0.;//cluster time is average for now ...
	fovclustamp->push_back(clustamp); clustamp=0.;
	fovclustsize->push_back(clustsize); clustsize=0;
	fovclustmaxamp->push_back(clustmaxamp); clustmaxamp=0.;
	fovclustmaxstrip->push_back(clustmaxstrip); clustmaxstrip=-1;
	fovclustn->push_back(clustN); clustN=-1;
      }
	
    
    fouttree->Fill();ndettr++;		      
  }//iev loop
  fouttree->Write();
}//DoAnalysis


//---------------------------------------------------
//-----------MAIN------------------------------------
void TBanalysisRD3_NoRays(int maxevn = -1) {//todo args...
  InitDataTree("workdir/output.root");
  InitHisto();
  ReadRMS("workdir/RMSPed.dat",false);
  DoAnalysis(maxevn);
  // Draw();
}

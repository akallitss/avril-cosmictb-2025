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
const double zdet = 709.1;//TPOT
//const double zdet = 210;//RD3
const double nsigma = 4.5;
const double PITCHZ = 2.;
const double PITCHRPHI = 1.;
const int NHOLES = 0;//holes allowed in clustering
const int NDET = 4;
const double pitchrd3 = 0.5;
const double THETA = 0.84;//0.85;//degrees

//for efficency
const double xdet = -137.6;
const double xres = 4.5;//mm


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
vector <double> * fovdetrms;
vector <int> * fovdetch;
vector <int> * fovstrip;
vector <int> * fovdet;
vector <int> * fovdetn;
vector <int> * fovsample;
vector <double> * fovdetamp;
vector <double> * fovdettime;
vector <double> * fovclustamp,* fovclustmaxamp, *fovclustmaxstrip,*fovclustfact,  * fovclusttime,  * fovclustx;
vector <int>  * fovclustsize, *fovclustn;

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
  foutfile = new TFile("outtree.root","recreate");
  fouttree = new TTree("Tout","Event");
  fouttree->Branch("EvIdData", &foDataEvId); // event number
  fouttree->Branch("EvIdRay", &foRayEvId); // event number
  
  //detector hits 
  fovdetx=new vector <double>;fouttree->Branch("DetX", &fovdetx);
  fovdetrms=new vector <double>;fouttree->Branch("DetRMS", &fovdetrms);
  fovdetch=new vector <int>;fouttree->Branch("DetCh", &fovdetch);
  fovdetn=new vector <int>;fouttree->Branch("DetN", &fovdetn);
  fovstrip=new vector <int>; fouttree->Branch("StripNb", &fovstrip);
  fovdet=new vector <int>; fouttree->Branch("Det", &fovdet);
  fovsample=new vector <int>; fouttree->Branch("DetMaxSample", &fovsample);
  fovdetamp=new vector <double>;fouttree->Branch("DetAmp", &fovdetamp);
  fovdettime=new vector <double>;fouttree->Branch("DetTime", &fovdettime);
  //detector clusters
  fovclustamp=new vector <double>;fouttree->Branch("ClustAmp", &fovclustamp);
  fovclustmaxamp=new vector <double>;fouttree->Branch("ClustMaxAmp", &fovclustmaxamp);
  fovclustmaxstrip=new vector <double>;fouttree->Branch("ClustMaxStrip", &fovclustmaxstrip);
  fovclustfact=new vector <double>;fouttree->Branch("ClustFact", &fovclustfact);
  fovclusttime=new vector <double>;fouttree->Branch("ClustTime", &fovclusttime);
  fovclustx=new vector <double>;fouttree->Branch("ClustX", &fovclustx);
  fovclustn=new vector <int>;fouttree->Branch("ClustN", &fovclustn);
  fovclustsize=new vector <int>;fouttree->Branch("ClustSize", &fovclustsize);

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

int getNdet(int idet){
  switch(idet)
    {
    case 0://Feu 134 TPOT Z
    case 1:
    case 2:
    case 3: return 1; break;
    case 4://Feu 134 TPOT PHI
    case 5:
    case 6:
    case 7: return 2; break;
    case 8://Feu 40 RD3-D1 and D2
    case 9:
    case 10:
    case 11: return 3; break;
    case 12:
    case 13:
    case 14:
    case 15: return 4; break;
    case 16: //Feu 28 RD3-D3
    case 17:
    case 18:
    case 19: return 5; break;
    default : return -1;
    };
  return -1;
}

double getXdet(int idet,int istr){//this is so beautiful
  int det = getNdet(idet);
  int newdet =idet;
  int offsetdet=0; 
  //offset in dream number, first dream of a detector should be 0
  switch(det){
  case 1:  offsetdet = 0; break;
  case 2:  offsetdet = 4; break;
  }

  if(true)
    switch(idet){
    case 0: newdet=0; break;
    case 1: newdet=1; break;
    case 2: newdet=2; break;
    case 3: newdet=3; break;

    case 4: newdet=7; break;
    case 5: newdet=6; break;
    case 6: newdet=5; break;
    case 7: newdet=4; break;
    }
  idet = newdet;
// if(idet == 22) idet=23;
//   else if(idet == 23) idet=22;//swap cables ...
    
  //inversion piste a piste
  if(true)//un inverted for 211130
    if(det==1 || det==2)
      {
	if(istr%2)//odd
	  istr--;
	else istr++;//even
      }
  double srd3 = 4*fNstrip*pitchrd3;
  switch(det)
    {
    case -1: return -1;
    case 1: return ((idet-offsetdet)*fNstrip+(fNstrip-istr))*PITCHZ; break;//inverted
      //return ((idet-offsetdet)*fNstrip+(istr))*PITCHZ; break;//non inverted
      //case 2: return ((idet-offsetdet)*fNstrip+(fNstrip-istr))*PITCHRPHI;break;
    case 2: return ((idet-offsetdet)*fNstrip+(istr))*PITCHRPHI;break;
    case 3: return srd3-((idet-8)*fNstrip+istr)*pitchrd3;break;
    case 4: return srd3-((idet-16)*fNstrip+istr)*pitchrd3;break;
    default : return -1;
    };
  return -1;
}


//-------------------------------------------------
//------------------ACTUAL ANALYSIS PART-----------
void DoAnalysis(int maxevn){
  int ntr=0, ndettr=0; bool seen = false, trindet = false;;
  int ncycle=0, ievdata=0; int dist2look = 10;
  bool lostsync = false;
  int nerrcount = 0;
  int nevd = treData->GetEntries();
  int nevdatamax=nevd;
  int nev = nevd;
  cout << "Analysis loop : " << nev << endl;
  for(int iev=0; iev<nev;iev++){
    treData->GetEntry(iev);
    if(iev%1000 == 0) cout << "." << flush;
    if(DEBUGCL) cout << "\n------------\nEVENT " << iev << endl;
    fovdetx->clear();fovdetch->clear();fovdetrms->clear(); fovdetn->clear(); fovdetamp->clear(); fovdettime->clear();
    fovclustx->clear(); fovclustamp->clear();fovclustmaxamp->clear(); fovclusttime->clear(); fovclustmaxstrip->clear(); 
    fovclustsize->clear(); fovclustn->clear();fovclustfact->clear();
    //double clustx = 0., clustamp=0., clustmaxamp=0., clustmaxstrip=0., clusttime=0.; int clustsize = 0, clustlaststrip = 0; int clustN=-1;
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
	  if(val>=nsigma*RMS[istr+fNstrip*idet] && val>50 && RMS[istr+fNstrip*idet]<40 && RMS[istr+fNstrip*idet]>10 )//ZS + cleanup
	    {
	      fovstrip->push_back(istr);
	      fovdet->push_back(idet);
	      fovdetn->push_back(getNdet(idet));
	      fovsample->push_back(maxsample);
	      fovdetx->push_back(getXdet(idet,istr));//contient le mapping
	      fovdetch->push_back(idet*fNstrip+istr);//raw dream channel
	      fovdetrms->push_back(RMS[istr+fNstrip*idet]);//noise from pedestal run
	      fovdetamp->push_back(val);
	      fovdettime->push_back(time);
	    }//Zero suppresion on val
	}//strip loop
    //NOW CLUSTERING
    if(DEBUGCL) cout << "["<<fovdetx->size()<<"] hits : " <<endl; 
    map <double, int> mhitnb[NDET];
    for(int ihit = 0 ; ihit< fovdetx->size(); ihit++)
      mhitnb[fovdetn->at(ihit)-1][fovdetx->at(ihit)] = ihit;//building index to have hits sorted by position
    //detector loop
    for(int clidet=0;clidet<NDET;clidet++){
      //clustring variables
      bool inclust = false;
      double clustamp = 0., clustx=0., clusttime=0., clustmaxamp =0.,  clustmaxstrip=0.;
      int clustsize=0, lasthit=0; double localpitch=0.;
      //hit loop per detectctor
      if(clidet==0) localpitch=PITCHZ;
      else if(clidet==1) localpitch=PITCHRPHI;
      for( map<double,int>::iterator ii=mhitnb[clidet].begin(); ii!=mhitnb[clidet].end(); ++ii)
	{
	  int ihit = (*ii).second;
	  if(DEBUGCL) cout << clidet << "  -  ["<< ihit << "] (" << (*ii).first << ") {" << fovdetamp->at(ihit) << "}" <<endl;
	  if(!inclust)//not in cluster, we create
	    {
	      inclust = true;lasthit=ihit;
	      clustamp = fovdetamp->at(ihit); 
	      clustx= fovdetx->at(ihit)*clustamp;//if(DEBUGCL) cout << "detx*amp = " << clustx <<endl; 
	      clusttime=fovdettime->at(ihit);  
	      clustmaxamp = fovdetamp->at(ihit); 
	      clustmaxstrip= fovdetx->at(ihit); 
	      clustsize=1;
	    }
	  else{
	    if(fovdetx->at(ihit)-fovdetx->at(lasthit)<(1.5+NHOLES)*localpitch)//consecutive strip ?
	      {
		clustamp += fovdetamp->at(ihit); 
		clustx += fovdetx->at(ihit)*fovdetamp->at(ihit);
		clusttime += fovdettime->at(ihit);  
		if(clustmaxamp<fovdetamp->at(ihit)){ 
		  clustmaxamp = fovdetamp->at(ihit); 
		  clustmaxstrip= fovdetx->at(ihit); }
		//clustn=fovdetn->at(ihit);
		clustsize++;
		lasthit=ihit; 
	      }
	    else {//save cluster + new one
	      fovclustx->push_back(clustx/clustamp);
	      fovclusttime->push_back(clusttime/clustsize);
	      fovclustamp->push_back(clustamp);
	      fovclustsize->push_back(clustsize);
	      fovclustmaxamp->push_back(clustmaxamp);
	      fovclustmaxstrip->push_back(clustmaxstrip);
	      fovclustn->push_back(clidet+1);
	      inclust = true;lasthit=ihit;
	      clustamp = fovdetamp->at(ihit); 
	      clustx= fovdetx->at(ihit)*clustamp;
	      clusttime=fovdettime->at(ihit);  
	      clustmaxamp = fovdetamp->at(ihit); 
	      clustmaxstrip= fovdetx->at(ihit);
	      //clustn=fovdetn->at(ihit); 
	      clustsize=1;
	    }
	  }//inclust = true
	}//for hit
      if(inclust){//last cluster needs to saved
	fovclustx->push_back(clustx/clustamp); //if(DEBUGCL) cout << clustx << "/" << clustamp << clustx/clustamp <<endl; 
	fovclusttime->push_back(clusttime/clustsize);
	fovclustamp->push_back(clustamp);
	fovclustsize->push_back(clustsize);
	fovclustmaxamp->push_back(clustmaxamp);
	fovclustmaxstrip->push_back(clustmaxstrip);
	fovclustn->push_back(clidet+1);
      }
      if(DEBUGCL){
	cout << clidet << "--> [" << fovclustx->size() << "] cluster(s) : "<< endl;
	for(int icl = 0; icl<fovclustx->size(); icl++)
	  if(fovclustn->at(icl)==clidet+1)  cout << "[" << icl << "] (" << fovclustx->at(icl)<< ") {" << fovclustamp->at(icl) << "}" << endl;
      }
    }	
    fouttree->Fill();//ndettr++;		      
  }//iev loop
  fouttree->Write();
  //cout << "EFF = " << ndettr << "/" << ntr << " = " << 1.*ndettr/ntr << endl;  
  
}//DoAnalysis


//---------------------------------------------------
//-----------MAIN------------------------------------
void TBanalysisTPOT_2022_Fe55(string filedata, string filerms, int maxevn = -1) {//todo args...
  InitDataTree(filedata);
  InitHisto();
  ReadRMS(filerms,false);
  DoAnalysis(maxevn);
  // Draw();//do spectrum



}

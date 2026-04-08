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
const double zdet = 225+20;//Plateau du bas + spacer pour p2
const double nsigma = 4.5;
const double MINAMP = 20; //adc to have a signal regardless of noise 

//alignament
const double THETA = 0; 

//for efficency
const double xdet = -137.6;
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
TFile * foutfile = 0;
TTree * fouttree = 0;
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

int InitRayTreeRead(TString filename){
  cout << "Read Ray tree : " << filename << endl; 
  fray = new TFile(filename.Data());
  treRay =  (TTree*)fray->Get("T");  
  cout << "Ray Tree = " << treRay->GetName() << "\t"<< treRay->GetEntries() << " events" << endl;
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

  switch(idet){
  case 0: newdet=3; break;
  case 1: newdet=2; break;
  case 2: newdet=1; break;
  case 3: newdet=0; break;
    
  case 4: newdet=7; break;
  case 5: newdet=6; break;
  case 6: newdet=5; break;
    case 7: newdet=4; break;
  }
  idet = newdet;
  return newdet*idet+istr; //not working for P2 ; to be done
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
  int nev = treRay->GetEntries();
  int nevd = treData->GetEntries();
  nev = (maxevn<nev && maxevn>0) ? maxevn:nev;
  nev = (nevd<nev) ? nevd:nev;
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
    treData->GetEntry(evn-1);//iev
    if(iev%1000 == 0) cout << "." << flush;
    if(DEBUGCL) cout << "\n------------\nEVENT " << iev << endl;
    if(turniev>=4095) {turn++; turniev=0;} else turniev++;
    int realevndata = fevnData ;//+ turn*4096 ;
    if(evn%4096 != realevndata) {
      if(skipped<100) cout << iev << " " << evn << " " << realevndata << endl;
      if(skipped<100) cout <<  evn%4096 << " " << realevndata << endl;
      skipped++; continue;}
        
    for(int iray=0;iray<rayN;iray++)
      {	
	fovdetx->clear();fovdetch->clear();fovdetrms->clear(); fovdetn->clear(); fovdetamp->clear(); fovdettime->clear();
	fovclustx->clear(); fovclustamp->clear();fovclustmaxamp->clear(); fovclusttime->clear(); fovclustmaxstrip->clear(); 
	fovclustsize->clear(); fovclustn->clear();fovclustfact->clear();
	//double clustx = 0., clustamp=0., clustmaxamp=0., clustmaxstrip=0., clusttime=0.; int clustsize = 0, clustlaststrip = 0; int clustN=-1;
	fovdet->clear(); fovstrip->clear();fovsample->clear();
	foDataEvId=fevnData;
	foRayEvId=evn;

	double rayx = getXray(zdet,iray);
	double rayy = getYray(zdet,iray);
	
	//rotation of tracker
	double theta = THETA*TMath::Pi()/180.;//in radian
	double rayxprime = rayx*cos(theta)-rayy*sin(theta);
	double rayyprime = rayx*sin(theta)+rayy*cos(theta);
	vector <double> cla, cls;	  
	
	rayx=rayxprime; rayy=rayyprime;

	hRayXY->Fill(rayx, rayy);
	// fovrayx->push_back(getXray(zdet,iray));//find a way to get a bijection; only one track??
	// fovrayy->push_back(getYray(zdet,iray));
	fovrayx = rayx;
	fovrayy = rayy;//X_Up->at(nray)-X_Down->at(nray))/(Z_Up-Z_Down);
	fovrayxup = X_Up->at(iray);	fovrayxdn = X_Down->at(iray);
	fovrayyup = Y_Up->at(iray);	fovrayydn = Y_Down->at(iray);
	fovrayzup = Z_Up;	fovrayzdn = Z_Down;
	fochi2x = Chi2X->at(iray);
	fochi2y = Chi2Y->at(iray);
	seen =false;
	if(fovrayy>-100 && fovrayy<100 && fovrayx>120 && fovrayx<220  && fochi2x<20 && fochi2y<20) // effiency test ok
	  {trindet=true; ntr++;}
	else trindet=false;
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
	      if(val>=nsigma*RMS[istr+fNstrip*idet] && val>MINAMP && RMS[istr+fNstrip*idet]>4 )//ZS + cleanup
		{
		  if(iev<10 && getXdet(idet,istr)==165) cout << "idet = " << idet << ", istr =" << istr << " ,rms=" << RMS[istr+fNstrip*idet] << endl;
		  //ndettr++;		      
		  //efficiency and hit stuff
		  if(val>MINAMP && 
		     getXdet(idet,istr)-fovrayx-xdet>-4.*xres && getXdet(idet,istr)-fovrayx-xdet<4.*xres && 
		     trindet && !seen) //track in detector and not seen yet
		    {
		      seen=true; 
		      ndettr++;		      
		      hRayXYseen->Fill(fovrayx,fovrayy);
		    }//efficiency cut
		  else hRayXYseen->Fill(fovrayx,fovrayy, 0.);
		  if(seen)
		    {
		      hRayXDetX->Fill(fovrayx,getXdet(idet,istr));
		      hRayYDetX->Fill(fovrayy,getXdet(idet,istr));
		    }
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
	fouttree->Fill();//ndettr++;		      
      }//nray
  }//iev loop
  fouttree->Write();
  cout << endl <<nev << "events, " << skipped << "skipped." << endl;
  //cout << "EFF = " << ndettr << "/" << ntr << " = " << 1.*ndettr/ntr << endl;  
  
}//DoAnalysis


//---------------------------------------------------
//-----------MAIN------------------------------------
void TBanalysisP2_2024(int maxevn = -1) {//todo args...
  InitRayTreeRead("workdir/run_rays.root");
  InitDataTree("workdir/output.root");
  InitHisto();
  ReadRMS("workdir/RMSPed.dat",false);
  DoAnalysis(maxevn);
  // Draw();
}

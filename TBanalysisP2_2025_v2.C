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
#include <map>

using namespace std;

#ifdef __MAKECINT__
#pragma link C++ class vector<double>+;
#endif

//DEBUG
bool DEBUGCL = false; //printout clustering

//Analysis constant
const double zdet = 225+20;//Plateau du bas + spacer pour p2
const double nsigma = 4.5;
const double MINAMP = 10; //adc to have a signal regardless of noise 

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
float *datatab;
//float data[16][64][32];

//RMS:
float *RMS;

//Out Tree
TFile * foutfile = 0;
TTree * fouttree = 0;
vector <double> * fovdetx;
vector <double> * fovdety;
vector <double> * fovdetr;
vector <double> * fovdetphi;
vector <double> * fovdetrms;
vector <int> * fovdetch;
vector <int> * fovstrip;
vector <int> * fovdet;
vector <int> * fovdetn;
vector <int> * fovsample;
vector <double> * fovdetamp;
vector <double> * fovdettime;
vector <double> * fovclustamp,* fovclustmaxamp, *fovclustmaxstrip,*fovclustfact,  * fovclusttime,  * fovclustx, *fovclusty, *fovclustmaxstripx, *fovclustmaxstripy ;
vector <int>  * fovclustsize, *fovclustn;
double fovlargestclust;


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

//struct maaping p2
struct mp2
{
  int conn;
  int chan;
  float X;
  float Y;
  float R;
  float Phi;
};
map<int, mp2> vmp2; //vector containing all mapping

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

void ReadMapping(string folderMapping){

  int Connector; string  PinName;int   Channel; int  PadIndex;string   PadName; float X,Y,Radius,Phi,DeltaPhi;
  ifstream inrms; string line;
  
  for(int b=0; b<10; b++)
    {
      char connf[20]; sprintf(connf, "/connector_%d_v2.txt",b);
      string fileMapping = connf;
      string fileconn = folderMapping + fileMapping;
      cout << "open mapping file : " << fileconn << endl; 
      inrms.open(fileconn.c_str());
      if (inrms.is_open()) 
	while (std::getline(inrms, line)) {
	  //std::cout << line << std::endl;
	  int nv = sscanf(line.c_str(),"%d %*s %d %d %*s %f %f %f %f %f",&Connector,&Channel,&PadIndex,&X, &Y,&Radius, &Phi, & DeltaPhi);
	  //cout << "["<< nv<<"]"<< Connector << " " << X << endl;
    if (Channel < 64){ Channel = 63 - Channel; }
    else{ Channel = 191 - Channel;}
	  if(nv>=7) {
	    mp2 currpad; 
	    currpad.conn = Connector;
	    currpad.chan = 128*Connector + Channel;
	    currpad.X = X;
	    currpad.Y = Y;
	    currpad.R = Radius;
	    currpad.Phi = Phi;
	    vmp2[currpad.chan] = currpad;
      // cout << "dans lecture mapping : connector : " << Connector << " chn : " << currpad.chan << " X, Y : " << X << " " << Y << endl;
	  }
	}
      else cout << "failed to open mapping file : " << fileconn << endl; 
      inrms.close();
    }//for
    std::cout << "finish reading mapping" << std::endl;
}//readmppaing


void InitHisto(){

  //ttree
  TFile * foutfile = new TFile("outtree.root","recreate");
  fouttree = new TTree("Tout","Event");
  fouttree->Branch("EvIdData", &foDataEvId); // event number
  fouttree->Branch("EvIdRay", &foRayEvId); // event number
  
  //detector hits 
  fovdetx=new vector <double>;fouttree->Branch("DetX", &fovdetx);
  fovdety=new vector <double>;fouttree->Branch("DetY", &fovdety);
  fovdetr=new vector <double>;fouttree->Branch("DetR", &fovdetr);
  fovdetphi=new vector <double>;fouttree->Branch("DetPhi", &fovdetphi);
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
  fovclustmaxstripx=new vector <double>;fouttree->Branch("ClustMaxStripX", &fovclustmaxstripx);
  fovclustmaxstripy=new vector <double>;fouttree->Branch("ClustMaxStripY", &fovclustmaxstripy);
  fovclustfact=new vector <double>;fouttree->Branch("ClustFact", &fovclustfact);
  fovclusttime=new vector <double>;fouttree->Branch("ClustTime", &fovclusttime);
  fovclustx=new vector <double>;fouttree->Branch("ClustX", &fovclustx);
  fovclusty=new vector <double>;fouttree->Branch("ClustY", &fovclusty);
  fovclustn=new vector <int>;fouttree->Branch("ClustN", &fovclustn);
  fovclustsize=new vector <int>;fouttree->Branch("ClustSize", &fovclustsize);
  fouttree->Branch("LargestClustAmp", &fovlargestclust);


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
  datatab = new float[fNdet*fNstrip*fNsample]();
 
  treData->SetBranchAddress("StripAmpl_corrped",datatab); 
  
  int iev=1;
  treData->GetEntry(iev);//[3][2][1]
  cout << "DATA EVENT "<<iev<<" (" << fevnData << "): data[0][0][0]= " << datatab[1*fNsample+2*(fNstrip+3*fNdet)] <<endl;
  for(int i=0;i<10;i++)
    cout << datatab[i+fNsample*(22+4*fNstrip)] << " ";
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

double getYray(double z, int nray){
  //solve Y=Az+B
  if(Z_Up-Z_Down == 0) {cout << "ERROR: Z_up == Z_Down" << endl; return -1;}  
  double A=(Y_Up->at(nray)-Y_Down->at(nray))/(Z_Up-Z_Down);
  double B=Y_Up->at(nray)-A*Z_Up;
  double y=A*z+B;//tadaaaa
  return y;
}


int getNdet(int idet){
  switch(idet)
    {
    case 0:return 0; break;
    case 1:return 1; break;
    case 2:return 2; break;
    case 3:return 3; break;
    case 4:return 4; break;
    case 5:return 5; break;
    case 6:return 6; break;
    case 7:return 7; break;
    case 8:return 8; break;
    case 9:return 9; break;
    case 10:return 10; break;
    case 11:return 11; break;
    case 12:return 12; break;
    case 13:return 13; break;
    case 14:return 14; break;
    case 15:return 15; break;
    case 16:return 16; break;
    case 17:return 17; break;
    case 18:return 18; break;
    case 19: return 19; break;;
    default : return -1;
    };
  return -1;
}

double getChan(int idet,int istr){//this is so beautiful
  int newdet = getNdet(idet); if(newdet<0) return -1; 
  int chan = newdet*fNstrip+istr;
  return chan;
}

void cluster(vector<int>* detchn, vector<double>* detx, vector<double>* dety, vector<double>* detamp, 
vector<double>* clusterAmplitudes,  vector<int>* clusterSize, vector<double>* clusterMaxAmp, vector<double>* clusterMaxStrip, 
vector<double>* clusterX, vector<double>* clusterY,vector<double>* clusterMaxStripX, vector<double>* clusterMaxStripY ){
  int n = detchn->size();
  vector<bool> visited(n,false);

  for(int i=0;i<n;i++){
    if(visited[i]) continue;

    double clusterAmp = 0.0;
    int Size = 0;
    double clustx = 0;
    double clusty = 0;
    vector<double> amps;
    vector<double> pads;
    vector<double> x;
    vector<double> y;
    vector<int> stack;
    stack.push_back(i);
    visited[i] = true;

    while(!stack.empty()){
      int k = stack.back();
      stack.pop_back();

      
      clusterAmp += (*detamp)[k];
      clustx += (*detx)[k]*(*detamp)[k];
      clusty += (*dety)[k]*(*detamp)[k];
      Size++;
      amps.push_back((*detamp)[k]); 
      pads.push_back((*detchn)[k]); 
      x.push_back((*detx)[k]);
      y.push_back((*dety)[k]);
      

      for(int j=0;j<n;j++){
        if(visited[j]) continue;

        double dx = (*detx)[k] - (*detx)[j];  
        double dy = (*dety)[k] - (*dety)[j];
        double dist = sqrt(dx*dx + dy*dy);

        if(dist < 14.15){
          visited[j] = true;
          stack.push_back(j);
          }
      }
    }

    if(Size > 0){

      double max = amps[0];
      int idx = 0;

      for(int l=1;l<amps.size();l++){
        if(amps[l] > max){
          max = amps[l];
          idx = l;
        }
      }

      clusterMaxAmp->push_back(max);
      clusterMaxStrip->push_back(pads[idx]);
      clusterMaxStripX->push_back(x[idx]);
      clusterMaxStripY->push_back(y[idx]);
      clusterAmplitudes->push_back(clusterAmp);
      clusterSize->push_back(Size);
      clusterX->push_back(clustx/clusterAmp);
      clusterY->push_back(clusty/clusterAmp);
    }
  }
}

//-------------------------------------------------
//------------------ACTUAL ANALYSIS PART-----------
void DoAnalysis(int maxevn){


  // ====== Lecture des EVN dans treRay ======
int evn_ray;
treRay->SetBranchAddress("evn", &evn_ray);

Long64_t nRay = treRay->GetEntries();
vector<int> rayID(nRay);

for (Long64_t i=0; i<nRay; i++) {
    treRay->GetEntry(i);
    rayID[i] = evn_ray;
}

// ====== Lecture des IDEvent dans treData ======
int evn_data;
treData->SetBranchAddress("IDEvent", &evn_data);

Long64_t nData = treData->GetEntries();

std::vector<int> dataID_abs(nData);

int last = -1;
int wrap = 0;

for (Long64_t i = 0; i < nData; i++) {
    treData->GetEntry(i);

    if (last != -1 && evn_data < last) {
        
        wrap++;
    }
    last = evn_data;

    dataID_abs[i] = evn_data + wrap * 4096;
}
// ====== Construction des maps ======
map<int,int> indexRay, indexData;

for (Long64_t i=0; i<nRay; i++)
    indexRay[ rayID[i] ] = i;

for (Long64_t i=0; i<nData; i++)
    indexData[ dataID_abs[i] ] = i;

cout << " construction des maps fini" << endl;

  int ntr=0, ndettr=0; bool seen = false, trindet = false;;
  int nev = treRay->GetEntries();
  int nevd = treData->GetEntries();
  nev = (maxevn<nev && maxevn>0) ? maxevn:nev;
  nev = (nevd<nev) ? nevd:nev;
  // cout << "building index table for tree synchronisation... "  << endl;
  // int turn = 0; int turniev = 0;
  // if(0)
  // for(int iev=0; iev<nev;iev++){
  //   treRay->GetEntry(iev);//nev
  //   treData->GetEntry(iev);
  //   if(turniev>=4095) {turn++; turniev=0;} else turniev++;
  //   int realevndata = fevnData + turn*4096 ;//
  //   cout << "[" << iev << "] = " << evn << " == " <<  realevndata << "\t" <<  fevnData << endl ;
  // }
  // turn = 0;turniev = 0; int skipped = 0; int desync = 0;
  // cout << "Analysis loop : " << nev << endl;
  // for(int iev=0; iev<nev;iev++){
  //   treRay->GetEntry(iev);
  //   treData->GetEntry(evn-1);//iev
  //   if(iev%1000 == 0) cout << "." << flush;
  //   if(DEBUGCL) cout << "\n------------\nEVENT " << iev << endl;
  //   if(turniev>=4095) {turn++; turniev=0;} else turniev++;
  //   int realevndata = fevnData ;//+ turn*4096 ;
    
  //   if(evn%4096 != realevndata) {
  //       bool recovered = false;
  //       desync++;
        
  //       for (int k = 1; k < 1000; k++) {
  //         int testIndex = iev + k;
  //         if (testIndex >= treRay->GetEntries()) break;

  //         treRay->GetEntry(testIndex);
  //         if (desync < 3) cout << evn % 4096 << " " << realevndata << endl;
  //         if ((evn % 4096) == realevndata) {
  //             iev = testIndex;   // resynchronisation
  //             recovered = true;
              
  //             break;
  //         }
  //     }

  //     if (!recovered) {
  //       skipped++;
  //       if(skipped<100) cout << iev << " " << evn << " " << realevndata << endl;
  //       if(skipped<100) cout <<  evn%4096 << " " << realevndata << endl;
  //       continue;
  //     }
  //     // if(skipped<100) cout << iev << " " << evn << " " << realevndata << endl;
  //     // if(skipped<100) cout <<  evn%4096 << " " << realevndata << endl;
  //     // skipped++; continue;
  //   }
  int skipped =0;
  int countevn = 0;
int errordet=0;
  map<int,int>::iterator it;

  for (it = indexRay.begin(); it != indexRay.end(); ++it) {

    int ev = it->first;
    if(ev%1000 == 0) cout << "." << flush;
    if (indexData.count(ev) == 0)
        {skipped++;continue;}

    Long64_t idxRay  = it->second;
    Long64_t idxData = indexData[ ev ];

    treRay->GetEntry(idxRay);
    treData->GetEntry(idxData);

    for(int iray=0;iray<rayN;iray++)
      {	
fovdetx->clear(); fovdety->clear();
  fovdetr->clear();fovdetphi->clear();
	fovdetch->clear();fovdetrms->clear(); fovdetn->clear(); fovdetamp->clear(); fovdettime->clear();
  fovclustx->clear();fovclusty->clear(); fovclustamp->clear();fovclustmaxamp->clear(); fovclusttime->clear(); fovclustmaxstrip->clear(); fovclustmaxstripx->clear(); fovclustmaxstripy->clear();
	fovclustsize->clear(); fovclustn->clear();fovclustfact->clear();
	fovdet->clear(); fovstrip->clear();fovsample->clear();
	foDataEvId=fevnData;
	foRayEvId=evn;

	double rayx = getXray(zdet,iray);
	double rayy = getYray(zdet,iray);
	
	//rotation of tracker
	double theta = THETA*TMath::Pi()/180.;//in radian
	double rayxprime = rayx*cos(theta)-rayy*sin(theta) ;
	double rayyprime = rayx*sin(theta)+rayy*cos(theta) ;
	vector <double> cla, cls;	  
	
	rayx=rayxprime; rayy=rayyprime;

	fovrayx = rayx;
	fovrayy = rayy;//X_Up->at(nray)-X_Down->at(nray))/(Z_Up-Z_Down);
	fovrayxup = X_Up->at(iray);	fovrayxdn = X_Down->at(iray);
	fovrayyup = Y_Up->at(iray);	fovrayydn = Y_Down->at(iray);
	fovrayzup = Z_Up;	fovrayzdn = Z_Down;
	fochi2x = Chi2X->at(iray);
	fochi2y = Chi2Y->at(iray);
	seen =false;
	
	for(int idet=0;idet<fNdet;idet++){//detloop
	if (ev==1) cout << "idet pour ev=1 : "<< idet << endl;
	  for(int istr=0;istr<fNstrip;istr++)//max = fNdet
	    {
	      double val=0; double maxsample = 0; double time = -1;       
	      for(int is=0;is<fNsample;is++)
		{
		  if(datatab[is+fNsample*(istr+fNstrip*idet)]>val) {
		    val = datatab[is+fNsample*(istr+fNstrip*idet)];//compute max
		    maxsample = is;
		    //compute time at maximum
		    if(is>0 && is<fNsample-1)
		      {
			//Fit
			double atemp = 0.5*(datatab[is+1+fNsample*(istr+fNstrip*idet)]-2.*datatab[is+fNsample*(istr+fNstrip*idet)]
			 		    + datatab[is-1+fNsample*(istr+fNstrip*idet)]);
			double btemp = (datatab[is+fNsample*(istr+fNstrip*idet)]-datatab[is-1+fNsample*(istr+fNstrip*idet)])-atemp*(2.*is-1.);
			time = -0.5*btemp/atemp;
			//time += fFineTimeStamp[idet/8]*1./6.;//Coorection avec le Fine Time Stamp
		      }
		    else time=-1;
		  }
		}
	      if(val>=nsigma*RMS[istr+fNstrip*idet] && val>MINAMP && RMS[istr+fNstrip*idet]>1 )//ZS + cleanup
		{
			if(ev==1 && idet ==19)cout<< " strip " << istr << " val = "<< val << endl;
		  fovstrip->push_back(istr);
		  fovdet->push_back(idet);
		  fovdetn->push_back(getNdet(idet));
			if ((idet==19 || idet == 18) && getNdet(idet)==-1){errordet++;}
		  fovsample->push_back(maxsample);
		  fovdetrms->push_back(RMS[istr+fNstrip*idet]);//noise from pedestal run
		  fovdetamp->push_back(val);
		  fovdettime->push_back(time);
		  //requires mapping
      //cout << "idet : "<<idet << endl;
		  int channel = getChan(idet,istr);
		  fovdetch->push_back(channel);//contains connector mapping
		  if(channel>=0 && channel<=vmp2.size())//channel seems ok
		    {
		      fovdetx->push_back(vmp2[channel].X);//contient le mapping
		      fovdety->push_back(vmp2[channel].Y);//contient le mapping
          // if (idet == 0){
          // cout << "chn = " << channel <<  " vs chn dans vmp2 = "<< vmp2[channel].chan << endl;
          // cout << "X = " << vmp2[channel].X << " Y = "<< vmp2[channel].Y << endl;
          // }

		      fovdetr->push_back(vmp2[channel].R);//contient le mapping
		      fovdetphi->push_back(vmp2[channel].Phi);//contient le mapping
		    }
		}//Zero suppresion on val
	    }//strip loop
		}//detloop
	fouttree->Fill();//ndettr++;		      
      }//nray
  }//iev loop
  fouttree->Write();
  cout << endl <<nev << "events, " << skipped << "skipped." << endl;
cout << " error det num : " << errordet << endl;
  //cout << "EFF = " << ndettr << "/" << ntr << " = " << 1.*ndettr/ntr << endl;  
  
}//DoAnalysis


//---------------------------------------------------
//-----------MAIN------------------------------------
void TBanalysisP2_2025_v2(string folder="", int maxevn = -1) {
  string filer = folder + "run_rays.root";
  cout << "opening : " << filer << endl;
  InitRayTreeRead(filer.data());
  filer = folder + "output.root";
  cout << "opening : " << filer << endl;
  InitDataTree(filer.data());
  InitHisto();
  filer = folder + "RMSPed.dat";
  cout << "opening : " << filer << endl;
  ReadRMS(filer.data(),false);
  ReadMapping("mappingP2_v2");
  if(vmp2.size())
    cout << "Mapping for " << vmp2.size() << " found. [" << vmp2[0].conn << "]=" << vmp2[0].X << endl;
  DoAnalysis(maxevn);
  // Draw();
}

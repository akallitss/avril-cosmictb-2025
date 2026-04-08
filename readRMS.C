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

void readRMS(){
    //RMS:
    float *RMS;

    //Out Tree
    TFile * foutfile = 0;
    TTree * fouttree = 0;
    
    string rmsfile = "";
    ReadRMS(rmsfile, true);
    fouttree->Fill();		      
    fouttree->Write();
}
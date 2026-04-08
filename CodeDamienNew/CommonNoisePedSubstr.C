#include "Riostream.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TFile.h"
#include "TBranch.h"

#include<set>

int CommonNoisePedSubstr(string filename) {
  gROOT->Reset();

  //const char* fname = gROOT->GetListOfFiles()->First()->GetName();  
  TFile f(filename.c_str(), "update");
  cout<<" root file to update: "<<filename<<endl;
  //TFile f(fname, "update");
  TTree *T = (TTree*)f.Get("T");  
  
  //get data size
  int Ndet=0;
  int Nsample=0;
  int Nstrip=0;
  const char * branchname =  T->GetBranch("StripAmpl")->GetTitle();  
  sscanf(branchname,"StripAmpl[%d][%d][%d]/F",&Ndet,&Nstrip,&Nsample);
  cout << "Data have " << Ndet << " DREAM, " << Nstrip << " strips, " 
       << Nsample << " samples." << endl;
  
  ifstream in;
  in.open("Pedestal.dat");
  int nlines=0;
  int det, strip;
  float pedestal_value;
  float Pedestal[Ndet][Nstrip];
  while (1) { // read the text file
    in >> det >> strip >> pedestal_value;
    Pedestal[det][strip] = pedestal_value;
    if (!in.good()) break;
    nlines++;
  }



//Declaration of leaves types
  Int_t           Nevent;
  int IDEvent = 0;  // ID of event in Feu header
  Int_t           StripAmpl[Ndet][Nstrip][Nsample];
  Int_t           TsampleNum[Nsample];
  Float_t StripAmpl_corrped[Ndet][Nstrip][Nsample];
  // Set branch addresses.
  T->SetBranchAddress("Nevent",&Nevent);
  T->SetBranchAddress("IDEvent", &IDEvent); // ID of the event
  T->SetBranchAddress("StripAmpl",StripAmpl);
  T->SetBranchAddress("TsampleNum",TsampleNum);
  char leefStripAmpl[100];
  sprintf(leefStripAmpl,"StripAmpl_corrped[%d][%d][%d]/F",Ndet,Nstrip,Nsample);
  TBranch *newBranch = T->Branch("StripAmpl_corrped", StripAmpl_corrped,leefStripAmpl);

  Long64_t nentries = T->GetEntries();
  Long64_t nbytes = 0;
  set<double> sStripAmplSort1, sStripAmplSort2;  // to separate thin and large strips
  int stripmin=0, stripmax=64, stripmid=32;
  bool splittedStrips = false;

  for (Long64_t i=0; i<nentries; i++) {
    nbytes += T->GetEntry(i);

    for ( int idet=0; idet<Ndet; idet++) {
      for ( int isample=0; isample<Nsample; isample++) {
	sStripAmplSort1.clear();
	sStripAmplSort2.clear();
	if (splittedStrips) {
	  for ( int istrip=stripmin; istrip<stripmid; istrip++)
	    sStripAmplSort1.insert(StripAmpl[idet][istrip][isample] - Pedestal[idet][istrip]);
	  for ( int istrip=stripmid; istrip<stripmax; istrip++)
	    sStripAmplSort2.insert(StripAmpl[idet][istrip][isample] - Pedestal[idet][istrip]);
	} else {
	  for ( int istrip=stripmin; istrip<stripmax; istrip++)
	    sStripAmplSort1.insert(StripAmpl[idet][istrip][isample] - Pedestal[idet][istrip]);
	}
	 int icount = sStripAmplSort1.size()/2;
	 double medAmpl1 = 0, medAmpl2 = 0;
	for (set<double>::const_iterator isci = sStripAmplSort1.begin(); isci != sStripAmplSort1.end(); isci++) {
	  if (icount == 0) {
	    medAmpl1 = *isci;
	    break;
              }
	  icount--;
	}
	if (splittedStrips) {
	  icount = sStripAmplSort2.size()/2;
	  for (set<double>::const_iterator isci = sStripAmplSort2.begin(); isci != sStripAmplSort2.end(); isci++) {
	    if (icount == 0) {
	      medAmpl2 = *isci;
	      break;
	    }
	    icount--;
	  }
	}
	
            if (splittedStrips) {
	      for ( int istrip=stripmin; istrip<stripmid; istrip++)  {
	        StripAmpl_corrped[idet][istrip][isample] = StripAmpl[idet][istrip][isample] - Pedestal[idet][istrip] - medAmpl1;
	      }
	      for ( int istrip=stripmid; istrip<stripmax; istrip++)  {
	        StripAmpl_corrped[idet][istrip][isample] = StripAmpl[idet][istrip][isample] - Pedestal[idet][istrip] - medAmpl2;
	      }
            } else {
	      for (int istrip=stripmin; istrip<stripmax; istrip++)  {
	        StripAmpl_corrped[idet][istrip][isample] = StripAmpl[idet][istrip][isample] - Pedestal[idet][istrip] - medAmpl1;
	      }
            }
	}
    }
    newBranch->Fill();
  }
  T->Write("", TObject::kOverwrite);
  return 1;
}

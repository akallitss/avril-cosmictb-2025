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

int CommonNoisePedSubstrZS(string filename, double nsigma = 0) {
  cout << "opening " << filename << endl; 
  cout << "Performing ZS with a threshold of " << nsigma << " sigma RMS." << endl;
  gROOT->Reset();
  
  TFile f(filename.c_str(), "update");
  cout<<" root file to update: "<<filename<<endl;
  TTree *T = (TTree*)f.Get("T");  
  
  //get data size
  int Ndet=0;
  int Nsample=0;
  int Nstrip=0;
  int Nfeu=0;
  const char * branchname =  T->GetBranch("StripAmpl")->GetTitle();  
  sscanf(branchname,"StripAmpl[%d][%d][%d]/F",&Ndet,&Nstrip,&Nsample);
  cout << "Data have " << Ndet << " DREAM, " << Nstrip << " strips, " 
       << Nsample << " samples." << endl;
  Nfeu=Ndet/8;
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
  
  ifstream inrms;  
  nlines=0;
  float rms_value;
  float RMS[Ndet][Nstrip];
  if(nsigma>0){//zero suppression activated
    inrms.open("RMSPed.dat");
    while (1) { // read the text file
      inrms >> det >> strip >> rms_value;
      //cout << det << " " <<  strip << " " <<  rms_value << endl;
      RMS[det][strip] = rms_value;
      if (!inrms.good()) break;
      nlines++;
    }
    }
  else 
    for(int idet=0;idet<Ndet;idet++)
      for(int istr=0;idet<Nstrip;istr++)
	RMS[idet][istr] = 0;
  const char* fname = gROOT->GetListOfFiles()->First()->GetName();
  
  //Declaration of leaves types
  Int_t           Nevent;
  int IDEvent = 0;  // ID of event in Feu header
  Int_t           FineTimeStamp[Nfeu];
  Int_t           StripAmpl[Ndet][Nstrip][Nsample];
  Int_t           TsampleNum[Nsample];
  Int_t           StripNum[Nstrip];
  // Set branch addresses.
  T->SetBranchAddress("Nevent",&Nevent);
  T->SetBranchAddress("IDEvent", &IDEvent); // ID of the event
  T->SetBranchAddress("StripAmpl",StripAmpl);
  T->SetBranchAddress("TsampleNum",TsampleNum);
  T->SetBranchAddress("StripNum", StripNum);
  T->SetBranchAddress("FineTimeStamp", FineTimeStamp);
  //signal tree
  Int_t Neventsig=0;
  Int_t IDEventsig=0;  // ID of event in Feu header
  Int_t TsampleNumsig[Nsample];
  Int_t FineTimeStampsig[Nfeu];
  Float_t StripAmpl_corrped[Ndet][Nstrip][Nsample];
  Int_t StripNumsig[Nstrip];

  TFile *fileout = new TFile("signal.root","RECREATE");
  fileout->cd();
  TTree *tsig = new TTree("T","event");
  tsig->Branch("Nevent", &Neventsig, "Nevent/I"); // event number
  tsig->Branch("IDEvent", &IDEventsig, "IDEvent/I"); // ID of the event
  char strtmp[50];
  sprintf(strtmp,"FineTimeStamp[%d]/I",Nfeu);
  tsig->Branch("FineTimeStamp", FineTimeStampsig, strtmp); // time sample number
  sprintf(strtmp,"TsampleNum[%d]/I",Nsample);
  tsig->Branch("TsampleNum", TsampleNumsig, strtmp); // time sample number
  sprintf(strtmp,"StripNum[%d]/I",Nstrip);
  tsig->Branch("StripNum", StripNumsig, strtmp); // time sample number
  char leefStripAmpl[100];
  sprintf(leefStripAmpl,"StripAmpl_corrped[%d][%d][%d]/F",Ndet,Nstrip,Nsample);
  TBranch *newBranch = tsig->Branch("StripAmpl_corrped", StripAmpl_corrped,leefStripAmpl);
  
  Long64_t nentries = T->GetEntries();
  Long64_t nbytes = 0;
  set<double> sStripAmplSort1, sStripAmplSort2;  // to separate thin and large strips
  int stripmin=0, stripmax=64, stripmid=32;
  bool splittedStrips = false;

  for (Long64_t i=0; i<nentries; i++) {
    nbytes += T->GetEntry(i);
    //copy form intree to outtree
    Neventsig=Nevent;
    IDEventsig=IDEvent;
    for(int is = 0;is<Nsample;is++) TsampleNumsig[is]=TsampleNum[is];
    for(int is = 0;is<Nstrip;is++) StripNumsig[is]=StripNum[is];
    for(int is = 0;is<Nfeu;is++) FineTimeStampsig[is]=FineTimeStamp[is];
    for (register int idet=0; idet<Ndet; idet++) {
      for (register int isample=0; isample<Nsample; isample++) {
	sStripAmplSort1.clear();
	sStripAmplSort2.clear();
	if (splittedStrips) {
	  for (register int istrip=stripmin; istrip<stripmid; istrip++)
	    sStripAmplSort1.insert(StripAmpl[idet][istrip][isample] - Pedestal[idet][istrip]);
	  for (register int istrip=stripmid; istrip<stripmax; istrip++)
	    sStripAmplSort2.insert(StripAmpl[idet][istrip][isample] - Pedestal[idet][istrip]);
	  } else {
	  for (register int istrip=stripmin; istrip<stripmax; istrip++)
	    sStripAmplSort1.insert(StripAmpl[idet][istrip][isample] - Pedestal[idet][istrip]);
            }
	register int icount = sStripAmplSort1.size()/2;
	register double medAmpl1 = 0, medAmpl2 = 0;
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
	  for (register int istrip=stripmin; istrip<stripmid; istrip++)  {
	    StripAmpl_corrped[idet][istrip][isample] = StripAmpl[idet][istrip][isample] - Pedestal[idet][istrip] - medAmpl1;
	  }
	  for (register int istrip=stripmid; istrip<stripmax; istrip++)  {
	    StripAmpl_corrped[idet][istrip][isample] = StripAmpl[idet][istrip][isample] - Pedestal[idet][istrip] - medAmpl2;
	  }
	} else {
	  for (register int istrip=stripmin; istrip<stripmax; istrip++)  {//do the actual suppression
	    double val = StripAmpl[idet][istrip][isample] - Pedestal[idet][istrip] - medAmpl1;
	    if(nsigma>0)
	      {
		//cout << val << "  " << nsigma*RMS[idet][istrip] << endl;
		if(val>nsigma*RMS[idet][istrip])
		  StripAmpl_corrped[idet][istrip][isample] = val;
		else 
		  StripAmpl_corrped[idet][istrip][isample] = 0;
	      }
	  }
	}//else stripsplit
      }
    }
    tsig->Fill();
  }
  tsig->Write("", TObject::kOverwrite);
  return 1;
}

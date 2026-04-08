#include "Riostream.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TFile.h"
#include "TSystem.h"

//const int Ndet=0;


int RMSPedestauxNZS(string filename) {
  //   gROOT->Reset();
  TFile f(filename.c_str(), "");
  TTree *T = (TTree*)f.Get("T");
  //TTree *T = (TTree*)gROOT->FindObject("T");


  //get data size
  int Ndet=0;
  int Nsample=0;
  int Nstrip=0;
  const char * branchname =  T->GetBranch("StripAmpl")->GetTitle();  
  sscanf(branchname,"StripAmpl[%d][%d][%d]/F",&Ndet,&Nstrip,&Nsample);
  cout << "Data have " << Ndet << " DREAM, " << Nstrip << " strips, " 
       << Nsample << " samples." << endl;

  TH2F* h2_bg[Ndet];
 
  //const int Nstrip=64;
  const Int_t LowerCut=0; 
  const Int_t HigherCut=32;
  //const Int_t Nsample=HigherCut-LowerCut;
  const Int_t nsamplemin=0;
  const Int_t nsamplemax=4;
  const int NEvtMax=1000;

   //Declaration of leaves types  
  // General variables
  //   int Nevent; 		 // event number
  //   int TsampleNum[Nsample]; // index of time sample
  // variables of the output root tree for the DR2
  Int_t StripAmpl[Ndet][Nstrip][Nsample];
  
  // Set branch addresses.
  //   T->SetBranchAddress("Nevent",&Nevent);
  //   T->SetBranchAddress("TsampleNum",TsampleNum);
  T->SetBranchAddress("StripAmpl",StripAmpl);
  
  double Ymin = -300;
  double Ymax = 700;
  char strtmp[200];
  for (register int ii=0; ii<Ndet; ii++) {
    sprintf(strtmp, "h2_bg_%2d", ii);
    h2_bg[ii] = new TH2F(strtmp,"",Nstrip,0,Nstrip,500,Ymin,Ymax);
    h2_bg[ii]->SetStats(0);
  }

  // For all detectors
  Long64_t nbytes = 0;
  for (register int i=1; i<NEvtMax; i++) {
    nbytes += T->GetEntry(i);
    for (register int idet=0; idet<Ndet;  idet++) {
      for (register int istrip=0; istrip<Nstrip; istrip++) {
	for (register int isample=nsamplemin; isample<nsamplemax; isample++) {
	  h2_bg[idet]->Fill(istrip,StripAmpl[idet][istrip][isample]);
	  //cout << " " << StripAmpl[idet][istrip][isample] << flush;
	}
      }
    }
  }
  
  FILE * RMSfile;
  RMSfile = fopen ("RMSPed_Raw.dat","w");
  for (register int idet=0; idet<Ndet; idet++) {
    h2_bg[idet]->FitSlicesY();
    sprintf(strtmp, "h2_bg_%2d_2", idet);
    TH1F *h2_2 = (TH1F*)gDirectory->Get(strtmp);
    h2_2->SetStats(0);
    h2_2->SetMaximum(Nstrip);
    h2_2->SetMinimum(0);
    for (register int i=1; i<Nstrip+1; i++) {
      fprintf (RMSfile, "%d %d %5.2f\n",idet,i-1,h2_2->GetBinContent(i));
    }
    h2_bg[idet]->Draw();
  }

  fclose (RMSfile);

  return 1;
}

//root -l 'Pedestal.C++("output.root")'

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

const int NsampleRead=4;  // first samples without signals

int Pedestal(string filename) {
  gROOT->Reset();
  TFile f(filename.c_str(), "r");
  TTree *T = (TTree*)f.Get("T");

  //get data size
  int Ndet=0;
  int Nsample=0;
  int Nstrip=0;
  const char * branchname =  T->GetBranch("StripAmpl")->GetTitle();  
  sscanf(branchname,"StripAmpl[%d][%d][%d]/F",&Ndet,&Nstrip,&Nsample);
  cout << "Data have " << Ndet << " DREAM, " << Nstrip << " strips, " 
       << Nsample << " samples." << endl;
  
  //Declaration of leaves types
  //   Int_t           Nevent;
  Int_t           StripAmpl[Ndet][Nstrip][Nsample];
  //   Int_t           TsampleNum[Nsample];
  
  // Set branch addresses.
  //   T->SetBranchAddress("Nevent",&Nevent);
  T->SetBranchAddress("StripAmpl",StripAmpl);
  //   T->SetBranchAddress("TsampleNum",TsampleNum);
  
  Long64_t nentries = T->GetEntries();
  Long64_t nbytes = 0;
  set<double> sStripAmplSort;
  double Pedestal[Ndet][Nstrip];
  register int idet=0, istrip=0;
  for (idet=0; idet<Ndet; idet++) for(istrip=0; istrip<Nstrip; istrip++) Pedestal[idet][istrip]=0; // initialization

  double Nevt=1000;
  for (Long64_t i=1; i<Nevt+1;i++) {
    nbytes += T->GetEntry(i);

    for (idet=0; idet<Ndet; idet++) {
        for (istrip=0; istrip<Nstrip; istrip++) {
	    sStripAmplSort.clear();
            for (register int isample=0; isample<NsampleRead; isample++) sStripAmplSort.insert(StripAmpl[idet][istrip][isample]);
            register int icount = sStripAmplSort.size()/2;
            register double valAmpl = 0;
            for (set<double>::const_iterator isci = sStripAmplSort.begin(); isci != sStripAmplSort.end(); isci++) {
              if (icount == 0) {
                valAmpl = *isci;
                break;
              }
              icount--;
            }

	    Pedestal[idet][istrip] += valAmpl;
	  }
      }
  }
  for (idet=0; idet<Ndet; idet++) for (istrip=0; istrip<Nstrip; istrip++)  Pedestal[idet][istrip] /= Nevt;
  //for(int idet=0;idet<Ndet;idet++) for(int istrip=0;istrip<Nstrip;istrip++) cout << idet << " " << istrip << " " << Pedestal[idet][istrip] << endl;
  FILE * pedfile;
  pedfile = fopen ("Pedestal.dat","w");
  for (idet=0;idet<Ndet;idet++) for (istrip=0;istrip<Nstrip;istrip++) fprintf (pedfile, "%d %d %-5.2f\n",idet,istrip,Pedestal[idet][istrip]);
  fclose (pedfile);


  return 1;
}

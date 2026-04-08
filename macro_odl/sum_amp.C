#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TMath.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1.h"
#include "TProfile.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "TPaveText.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TCrown.h"
#include "TBox.h"
#include "TPolyLine.h"
#include "TText.h"
#include "TLine.h"
#include "TArc.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TColor.h"
#include "TStyle.h"

#include <iomanip>
#include <sstream>


using namespace std;


void sum_amp() {

   // --- Ouvrir le fichier ---
   TFile *f = new TFile("outtree_clust.root");
   if (!f || f->IsZombie()) {
      cout << "Impossible d'ouvrir outtree.root" << endl;
      return;
   }

   // --- Charger le TTree ---
   TTree *t = (TTree*)f->Get("Tout");
   if (!t) {
      cout << "Impossible de trouver le TTree 'Tout'" << endl;
      return;
   }

     // --- Branches ---
 
  
   vector<double> *DetAmp = 0;

  

   t->SetBranchAddress("DetAmp", &DetAmp);  
  
  vector<double> amp;
  Long64_t N = t->GetEntries();

  for (Long64_t iev = 0; iev < N; iev++) {
    double sum = 0;
    t->GetEntry(iev);  
    int nhits = DetAmp->size();
    for (int h = 0; h < nhits; h++) {
        sum += DetAmp->at(h);
    }
    amp.push_back(sum);
   }


    TH1F *h = new TH1F("h", "amp sum", 300, 0, 100000);

    for(int v = 0; v < amp.size(); v++){
        h->Fill(amp[v]);
    }
    
    TCanvas *c = new TCanvas("c","c",800,600);
    h->Draw();


}
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
#include "TGraph.h"
#include "TF1.h"
#include "TH1.h"
#include "TProfile.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

using namespace std;

void alignement(){

    
    TFile *f = new TFile("/data/cosmic_data/2025/W43/workdir_CosTb_P2_mapping_M400V_v2/outtree.root");
    if (!f || f->IsZombie()) {
        cout << "Impossible d'ouvrir outtree.root" << endl;
        return;
    }

    TTree *t = (TTree*)f->Get("Tout");
    if (!t) {
        cout << "Impossible de trouver le TTree 'Tout'" << endl;
        return;
    }

    // --- Branches ---
    double Chi2X; 
    double RayX, RayY;           
    vector<double> *DetX = 0;    
    vector<double> *DetY = 0;
    vector<double> *DetAmp = 0;

    t->SetBranchAddress("RayX", &RayX);
    t->SetBranchAddress("RayY", &RayY);
    t->SetBranchAddress("DetX", &DetX);
    t->SetBranchAddress("DetY", &DetY);

    t->SetBranchAddress("DetAmp", &DetAmp);  
    t->SetBranchAddress("Chi2X", &Chi2X);


  double bestChi2 = 1e30;
double bestTheta = 0, bestTx = 0, bestTy = 0;

Long64_t nEntries = t->GetEntries();

for (double theta_deg = -73; theta_deg <= -71; theta_deg += 0.1) {

double theta = theta_deg * TMath::DegToRad();

  double c = std::cos(theta);
  double s = std::sin(theta);

  double tx = 0.0;
  double ty = 0.0;
  int nUsed = 0;

  
  for (Long64_t ie = 0; ie < nEntries; ++ie) {
    t->GetEntry(ie);
if (!(Chi2X  < 0.1))  continue;

    if (!DetX || DetX->empty()) continue;


    int nPadsUsed = 0;
double xHit = 0.0, yHit = 0.0;

for (size_t i = 0; i < DetX->size(); ++i) {
    if (DetAmp->at(i) <= 500) continue;

    xHit += DetX->at(i);
    yHit += DetY->at(i);
    nPadsUsed++;
}

if (nPadsUsed == 0) continue;

xHit /= nPadsUsed;
yHit /= nPadsUsed;

    tx += xHit - (c*RayX - s*RayY);
    ty += yHit - (s*RayX + c*RayY);
    nUsed++;
  }


  tx /= nUsed;
  ty /= nUsed;

  
  double chi2 = 0.0;

  for (Long64_t ie = 0; ie < nEntries; ++ie) {
    t->GetEntry(ie);
if (!(Chi2X  < 0.1))  continue;

    if (!DetX || DetX->empty()) continue;

   int nPadsUsed = 0;
double xHit = 0.0, yHit = 0.0;

for (size_t i = 0; i < DetX->size(); ++i) {
    if (DetAmp->at(i) <= 500) continue;

    xHit += DetX->at(i);
    yHit += DetY->at(i);
    nPadsUsed++;
}

if (nPadsUsed == 0) continue;

xHit /= nPadsUsed;
yHit /= nPadsUsed;

    double xPred = c*RayX - s*RayY + tx;
    double yPred = s*RayX + c*RayY + ty;

    double dx = xHit - xPred;
    double dy = yHit - yPred;

    chi2 += dx*dx + dy*dy;
  }


  if (chi2 < bestChi2) {
    bestChi2  = chi2;
    bestTheta = theta_deg;
    bestTx    = tx;
    bestTy    = ty;
  }
}

std::cout << "Best alignment:\n";
std::cout << "  theta = " << bestTheta <<
" deg\n";
std::cout << "  tx = " << bestTx << "\n";
std::cout << "  ty = " << bestTy << "\n";

}

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

void translation() {

    // --- Ouvrir le fichier ---
    TFile *f = new TFile("outtree.root");
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
    double Chi2X; 
    double RayX, RayY;           // simples doubles
    vector<double> *DetX = 0;    // vecteurs
    vector<double> *DetY = 0;
    vector<double> *DetAmp = 0;

    t->SetBranchAddress("RayX", &RayX);
    t->SetBranchAddress("RayY", &RayY);
    t->SetBranchAddress("DetX", &DetX);
    t->SetBranchAddress("DetY", &DetY);

    t->SetBranchAddress("DetAmp", &DetAmp);  
    t->SetBranchAddress("Chi2X", &Chi2X);

    double sum_xr = 0, sum_yr = 0;
double sum_xm = 0, sum_ym = 0;
long count = 0;
double theta_deg = 0;
double theta = theta_deg * TMath::DegToRad();
Long64_t N = t->GetEntries();
for (Long64_t i=0; i<N; i++) {
    t->GetEntry(i);
    if (!(Chi2X  < 0.1))  continue;
    int nh = DetX->size();
    for (int h=0; h<nh; h++) {
        if (DetAmp->at(h) <= 500) continue;
        double xm = DetX->at(h);
        double ym = DetY->at(h);

        double xm_rot = xm*cos(theta) - ym*sin(theta);
        double ym_rot = xm*sin(theta) + ym*cos(theta);

        sum_xm += xm_rot;
        sum_ym += ym_rot;

        sum_xr += RayX;
        sum_yr += RayY;

        count++;
    }
}

double mean_xr = sum_xr / count;
double mean_yr = sum_yr / count;
double mean_xm = sum_xm / count;
double mean_ym = sum_ym / count;

double tx = mean_xr - mean_xm;
double ty = mean_yr - mean_ym;

cout << "Translation optimale :" << endl;
cout << "  t_x = " << tx << endl;
cout << "  t_y = " << ty << endl;

}

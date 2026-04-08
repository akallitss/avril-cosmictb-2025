#include <vector>
#include <string>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH2D.h"

using namespace std;


void efficacite_no_mapping() {

    // --- Ouvrir le fichier ---
    TFile *f = new TFile("/data/cosmic_data/2026/W2/workdir_CosTb_det_480_800/outtree.root");
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

    vector<double>* amp = 0;
    double rayx = 0;
    double rayy = 0;
    double chi2x = 0;
    double chi2y = 0;
    Long64_t nseen = 0;
    Long64_t ngood = 0;
    Long64_t nchi = 0;
    Long64_t nentries = t->GetEntries();
    
    t->SetBranchAddress("DetAmp", &amp);
    t->SetBranchAddress("RayX", &rayx);
    t->SetBranchAddress("RayY", &rayy);
    t->SetBranchAddress("Chi2X", &chi2x); 
    t->SetBranchAddress("Chi2Y", &chi2y);



TH2D *hRayPassed = new TH2D(
    "hRayPassed",";RayX;RayY",
    100, -300, 300,
    100, -300, 300
);

TH2D *hRaySeen = new TH2D(
    "hRaySeen",";RayX;RayY",
    100, -300, 300,
    100, -300, 300
);


for (Long64_t k = 0; k < nentries; k++) {
        t->GetEntry(k);
        if (!amp){ cout<<"pas de Amp"<< endl;continue;} 
        // if (chi2x > 0.5) continue;
        // if (chi2y > 0.5) continue;
        nchi ++;
        // if (rayx < -40 || rayx > 40) continue;
        // if (rayy < 10  || rayy > 90) continue;
        
        ngood++; 
        hRayPassed->Fill(rayx, rayy);
        bool seen = false;
        for (size_t j = 0; j < amp->size(); j++) {
            if (amp->at(j) > 50) {
                seen = true;
                break;
            }
        }
        if (seen) {nseen++ ;hRaySeen->Fill(rayx, rayy);}
    }
    double goodchi = double(nchi)/double(nentries) ;
     cout << goodchi <<  endl;
    double efficiency = (nentries > 0) ? double(nseen)/double(ngood) : 0;
     cout << "efficiency = " << efficiency <<  endl;
     cout << hRaySeen->GetEntries() << endl;


TH2D *hRayEff = (TH2D*)hRaySeen->Clone("hRayEff");
hRayEff->Divide(hRayPassed);
hRayEff->SetTitle("efficiency (DetAmp > 50)");

TCanvas *cEff = new TCanvas("cEff","Ray efficiency",1000,1000);
hRayEff->SetMinimum(0);
hRayEff->SetMaximum(1);
hRayEff->Draw("COLZ");

}
#include <vector>
#include <string>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"

void eff_vs_volt() {

std::vector<std::string> filenames;  
filenames.push_back("/data/cosmic_data/2025/W43/workdir_CosTb_p2_mesh_420V_datrun_251025_00H05/outtree_420V.root");
filenames.push_back("/data/cosmic_data/2025/W43/workdir_CosTb_p2_mesh_415V_datrun_251025_04H10/outtree_415V.root");
filenames.push_back("/data/cosmic_data/2025/W43/workdir_CosTb_p2_mesh_410V_datrun_251025_08H15/outtree_410V.root");
filenames.push_back("/data/cosmic_data/2025/W43/workdir_CosTb_p2_mesh_405V_datrun_251025_12H20/outtree_405V.root");                  
filenames.push_back("/data/cosmic_data/2025/W43/workdir_CosTb_p2_mesh_400V_datrun_251025_16H24/outtree_400V.root");
filenames.push_back("/data/cosmic_data/2025/W43/workdir_CosTb_p2_mesh_395V_datrun_251025_20H28/outtree_395V.root");
filenames.push_back("/data/cosmic_data/2025/W43/workdir_CosTb_p2_mesh_390V_datrun_251026_00H33/outtree_390V.root");
filenames.push_back("/data/cosmic_data/2025/W43/workdir_CosTb_p2_mesh_385V_datrun_251026_03H38/outtree_385V.root");
filenames.push_back("/data/cosmic_data/2025/W43/workdir_CosTb_p2_mesh_380V_datrun_251026_07H42/outtree_380V.root");
filenames.push_back("/data/cosmic_data/2025/W43/workdir_CosTb_p2_mesh_375V_datrun_251026_11H46/outtree_375V.root");
filenames.push_back("/data/cosmic_data/2025/W43/workdir_CosTb_p2_mesh_370V_datrun_251026_15H51/outtree_370V.root");

std::vector<double> volt;
volt.push_back(420);
volt.push_back(415);
volt.push_back(410);
volt.push_back(405);
volt.push_back(400);
volt.push_back(395);
volt.push_back(390);
volt.push_back(385);
volt.push_back(380);
volt.push_back(375);
volt.push_back(370);

std::vector<double> eff;
double amp_min = 100;

for (size_t i = 0; i < filenames.size(); i++) {
    TFile* f = TFile::Open(filenames[i].c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "Impossible d'ouvrir " << filenames[i] << std::endl;
        eff.push_back(0);
        continue;
    }

    TTree* t = (TTree*)f->Get("Tout");
    if (!t) {
        std::cerr << "Arbre 'Tout' introuvable dans " << filenames[i] << std::endl;
        eff.push_back(0);
        f->Close();
        continue;
    }

    std::vector<double>* amp = 0;
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

    for (Long64_t k = 0; k < nentries; k++) {
        t->GetEntry(k);
        if (!amp){std::cout<<"pas de Amp"<<std::endl;continue;} 
        if (chi2x > 100) continue;
        if (chi2y > 100) continue;
        nchi ++;
        // if (rayx < 30 || rayx > 130) continue;
        // if (rayy < 0  || rayy > 100) continue;
        
        ngood++; 
        bool seen = false;
        for (size_t j = 0; j < amp->size(); j++) {
            if (amp->at(j) > amp_min) {
                seen = true;
                break;
            }
        }
        if (seen) nseen++;
    }
    double goodchi = double(nchi)/double(nentries) ;
    std::cout << goodchi << std::endl;
    double efficiency = (nentries > 0) ? double(nseen)/double(ngood) : 0;
    std::cout << "efficiency = " << efficiency << std::endl;
    eff.push_back(efficiency);

    f->Close();
}

TCanvas* c1 = new TCanvas("c1","Efficiency vs tension (chi2 selection)",800,600);
TGraph* g = new TGraph(volt.size(), &volt[0], &eff[0]);
g->SetTitle("Efficiency vs tension (chi2 selection);Voltage [V];Efficiency");
g->SetMarkerStyle(20);
g->Draw("APL");
c1->Update();
}

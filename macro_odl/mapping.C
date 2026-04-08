#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TText.h"

std::vector<double> clean_vector(const std::vector<double> &v) {
    if (v.size() < 3) return v; 
    std::vector<double> tmp = v;
    std::sort(tmp.begin(), tmp.end());
    double median = tmp[tmp.size() / 2];

    
    std::vector<double> deviations;
    deviations.reserve(v.size());
    for (size_t i = 0; i < v.size(); ++i)
        deviations.push_back(fabs(v[i] - median));
    std::sort(deviations.begin(), deviations.end());
    double mad = deviations[deviations.size() / 2];
    if (mad == 0) return v; 

    std::vector<double> vclean;
    vclean.reserve(v.size());
    for (size_t i = 0; i < v.size(); ++i)
        if (fabs(v[i] - median) < 3 * mad)
            vclean.push_back(v[i]);

    return vclean;
}


void mapping() {

    std::string filename = "outtree.root";
    TFile *f = TFile::Open(filename.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "Erreur: impossible d’ouvrir " << filename << std::endl;
        return;
    }

    TTree *t = (TTree*)f->Get("Tout");
    if (!t) {
        std::cerr << "Erreur: arbre 'Tout' introuvable dans " << filename << std::endl;
        return;
    }

    double rayx = 0;
    double rayy = 0;
    std::vector<int> *detch = 0;
    std::vector<double> *amp = 0;

    t->SetBranchAddress("RayX", &rayx);
    t->SetBranchAddress("RayY", &rayy);
    t->SetBranchAddress("DetCh", &detch);
    t->SetBranchAddress("DetAmp", &amp);

    std::map<int, std::vector<double> > xvals, yvals;
    Long64_t nentries = t->GetEntries();

    for (Long64_t i = 0; i < nentries; i++) {
        t->GetEntry(i);
        if (!amp || !detch) continue;
        if (amp->size() != detch->size()) continue;

        for (size_t j = 0; j < amp->size(); j++) {
            int ch = detch->at(j);
            if (ch < 0) continue;
            if (amp->at(j) > 100) {
                xvals[ch].push_back(rayx);
                yvals[ch].push_back(rayy);
            }
        }
    }

    // --- Calcul barycentres propres avec rejet d’outliers
    std::vector<double> bx, by;
    std::vector<int> channels;

    for (std::map<int, std::vector<double> >::iterator it = xvals.begin(); it != xvals.end(); ++it) {
        int ch = it->first;
        std::vector<double> &vx = xvals[ch];
        std::vector<double> &vy = yvals[ch];

        if (vx.size() < 5) continue; 


        std::vector<double> cx = clean_vector(vx);
        std::vector<double> cy = clean_vector(vy);
        if (cx.empty() || cy.empty()) continue;

        double meanx = std::accumulate(cx.begin(), cx.end(), 0.0) / cx.size();
        double meany = std::accumulate(cy.begin(), cy.end(), 0.0) / cy.size();

        bx.push_back(meanx);
        by.push_back(meany);
        channels.push_back(ch);
    }

    // --- Affichage graphique
    TCanvas *c = new TCanvas("c", "Barycentres ", 800, 700);
    TGraph *g = new TGraph(bx.size(), &bx[0], &by[0]);
    g->SetTitle("Barycentres par channel ;RayX;RayY");
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1);
    g->Draw("AP");

    // Annoter les numéros de channels
    for (size_t i = 0; i < bx.size(); i++) {
        TText *label = new TText(bx[i], by[i], Form("%d", channels[i]));
        label->SetTextSize(0.015);
        label->Draw("same");
    }
    
    f->Close();
}

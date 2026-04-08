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

void sigma_rotation() {

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

    vector<double> angles;
    vector<double> sigmas;
    vector<double> slopes;
vector<double> sigmasy;
    vector<double> slopesy;

for (double theta_deg = 73; theta_deg <= 73.5; theta_deg += 0.1) {

    double theta = theta_deg * TMath::DegToRad();

        TH1F hres("hres",
            Form("Residual X (theta = %.1f deg); RayX - DetX_rot; Counts",
            theta_deg),
            200, -1000, 600);
        TH1F hresy("hresy",
            Form("Residual X (theta = %.1f deg); RayX - DetX_rot; Counts",
            theta_deg),
            200, -1000, 600);

        TH2F h2("h2",
            Form("Residual vs DetX (theta = %.1f deg); DetX; RayX - DetX_rot",
            theta_deg),
            200, 50, 600,
            200, -1000, 600);
        TH2F h2y("h2y",
            Form("Residual vs DetX (theta = %.1f deg); DetX; RayX - DetX_rot",
            theta_deg),
            200, 0, 550,
            200, -1000, 600);
hres.SetDirectory(nullptr);
hresy.SetDirectory(nullptr);
h2.SetDirectory(nullptr);
h2y.SetDirectory(nullptr);

        

        Long64_t N = t->GetEntries();

        for (Long64_t i = 0; i < N; i++) {

            t->GetEntry(i);
            
            if (!(Chi2X  < 0.1))  continue;
            int nhits = DetX->size();

            for (int h = 0; h < nhits; h++) {
                if (DetAmp->at(h) <= 500) continue;
                double xm = DetX->at(h);
                double ym = DetY->at(h);

                // --- Rotation complète 2D ---
                double xm_rot = xm*cos(theta) - ym*sin(theta);
                double ym_rot = xm*sin(theta) + ym*cos(theta);

                double dx = RayX - xm_rot;
                double dy = RayY - ym_rot;

                hres.Fill(dx);
                h2.Fill(xm, dx);
                hresy.Fill(dy);
                h2y.Fill(ym, dy);
                
            }
        }


    // RMS du résidu
TF1 gfit("gfit","gaus",300,500);
hres.Fit(&gfit,"Q");
double sigma = gfit.GetParameter(2);
  

    TF1 f("f","pol1");
    TProfile *p = h2.ProfileX();
    p->Fit(&f,"Q");  // Q = silencieux

    double slope = f.GetParameter(1);   // pente p1

   
    slopes.push_back(slope);
    angles.push_back(theta_deg);
    sigmas.push_back(sigma);



TF1 gfity("gfity","gaus",300,500);
hresy.Fit(&gfity,"Q");
double sigmay = gfity.GetParameter(2);


    TF1 fy("fy","pol1");
    TProfile *py = h2y.ProfileX();
    py->Fit(&fy,"Q");  // Q = silencieux

    double slopey = fy.GetParameter(1);   // pente p1

    slopesy.push_back(slopey);
    sigmasy.push_back(sigmay);
    
}
int n = angles.size();
TGraph *g = new TGraph(n, &angles[0], &sigmas[0]);
g->SetTitle("Sigma(residual X) vs rotation angle; angle (deg); sigma");
g->SetMarkerStyle(20);

double best_angle = angles[TMath::LocMin(sigmas.size(), &sigmas[0])];
cout << "Best angle = " << best_angle << " deg" << endl;

TGraph *g2 = new TGraph(n, &angles[0], &slopes[0]);
g2->SetTitle("Slope (residual X vs DetX) vs rotation angle; angle (deg); slope");
g2->SetMarkerStyle(20);

TGraph *g3 = new TGraph(n, &angles[0], &sigmasy[0]);
g3->SetTitle("Sigma(residual Y) vs rotation angle; angle (deg); sigma");
g3->SetMarkerStyle(20);

double best_angle_y = angles[TMath::LocMin(sigmasy.size(), &sigmasy[0])];
cout << "Best angle y = " << best_angle_y << " deg" << endl;

TGraph *g4 = new TGraph(n, &angles[0], &slopesy[0]);
g4->SetTitle("Slope (residual Y vs DetY) vs rotation angle; angle (deg); slope");
g4->SetMarkerStyle(20);

TCanvas *c = new TCanvas("c","Results",1200,500);
c->Divide(2,2);

c->cd(1);
g->Draw("APL");

c->cd(2);
g2->Draw("APL");

c->cd(3);
g3->Draw("APL");

c->cd(4);
g4->Draw("APL");

// Trouver l’angle où |slope| est minimum
double best_angle_slope = angles[0];
double best_val = fabs(slopes[0]);

for (int i=1; i<n; i++) {
    double val = fabs(slopes[i]);
    if (val < best_val) {
        best_val = val;
        best_angle_slope = angles[i];
    }
}

cout << "Best angle = " << best_angle_slope << " deg (slope closest to 0)" << endl;   


double best_angle_slope_y = angles[0];
double best_val_y = fabs(slopesy[0]);

for (int i=1; i<n; i++) {
    double valy = fabs(slopesy[i]);
    if (valy < best_val_y) {
        best_val_y = valy;
        best_angle_slope_y = angles[i];
    }
}

cout << "Best angle y = " << best_angle_slope_y << " deg (slope closest to 0)" << endl;   
}

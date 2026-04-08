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
#include "TPaveText.h"

using namespace std;

void angle_loop(double min_angle, double max_angle, double pas, double histSize, double &best_angle, double &Tx, double &Ty, bool stop, TTree* t, double &RayX, double &RayY, double &Chi2X, double &Chi2Y, vector<double>* &DetX, vector<double>* &DetY, vector<double>* &DetAmp){
// if bool stop is true will have to press enter between each plots
    // --- Branches ---
    // double Chi2X, Chi2Y; 
    // double RayX, RayY;           
    // vector<double> *DetX = 0;   
    // vector<double> *DetY = 0;
    // vector<double> *DetAmp = 0;

    // t->SetBranchAddress("RayX", &RayX);
    // t->SetBranchAddress("RayY", &RayY);
    // t->SetBranchAddress("DetX", &DetX);
    // t->SetBranchAddress("DetY", &DetY);

    // t->SetBranchAddress("DetAmp", &DetAmp);  
    // t->SetBranchAddress("Chi2X", &Chi2X);
    // t->SetBranchAddress("Chi2Y", &Chi2Y);

    //TString cuts = "DetAmp > 100 && Chi2X < 0.1 && Chi2Y < 0.1";
  
	vector<double> sigma;
	vector<double> angles;
	vector<double> mean_x;
	vector<double> mean_y;

    // --- Boucle sur les angles ---
    for (double theta_deg = min_angle; theta_deg <= max_angle; theta_deg += pas) {

        double theta = theta_deg * TMath::DegToRad();
	    angles.push_back(theta_deg);

        TH1F *hres = new TH1F("hres",
            Form("Residual X (theta = %.1f deg); RayX - DetX_rot; Counts",
            theta_deg),
            200, Tx-histSize, Tx+histSize);
        
        TH1F *hresy = new TH1F("hresy",
            Form("Residual Y (theta = %.1f deg); RayY - DetY_rot; Counts",
            theta_deg),
            200, Ty-histSize,Ty+histSize);

        TH2F *h2 = new TH2F("h2",
            Form("Residual vs DetX (theta = %.1f deg); DetX; RayX - DetX_rot",
            theta_deg),
            200, 50, 600,
            200, -600, 600);

        TH2F *h2y = new TH2F("h2y",
            Form("Residual vs DetY (theta = %.1f deg); DetY; RayY - DetY_rot",
            theta_deg),
            200, 0, 550,
            200, -600, 600);

        Long64_t N = t->GetEntries();

        for (Long64_t i = 0; i < N; i++) {

            t->GetEntry(i);
              if (!DetX || !DetY || !DetAmp){cout << " detx vide "<< endl;continue; } 
            if (!(Chi2X  < 1))  continue;
	    if (!(Chi2Y  < 1))  continue;
            int nhits = DetX->size();

            for (int h = 0; h < nhits; h++) {
                if (DetAmp->at(h) <= 100) continue;
                double xm = DetX->at(h);
                double ym = DetY->at(h);

                // --- Rotation complète 2D ---
                double xm_rot = xm*cos(theta) - ym*sin(theta);
                double ym_rot = xm*sin(theta) + ym*cos(theta);

                // résidu uniquement en X
                double dx = RayX - xm_rot;
                double dy = RayY - ym_rot;

                hres->Fill(dx);
                hresy->Fill(dy);
                h2->Fill(xm, dx);
                h2y->Fill(ym, dy);
            }
        }
TCanvas *c = new TCanvas(Form("c_%d",(int)theta_deg),
                         Form("theta = %.1f deg", theta_deg),
                         1100, 900);

c->Divide(2,2);

// PAD 1 : hres + fit gauss
c->cd(1);
hres->Draw();
double mean = hres->GetMean();
mean_x.push_back(mean);
// Fit gaussien
TF1 *gfit = new TF1("gfit","gaus",mean-100,mean+100);
hres->Fit(gfit, "QR");

// paramètres du fit + RMS
double sigma_fit = gfit->GetParameter(2);

// Affichage du texte
TPaveText *pt = new TPaveText(0.60,0.60,0.880,0.8880,"NDC");
pt->AddText(Form("Sigma = %.2f", sigma_fit));

pt->SetFillColor(0);
pt->Draw();

// PAD 2 : 2D X
c->cd(2);
h2->Draw();

// PAD 3 : residual Y
c->cd(3);
hresy->Draw();

double meany = hresy->GetMean();
mean_y.push_back(meany);
// Fit gaussien
TF1 *gfity = new TF1("gfity","gaus",meany-100,meany+100);
hresy->Fit(gfity, "QR");

// paramètres du fit + RMS

double sigma_fity = gfity->GetParameter(2);

// Affichage du texte
TPaveText *pty = new TPaveText(0.60,0.60,0.880,0.8880,"NDC");
pty->AddText(Form("Sigma = %.2f", sigma_fity));

pty->SetFillColor(0);
pty->Draw();

// PAD 4 : 2D Y
c->cd(4);
h2y->Draw();

sigma.push_back(sigma_fit*sigma_fit + sigma_fity*sigma_fity);
c->Update();

        if(stop){
        cout << "Angle = " << theta_deg << " deg → appuie sur ENTREE" << endl;
        cin.get();}

        delete hres;
        delete h2;
        delete hresy;
        delete h2y;
    }
if (sigma.size() == 0) {
    cout << "ERREUR : sigma vide !" << endl;
    return;
}
best_angle = angles[TMath::LocMin(sigma.size(), &sigma[0])];
cout << "best angle = " << best_angle << endl;
//  double Tx = 0;
// double Ty= 0;
for (size_t i = 0; i < angles.size(); ++i) {
    if (angles[i] == best_angle) {
        Tx = mean_x[i];
	Ty= mean_y[i];

        break;
    }

}
}

void rotation() {

      // --- Ouvrir le fichier ---
    TFile *f = new TFile("outtree_final.root");
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
    double Chi2X, Chi2Y; 
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
    t->SetBranchAddress("Chi2Y", &Chi2Y);
    
double best_angle = 0; double Tx =0; double Ty=0;
angle_loop(50, 90, 5, 600, best_angle, Tx, Ty, false, t,RayX, RayY,Chi2X, Chi2Y,DetX, DetY, DetAmp);
angle_loop(best_angle - 5, best_angle + 5, 1, 200, best_angle, Tx, Ty, false, t,RayX, RayY,Chi2X, Chi2Y,DetX, DetY, DetAmp);
angle_loop(best_angle - 1, best_angle + 1, 0.1, 50, best_angle, Tx, Ty, false, t, RayX, RayY,Chi2X, Chi2Y,DetX, DetY, DetAmp);

cout << "tx et ty : " << Tx <<" "<<  Ty << endl;

TFile* fout = new TFile("outtree_aligned.root", "RECREATE");
cout << " file outtre_aligned created" << endl;

TTree* t_out = t->CloneTree(0);
cout << " ttree cloned" << endl;
  
double RayX_corr = 0.0;
double RayY_corr = 0.0;

TBranch* bRayXcorr = t_out->Branch("RayX_corr", &RayX_corr, "RayX_corr/D");
TBranch* bRayYcorr = t_out->Branch("RayY_corr", &RayY_corr, "RayY_corr/D");

double theta = best_angle * TMath::DegToRad();
double c = std::cos(theta);
double s = std::sin(theta);
cout << " tx et ty : " << Tx << " " << Ty << endl;
Long64_t nEntries = t->GetEntries();

for (Long64_t i = 0; i < nEntries; ++i) {

    t->GetEntry(i);

    // correction de RayX, RayY
   RayX_corr =  c * (RayX - Tx) + s * (RayY - Ty);
RayY_corr = -s * (RayX - Tx) + c * (RayY - Ty);


    t_out->Fill();
}
fout->cd();
t_out->Write();
fout->Close();

std::cout << "Fichier outtree_aligned.root écrit avec succès" << std::endl;
}
